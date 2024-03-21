#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <argp.h>
#include "seq_utils.h"
#include "mcmc.h"

#define MAX_MUT 4
// default mode is standard MCMC
Bool OPT = false;

char amino_acids_char[] = {'A' ,'C' ,'T' ,'E' ,'D' ,'F' ,'W' ,'I' ,'V' ,'L' ,'K' ,'M' ,'N' ,'Q' ,'S' ,'R' ,'Y' ,'H' ,'P' ,'G', '*', '-'};
char nuc_2_char[] = {'A','C','T','G'};

// Parse arguments -------------------------------------------------------------
const char *argp_program_version = "Version 1.0";
const char *argp_program_bug_address = "@polytechnique.edu";

/* Program documentation. */
static char doc[] = "Standard MCMC procedure";

/* A description of the arguments we accept. */
static char args_doc[] = "<seq> <dca>";

/* The options we understand. */
static struct argp_option options[] = {
                                       {"nb_steps",   'n', "INT", 0, "Nb of MC steps", 0 },
                                       {"temp",   't', "FLOAT", 0, "temperature KT", 0 },
                                       {"mode",   'm', "<MC|OPT>", 0, "Normal MCMC or annealing procedure", 0 },
                                       { 0 }
};

/* Used by main to communicate with parse_opt. */
struct arguments
{
  char* args[2];
  char* mode;
  float temp;
  int nb_steps;
};

/* Parse a single option. */
static error_t parse_opt (int key, char *arg, struct argp_state *state) {
  /* Get the input argument from argp_parse, which we
     know is a pointer to our arguments structure. */
  struct arguments *arguments = state->input;

  switch (key)
    {
    case 't':
      arguments->temp = atof(arg);
      break;
    case 'n':
      arguments->nb_steps = atoi(arg);
      break;
    case 'm':
      arguments->mode = arg;
      break;

    case ARGP_KEY_ARG:
      if (state->arg_num >= 2)
        /* Too many arguments. */
        argp_usage (state);

      arguments->args[state->arg_num] = arg;

      break;

    case ARGP_KEY_END:
      if (state->arg_num < 2)
        /* Not enough arguments. */
        argp_usage (state);
      break;

    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc };
// Parse arguments -------------------------------------------------------------

// Initialize the conf
mc_conf_p initialize(char* seq_str) {
  int i;
  int len = strlen(seq_str);
  mc_conf_p init_conf = (mc_conf_p)malloc(sizeof(mc_conf));
  Aminoacid* seq = (Aminoacid*)malloc(sizeof(Aminoacid)*len);
  init_conf->seq = seq;

  // convert the strings
  for (i = 0; i < len; i++)
    init_conf->seq[i] = convert_aa(seq_str[i]);

  init_conf->len = len;
  init_conf->count = 1;
  return init_conf;
}

// free configuration
void free_conf(mc_conf_p conf) {
  free(conf->seq);
  free(conf);
}

// Propose a random mutation in the dna sequence return the position to change
// in both X and Y the position in the dna seq, and the aa and nuc to change
mutation_mc_p* prop_mutation(mc_conf_p old_conf, mutation_mc_p* tmp_mutation) {
  int nb_mut = rand_int(1, MAX_MUT), i, pmut, j;
  Aminoacid aa;

  // reinitialize the mutations
  for (i = 0; i < MAX_MUT; i++)
    tmp_mutation[i] = NULL;

  for (i = 0; i < nb_mut; i++) {
    aa = (Aminoacid)rand_int(0, 18);
    pmut = rand_int(0, old_conf->len-i);
    for (j = 0; j < i; j++) {
      if (pmut >= tmp_mutation[j]->pos)
        pmut++;
    }
    if (aa >= old_conf->seq[pmut])
      aa++;
    // save the proposed mutation
    tmp_mutation[i] = (mutation_mc_p)malloc(sizeof(mutation_mc));
    tmp_mutation[i]->pos = pmut;
    tmp_mutation[i]->aa = aa;
  }
  return tmp_mutation;
}

// apply the propose mutation in the protein
void apply_mutation(mc_conf_p old_conf, mutation_mc_p* prop_mut) {
  int i;
  for (i = 0; i < MAX_MUT && prop_mut[i] != NULL; i++) {
    old_conf->seq[prop_mut[i]->pos] = prop_mut[i]->aa;
  }
  old_conf->count = 1;
}

void print_conf(mc_conf_p conf) {
  int i;
  printf("\\SEQ  %8d %8.3f ", conf->count, conf->nrj);
  for (i = 0; i < conf->len; i++)
    printf("%c", amino_acids_char[conf->seq[i]]);
  printf("\n");
}

// compute the delta nrj for a given mutation
double update_nrj(Aminoacid* seq, Aminoacid* tmp_seq, int len_seq, float**** current_score, mutation_mc_p* prop_mut) {
  int i, imut, cpos;
  Aminoacid aai, iaai, maai, imaai;

  double delta_nrj = 0.0;
  // re-Initialize to the current sequence
  for (i = 0; i < len_seq; i++)
    tmp_seq[i] = seq[i];

  // Change to the new sequence
  for (imut = 0; imut < MAX_MUT && prop_mut[imut] != NULL; imut++)
    tmp_seq[prop_mut[imut]->pos] = prop_mut[imut]->aa;

  for (imut = 0; imut < MAX_MUT && prop_mut[imut] != NULL; imut++) {
    cpos = prop_mut[imut]->pos;
    aai = seq[cpos];
    maai = tmp_seq[cpos];
    for (i = 0; i < len_seq; i++) {
      iaai = seq[i];
      imaai = tmp_seq[i];
      delta_nrj += current_score[cpos][i][maai][imaai] - current_score[cpos][i][aai][iaai];
    }
  }
  return delta_nrj;
}

// compute the total score nrj
double compute_nrj(Aminoacid* seq, int len_seq, float**** current_score) {
  double tot_nrj = 0.0;
  int i, j, ui=0, uj=0;
  for (i = 0; i < len_seq; i++) {
    if (seq[i] != GAP) {
      uj = ui;
      for (j = i; j < len_seq && seq[j] != GAP; j++) {
        if (seq[j] != GAP) {
          if ((uj - ui) < 1) {
            if (current_score[ui][uj][seq[i]][seq[j]] > 10000)
              printf("%d %d %d %d %c %c %f\n", i, j, ui, uj, amino_acids_char[seq[i]], amino_acids_char[seq[j]], current_score[ui][uj][seq[i]][seq[j]]);
            tot_nrj += current_score[ui][uj][seq[i]][seq[j]];
          }
          uj++;
        }
      }
      ui++;
    }
  }
  return tot_nrj;
}

// MCMC in DNA sequence space
// - For each step, perform a mutation in the DNA sequence, then compute the nrj
// - We start from the sequence optimized with diagonal terms only
void mcmc_run(int nb_steps, float temp, mc_conf_p old_conf, float**** scores) {
  mutation_mc_p* current_mut = (mutation_mc_p*)malloc(sizeof(mutation_mc)*MAX_MUT);
  Aminoacid* tmp_seq = (Aminoacid*)malloc(sizeof(Aminoacid)*old_conf->len);
  double beta = 1.0/temp;
  double nrj_x = 0.0, nrj_y = 0.0;
  int i, j;
  double delta_nrj = 0.0, rand_val;
  // from initial temp to final temp => cooling scheme
  double final_temp = 0.00001, init_temp = temp;
  double mult;
  // If standard MC, don't change the temperature
  if (OPT == true)
    mult = exp(log(final_temp/init_temp)/(float)nb_steps);
  else
    mult = 1.0;

  old_conf->nrj = compute_nrj(old_conf->seq, old_conf->len, scores);
  print_conf(old_conf);

  // main mcmc loop
  for (i = 0; i < nb_steps; i++) {
    // propose a mutation
    prop_mutation(old_conf, current_mut);

    // update the nrj
    delta_nrj = update_nrj(old_conf->seq, tmp_seq, old_conf->len, scores, current_mut);

    // metropolis scheme
    if (delta_nrj < 0.0 || exp(-(1.0/temp) * delta_nrj) >= rand_float(0, 1)) {
      old_conf->nrj += delta_nrj;

      // save the mutation in the actual configuration
      apply_mutation(old_conf, current_mut);
      print_conf(old_conf);
    } else {
      // update the count of the current conf
      old_conf->count++;
    }
    temp = temp * mult;
  }
  print_conf(old_conf);
  free(current_mut);
}

int main(int argc, char *argv[]) {
  struct arguments arguments;
  int* frames;
  int nb_steps;
  fasta_seq_p seq_el;
  char *seq_str, *dca_file;
  arguments.temp = 0.6;
  arguments.nb_steps = 1000;
  arguments.mode = "MC";
  argp_parse(&argp, argc, argv, 0, 0, &arguments);

  if (strcmp(arguments.mode, "OPT") == 0)
    OPT = true;

  seq_el = read_fasta(arguments.args[0]);
  seq_str = seq_el->seq;
  dca_file = arguments.args[1];

  mc_conf_p init_conf = initialize(seq_str);
  float**** scores = read_score(init_conf->len, dca_file, 1);

  mcmc_run(arguments.nb_steps, arguments.temp, init_conf, scores);

  free(seq_str);
  free_score(scores, init_conf->len);
  free_conf(init_conf);
  return 0;
}

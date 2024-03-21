#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <argp.h>
#include "seq_utils.h"
#include "design.h"

#define MAX_CHAR 3000

char amino_acids_char[] = {'A' ,'C' ,'T' ,'E' ,'D' ,'F' ,'W' ,'I' ,'V' ,'L' ,'K' ,'M' ,'N' ,'Q' ,'S' ,'R' ,'Y' ,'H' ,'P' ,'G', '*', '-'};
char nuc_2_char[] = {'A','C','T','G'};

// Parse arguments -------------------------------------------------------------
const char *argp_program_version = "Version 1.0";
const char *argp_program_bug_address = "@polytechnique.edu";

/* Program documentation. */
static char doc[] = "Overlapping genes pair optimization";

/* A description of the arguments we accept. */
static char args_doc[] = "<seq_x> <seq_y> <dna> <score_x> <score_y> <frame>";

/* The options we understand. */
static struct argp_option options[] = {
                                       {"phase", 'f', "<-2|-1|0|1|2>", 0, "Overlapping frame, negative frames are n2 and n1", 0},
                                       {"nb_steps",   'n', "INT", 0, "Nb of MC steps", 0 },
                                       {"temp",   't', "FLOAT", 0, "temperature KT", 0 },
                                       { 0 }
};

/* Used by main to communicate with parse_opt. */
struct arguments
{
  char* args[6];
  float temp;
  char *frame;
  int nb_steps;
};

/* Parse a single option. */
static error_t parse_opt (int key, char *arg, struct argp_state *state) {
  /* Get the input argument from argp_parse, which we
     know is a pointer to our arguments structure. */
  struct arguments *arguments = state->input;

  switch (key)
    {
    case 'f': arguments->frame = (char*)malloc(sizeof(char)*MAX_CHAR); strcpy(arguments->frame, arg); break;
    case 't': arguments->temp = atof(arg); break;
    case 'n': arguments->nb_steps = atoi(arg); break;

    case ARGP_KEY_ARG:
      if (state->arg_num >= 5)
        /* Too many arguments. */
        argp_usage (state);

      arguments->args[state->arg_num] = arg;

      break;

    case ARGP_KEY_END:
      if (state->arg_num < 5)
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
over_conf_p initialize(char* seq_x_str, char* seq_y_str, char* dna_str, int* frames) {
  int i;
  int len_x = strlen(seq_x_str), len_y = strlen(seq_y_str), len_dna = strlen(dna_str);
  over_conf_p init_conf = (over_conf_p)malloc(sizeof(over_conf));
  Aminoacid* seq_x = (Aminoacid*)malloc(sizeof(Aminoacid)*len_x);
  Aminoacid* seq_y = (Aminoacid*)malloc(sizeof(Aminoacid)*len_y);
  Nucleotide* dna = (Nucleotide*)malloc(sizeof(Nucleotide)*len_dna);
  init_conf->seq_x = seq_x;
  init_conf->seq_y = seq_y;
  init_conf->dna = dna;

  // convert the strings
  for (i = 0; i < len_x; i++)
    init_conf->seq_x[i] = convert_aa(seq_x_str[i]);

  for (i = 0; i < len_y; i++)
    init_conf->seq_y[i] = convert_aa(seq_y_str[i]);

  for (i = 0; i < len_dna; i++)
    init_conf->dna[i] = convert_nuc(dna_str[i]);

  init_conf->ugap_x = ungapped_seq_id(init_conf->seq_y, len_x);
  init_conf->ugap_y = ungapped_seq_id(init_conf->seq_x, len_y);
  init_conf->frame_x = frames[0];
  init_conf->frame_y = frames[1];
  init_conf->len_x = len_x;
  init_conf->len_y = len_y;
  init_conf->len_dna = len_dna;
  return init_conf;
}

// free configuration
void free_conf(over_conf_p conf) {
  free(conf->seq_x);
  free(conf->seq_y);
  free(conf->dna);
  free(conf->ugap_x);
  free(conf->ugap_y);
  free(conf);
}

// Propose a random mutation in the dna sequence return the position to change
// in both X and Y the position in the dna seq, and the aa and nuc to change
mutation_p prop_mutation(over_conf_p old_conf, mutation_p tmp_mutation) {
  Nucleotide nuc = (Nucleotide)rand_int(0, 3);
  int pnuc = rand_int(0, old_conf->len_dna);
  if (nuc >= old_conf->dna[pnuc])
    nuc++;
  int mfx = abs(old_conf->frame_x) - 1, mfy = abs(old_conf->frame_y) - 1;
  int codonxi, codonyi, icodx, icody;
  Aminoacid aax, aay;

  // index of the changed codon
  codonxi = (pnuc - mfx)/3;
  codonyi = (pnuc - mfy)/3;

  // index in both codon
  icodx = pnuc - (mfx + codonxi * 3);
  icody = pnuc - (mfy + codonyi * 3);

  // get the mutated types
  aax = get_mut_codon(old_conf->dna, old_conf->frame_x, codonxi, icodx, nuc);
  aay = get_mut_codon(old_conf->dna, old_conf->frame_y, codonyi, icody, nuc);

  // save the proposed mutation
  tmp_mutation->posx = codonxi; tmp_mutation->posy = codonyi;
  tmp_mutation->icodx = icodx; tmp_mutation->icody = icody;
  tmp_mutation->aax = aax; tmp_mutation->aay = aay;
  tmp_mutation->pnuc = pnuc; tmp_mutation->nuc = nuc;

  return tmp_mutation;
}

/* // Propose a random mutation in the dna sequence return the position to change */
/* // in both X and Y the position in the dna seq, and the aa and nuc to change */
/* mutation_p prop_mutation(over_conf_p old_conf, mutation_p tmp_mutation) { */
/*   // Set the number of mutations */
/*   int nb_mutation = rand_int(0, 3); */
/*   mutation_p* list_mutation = (mutation_p*)malloc(sizeof(mutation)*nb_mutation); */
/*   int imut; */
/*   Nucleotide nuc; */
/*   int pnuc; */
/*   int mfx = abs(old_conf->frame_x) - 1, mfy = abs(old_conf->frame_y) - 1; */
/*   int codonxi, codonyi, icodx, icody; */
/*   Aminoacid aax, aay; */

/*   for (imut = 0; imut < nb_mutation; imut++) { */
/*     nuc = (Nucleotide)rand_int(0, 3); */
/*     pnuc = rand_int(0, old_conf->len_dna); */

/*     // index of the changed codon */
/*     codonxi = (pnuc - mfx)/3; */
/*     codonyi = (pnuc - mfy)/3; */

/*     // index in both codon */
/*     icodx = pnuc - (mfx + codonxi * 3); */
/*     icody = pnuc - (mfy + codonyi * 3); */

/*     // get the mutated types */
/*     aax = get_mut_codon(old_conf->dna, old_conf->frame_x, codonxi, icodx, nuc); */
/*     aay = get_mut_codon(old_conf->dna, old_conf->frame_y, codonyi, icody, nuc); */

/*     // save the proposed mutation */
/*     tmp_mutation->posx = codonxi; tmp_mutation->posy = codonyi; */
/*     tmp_mutation->icodx = icodx; tmp_mutation->icody = icody; */
/*     tmp_mutation->aax = aax; tmp_mutation->aay = aay; */
/*     tmp_mutation->pnuc = pnuc; tmp_mutation->nuc = nuc; */
/*     list_mutation[imut] = tmp_mutation; */
/*   } */

/*   return list_mutation; */
/* } */

// apply the propose mutation in dna and proteins
void apply_mutation(over_conf_p old_conf, mutation_p prop_mut) {
  int px = prop_mut->posx;
  int py = prop_mut->posy;
  int pnuc = prop_mut->pnuc;
  Aminoacid ax = prop_mut->aax;
  Aminoacid ay = prop_mut->aay;
  Nucleotide nuc = prop_mut->nuc;
  old_conf->seq_x[px] = ax;
  old_conf->seq_y[py] = ay;
  old_conf->dna[pnuc] = nuc;
}

void print_conf(over_conf_p conf) {
  int i;
  printf("\\SCORE %.3f\n", conf->nrj);
  /* printf("\\DNA   "); */
  /* for (i = 0; i < conf->len_dna; i++) */
  /*   printf("%c", nuc_2_char[conf->dna[i]]); */
  /* printf("\n"); */
  printf("\\SEQX  ");
  for (i = 0; i < conf->len_x; i++)
    printf("%c", amino_acids_char[conf->seq_x[i]]);
  printf("\n");
  printf("\\SEQY  ");
  for (i = 0; i < conf->len_y; i++)
    printf("%c", amino_acids_char[conf->seq_y[i]]);
  printf("\n");
}

// compute the delta nrj for a given mutation
float update_nrj(Aminoacid* seq, int len_seq, float**** current_score, Aminoacid prop_mut, int posi, int* ugap) {
  int i, ui, uj;
  Aminoacid aai, aa_old;
  float delta_nrj = 0.0;
  aa_old = seq[posi];
  ui = ugap[posi];

  if (posi < len_seq) {
    delta_nrj = current_score[ui][ui][prop_mut][prop_mut] - current_score[ui][ui][aa_old][aa_old];
    for (i = 0; i < len_seq; i++) {
      if (i != posi) {
        if (seq[i] != GAP) {
          aai = seq[i];
          // get the ungapped id
          uj = ugap[i];
          delta_nrj += current_score[ui][uj][prop_mut][aai] - current_score[ui][uj][aa_old][aai];
        }
      }
    }
  }
  return delta_nrj;
}

// update the DCA score
// need to save the mutated positions (in the DNA and protein)
// always at least two positions
float update_pair_nrj(over_conf_p current_conf, mutation_p prop_mut, float**** scores_x, float**** scores_y) {
  float new_nrj_conf = 0.0;
  float old_nrj_conf = 0.0;
  int i;

  // compute the change in energy for the mutation in each sequence
  new_nrj_conf = update_nrj(current_conf->seq_x, current_conf->len_x, scores_x, prop_mut->aax, prop_mut->posx, current_conf->ugap_x);
  new_nrj_conf += update_nrj(current_conf->seq_y, current_conf->len_y, scores_y, prop_mut->aay, prop_mut->posy, current_conf->ugap_y);

  // return the change of energy
  return new_nrj_conf;
}

// compute the total score nrj
float compute_nrj(Aminoacid* seq, int len_seq, float**** current_score, int* ugap) {
  float tot_nrj = 0.0;
  int i, j, ui, uj;
  for (i = 0; i < len_seq; i++) {
    if (seq[i] != GAP) {
      ui = ugap[i];
      for (j = i; j < len_seq; j++) {
        if (seq[j] != GAP) {
          uj = ugap[j];
          tot_nrj += current_score[ui][uj][seq[i]][seq[j]];
        }
      }
    }
  }
  return tot_nrj;
}

// MCMC in DNA sequence space
// - For each step, perform a mutation in the DNA sequence, then compute the nrj
// - We start from the sequence optimized with diagonal terms only
void mcmc_run(int nb_steps, float temp, over_conf_p old_conf, float**** scores_x, float**** scores_y) {
  mutation_p current_mut = (mutation_p)malloc(sizeof(mutation));
  float beta = 1.0/temp;
  float nrj_x = 0.0, nrj_y = 0.0;
  int i;
  float delta_nrj = 0.0, rand_val;
  // from initial temp to final temp => cooling scheme
  double final_temp = 0.001, init_temp = temp;
  double mult = exp(log(final_temp/init_temp)/(float)nb_steps);
  nrj_x = compute_nrj(old_conf->seq_x, old_conf->len_x, scores_x, old_conf->ugap_x);
  nrj_y = compute_nrj(old_conf->seq_y, old_conf->len_y, scores_y, old_conf->ugap_y);
  old_conf->nrj = nrj_x + nrj_y;
  printf("# %f init_conf\n", old_conf->nrj);
  print_conf(old_conf);

  // main mcmc loop
  for (i = 0; i < nb_steps; i++) {
    // propose a mutation
    prop_mutation(old_conf, current_mut);

    // update the nrj
    delta_nrj = update_pair_nrj(old_conf, current_mut, scores_x, scores_y);

    // metropolis scheme
    if (delta_nrj < 0.0 || exp(-(1.0/temp) * delta_nrj) >= rand_float(0, 1)) {
      old_conf->nrj += delta_nrj;

      // save the mutation in the actual configuration
      apply_mutation(old_conf, current_mut);
      /* printf("# %d %f %f\n", i, old_conf->nrj, compute_nrj(old_conf->seq_x, old_conf->len_x, scores_x, old_conf->ugap_x) + compute_nrj(old_conf->seq_y, old_conf->len_y, scores_y, old_conf->ugap_y)); */
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
  fasta_seq_p seq_x_el, seq_y_el, dna_el;
  char* seq_x, *seq_y, *dna, *dca_x_file, *dca_y_file, *frame_str;
  arguments.temp = 0.6;
  arguments.nb_steps = 1000;
  argp_parse(&argp, argc, argv, 0, 0, &arguments);
  seq_x_el = read_fasta(arguments.args[0]);
  seq_y_el = read_fasta(arguments.args[1]);
  dna_el = read_fasta(arguments.args[2]);
  seq_x = seq_x_el->seq;
  seq_y = seq_y_el->seq;
  dna = dna_el->seq;
  dca_x_file = arguments.args[3];
  dca_y_file = arguments.args[4];
  frames = convert_frame(arguments.frame);

  over_conf_p init_conf = initialize(seq_x, seq_y, dna, frames);
  float**** scores_x = read_score(init_conf->len_x, dca_x_file, init_conf->frame_x);
  float**** scores_y = read_score(init_conf->len_y, dca_y_file, init_conf->frame_y);

  printf("# FRAME = %s\n# TEMP = %f\n# NB_STEPS = %d\n", arguments.frame, arguments.temp, arguments.nb_steps);
  mcmc_run(arguments.nb_steps, arguments.temp, init_conf, scores_x, scores_y);

  free_score(scores_x, init_conf->len_x);
  free_score(scores_y, init_conf->len_y);
  free_conf(init_conf);
  return 0;
}

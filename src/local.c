/**
 *   \file local.c
 *   \brief Compute overlapping genes
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <argp.h>
#include "seq_utils.h"

// Gap symbol
#define MINF -9999999999.0
#define MAX_CHAR 3000

// macro max
#define MAX(x, y) (x->score >= y->score ? x:y)

// Store all the kmers with the corresponding translated residues
Bool SIM = true;
// perform global search
float GAP_PEN = -0.5;
float GAP_OPN = -10.0;

// score threshold value
float THR_SCR = -10.0;
int MAT_TYPE = 0; // Means Identity matrix:: 1=blosum62, 2=blosum80, 3=blosum90, 4=vtml
quartet_p**** KMER_LIST;
char amino_acids_char[] = {'A' ,'C' ,'T' ,'E' ,'D' ,'F' ,'W' ,'I' ,'V' ,'L' ,'K' ,'M' ,'N' ,'Q' ,'S' ,'R' ,'Y' ,'H' ,'P' ,'G', '*', '-'};
char nuc_2_char[] = {'A','C','T','G'};

// Parse arguments -------------------------------------------------------------
const char *argp_program_version = "Version 1.0";
const char *argp_program_bug_address = "@polytechnique.edu";

/* Program documentation. */
static char doc[] = "Overlapping gene search";

/* A description of the arguments we accept. */
static char args_doc[] = "<seq_x> <seq_y>";

/* The options we understand. */
static struct argp_option options[] = {
                                       {"phase", 'f', "<-2|-1|0|1|2>", 0, "Overlapping frame, negative frames are n2 and n1", 0},
                                       {"gap_cost", 'g', "INT", 0, "Gap cost value", 0 },
                                       {"gap_open", 'o', "INT", 0, "Gap open value", 0 },
                                       {"thres_scr", 't', "INT", 0, "Threshold score value", 0 },
                                       {"matrix",   's', "<ident|blosum62|blosum80|blosum90|vtml200>", 0, "Substitution matrix to use", 0 },
                                       { 0 }
};

/* Used by main to communicate with parse_opt. */
struct arguments {
  char *args[2]; char *frame;
  char *matrix; float gap_cost;
  float gap_open; float thr_scr;
};


/* Parse a single option. */
error_t parse_opt (int key, char *arg, struct argp_state *state) {
  /* Get the input argument from argp_parse, which we
     know is a pointer to our arguments structure. */
  struct arguments *arguments = state->input;

  switch (key)
    {
    case 'f': arguments->frame = (char*)malloc(sizeof(char)*MAX_CHAR); strcpy(arguments->frame, arg); break;
    case 'g': arguments->gap_cost = atof(arg); break;
    case 'o': arguments->gap_open = atof(arg); break;
    case 't': arguments->thr_scr = atof(arg); break;
    case 's': arguments->matrix = (char*)malloc(sizeof(char)*MAX_CHAR); strcpy(arguments->matrix, arg); break;

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

static struct argp argp = { options, parse_opt, args_doc, doc, 0, 0, 0};
// Parse arguments -------------------------------------------------------------

/**
 *  \brief connect one matrix entry to the previous entry
 *  \param previous quartet_entry, current quartet_entry, is there a gap
 *
 *  \details: we connect the current quartet_entry with the previous one.
 *  Therefore, we need to loop in the kmer list and find the best quartet for
 *  the current pair of aa.
 *
 *  \return quartet_unit
 */
quartet_p* connect(matrix_entry_p prev_matrix, matrix_entry_p cur_matrix,
                   Bool gap, pair_p cur_pair) {
  int i, j, k, l;
  quartet_p* max_quartet = (quartet_p*)malloc(sizeof(quartet_p)*4);
  quartet_p prev_pos = NULL, tmp_kmer = NULL, cur_pos = NULL;
  float new_score;

  // Initialize the maximum quartets
  for (i = A; i <= G; ++i) {
    max_quartet[i] = (quartet_p)malloc(sizeof(quartet));
    max_quartet[i]->score = -INFINITY;
  }

  // XXX: surely a better way to do this
  for (i = A; i <= G; ++i) {
    for (j = A; j <= G; ++j) {
      for (k = A; k <= G; ++k) {
        for (l = A; l <= G; ++l) {
          // save the kmer translated pair
          tmp_kmer = KMER_LIST[i][j][k][l];
          cur_pos = cur_matrix->quartet_entries[l];

          new_score = compute_score_sim_local(cur_pair, tmp_kmer->prime, MAT_TYPE, THR_SCR);

          // if previous entry is not NULL (0, 0)
          if (prev_matrix) {
            // matching i —> previous l
            prev_pos = prev_matrix->quartet_entries[i];

            // add gap cost and history
            if (gap == true) {
              new_score += GAP_PEN;
              if (prev_pos->gap == false) {
                new_score += GAP_OPN;
              }
            }
            new_score += prev_pos->score;

          } else
            // no previous position
            prev_pos = NULL;

          // If score is better, save the quartet
          if (max_quartet[l]->score < new_score) {
            // create a new quartet object
            max_quartet[l]->kmer = tmp_kmer->kmer;
            max_quartet[l]->prime = tmp_kmer->prime;
            max_quartet[l]->init = cur_pair;
            // set negative scores to zero
            max_quartet[l]->score = new_score * (new_score > 0.0);
            max_quartet[l]->previous = (struct quartet*)prev_pos;
            max_quartet[l]->gap = gap;
          }
        }
      }
    }
  }

  return max_quartet;
}

/**
 *  \brief Takes two sequences, compute the non-continuous overlapping optimization
 *  \param seq_a = string, seq_b = string, frames = frame_list
 *  \details: We take as an input, matrix already filled up with matrix entries
 *  \return quartet_p
 */
result_quartet_p overlapping_nc(matrix_qp matrix) {
  int i, j, nuc, imax, jmax;
  quartet_p* substitution;
  quartet_p* deletion;
  quartet_p* insertion;
  quartet_p max_quartet = NULL;
  quartet_p tmp_max;
  // array of possibilities
  pair_p subs_pair, delt_pair, insr_pair;
  result_quartet_p solution = malloc(sizeof(result_quartet));

  // connect all the matrix_entries
  for (i = 1; i <= matrix->len_x; ++i) {
    for (j = 1; j <= matrix->len_y; ++j) {

      subs_pair = (pair_p)malloc(sizeof(pair)); delt_pair = (pair_p)malloc(sizeof(pair)); insr_pair = (pair_p)malloc(sizeof(pair));
      subs_pair->x = matrix->matrix_entries[i][j]->quartet_entries[A]->init->x;
      subs_pair->y = matrix->matrix_entries[i][j]->quartet_entries[A]->init->y;

      delt_pair->x = GAP;
      delt_pair->y = matrix->matrix_entries[i][j]->quartet_entries[A]->init->y;

      insr_pair->x = matrix->matrix_entries[i][j]->quartet_entries[A]->init->x;
      insr_pair->y = GAP;

      // Mutation, deletion or substitution move
      substitution = connect(matrix->matrix_entries[i - 1][j - 1], matrix->matrix_entries[i][j], false, subs_pair);
      deletion = connect(matrix->matrix_entries[i][j - 1], matrix->matrix_entries[i][j], true, delt_pair);
      insertion = connect(matrix->matrix_entries[i - 1][j], matrix->matrix_entries[i][j], true, insr_pair);

      // for each nuc, find the best previous matrix entry
      // to do so, we allocate temporarily matrix entries for substitution, deletion and insertion
      for (nuc = A; nuc <= G; nuc++) {
        tmp_max = get_max_quartet(substitution[nuc], deletion[nuc], insertion[nuc]);

        // copy values values
        copy_to_max(matrix->matrix_entries[i][j]->quartet_entries[nuc], tmp_max);
        /* printf("\\\\ %d %d %c %c %c %f\n", i, j, nuc_2_char[nuc], amino_acids_char[matrix->matrix_entries[i][j]->quartet_entries[nuc]->init->x], amino_acids_char[matrix->matrix_entries[i][j]->quartet_entries[nuc]->init->y], matrix->matrix_entries[i][j]->quartet_entries[nuc]->score); */

        if (max_quartet == NULL || max_quartet->score < matrix->matrix_entries[i][j]->quartet_entries[nuc]->score) {
          imax = i; jmax = j;
          max_quartet = matrix->matrix_entries[i][j]->quartet_entries[nuc];
        }
        free(substitution[nuc]); free(deletion[nuc]); free(insertion[nuc]);
      }
      free(substitution); free(deletion); free(insertion);
      free(subs_pair); free(delt_pair); free(insr_pair);
    }
  }
  solution->result = max_quartet;
  solution->posx = imax;
  solution->posy = jmax;
  return solution;
};

int main(int argc, char *argv[]) {
  struct arguments arguments;
  int* frames;
  char *seq_x, *seq_y;
  fasta_seq_p seq_x_el, seq_y_el;
  // default arguments
  arguments.frame = "-2";
  arguments.gap_cost = -1.0;
  arguments.gap_open = -2.0;
  arguments.thr_scr = 0.0;
  arguments.matrix = "ident";

  // Parse the command line
  argp_parse(&argp, argc, argv, 0, NULL, &arguments);
  GAP_PEN = arguments.gap_cost;
  GAP_OPN = arguments.gap_open;
  MAT_TYPE = convert_matrix(arguments.matrix);
  THR_SCR = arguments.thr_scr;
  frames = convert_frame(arguments.frame);
  /* printf("# %d %d %d \n" ,frames[0], frames[1], MAT_TYPE); */

  seq_x_el = read_fasta(arguments.args[0]);
  seq_y_el = read_fasta(arguments.args[1]);
  seq_x = seq_x_el->seq;
  seq_y = seq_y_el->seq;

  if (frames[0] < 0)
    reverse(seq_x, 0, strlen(seq_x)-1);
  if (frames[1] < 0)
    reverse(seq_y, 0, strlen(seq_y)-1);

  // Prepare the system
  matrix_qp matrix = init_matrix(seq_x, seq_y);
  result_quartet_p results;

  KMER_LIST = get_all_kmers(frames);

  results = overlapping_nc(matrix);

  printf("# GAP_PEN = %.1f\n# GAP_OPEN = %.1f\n# MAT = %s\n# FRAME = %s\n", arguments.gap_cost, arguments.gap_open, arguments.matrix, arguments.frame);
  printf("# NAMEX = %s\n# NAMEY = %s\n", seq_x_el->name, seq_y_el->name);
  printf("# SCORE = %.1f\n", results->result->score);
  if (frames[0] < 0)
    reverse(seq_x, 0, strlen(seq_x)-1);
  if (frames[1] < 0)
    reverse(seq_y, 0, strlen(seq_y)-1);
  backtrace(results, seq_x, seq_y, frames, false);

  free(seq_x); free(seq_y);
  free_kmer_list(KMER_LIST);
  free_matrix(matrix);
  free(results->result);
  free(results);
  return 0;
}

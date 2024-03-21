/**
 *   \file global.c
 *   \brief Compute overlapping genes
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <argp.h>
#include "seq_utils.h"

// Gap symbol
#define MINF -1.0/0.0
#define MAX_CHAR 3000

// macro max
#define MAX(x, y) (x->score >= y->score ? x:y)

// Store all the kmers with the corresponding translated residues
Bool SIM = true;
// perform global search
Bool DCA = false;
float GAP_PEN = -0.5;
float GAP_OPN = -10.0;
int MAT_TYPE = 0; // Means Identity matrix:: 1=blosum62, 2=blosum80, 3=blosum90, 4=vtml
quartet_p**** KMER_LIST;
float**** dca_scores_x; // dca scores
float**** dca_scores_y; // dca scores
char amino_acids_char[] = {'A' ,'C' ,'T' ,'E' ,'D' ,'F' ,'W' ,'I' ,'V' ,'L' ,'K' ,'M' ,'N' ,'Q' ,'S' ,'R' ,'Y' ,'H' ,'P' ,'G', '*', '-'};
char nuc_2_char[] = {'A','C','T','G'};

// Parse arguments -------------------------------------------------------------
const char *argp_program_version = "Version 1.0";
const char *argp_program_bug_address = "@polytechnique.edu";

/* Program documentation. */
static char doc[] = "Overlapping gene search or design";

/* A description of the arguments we accept. */
static char args_doc[] = "<seq_x> <seq_y>";

/* The options we understand. */
static struct argp_option options[] = {
                                       {"phase", 'f', "<-2|-1|0|1|2>", 0, "Overlapping frame, negative frames are n2 and n1", 0},
                                       {"score_x", 'x', "FILE", 0, "seq_x dca scores", 0 },
                                       {"score_y", 'y', "FILE", 0, "seq_y dca scores", 0 },
                                       {"gap_cost", 'g', "INT", 0, "Gap cost value", 0 },
                                       {"gap_open", 'o', "INT", 0, "Gap open value", 0 },
                                       {"matrix",   's', "<ident|blosum62|blosum80|blosum90|vtml200>", 0, "Substitution matrix to use", 0 },
                                       { 0 }
};

/* Used by main to communicate with parse_opt. */
struct arguments {
  char *args[2];                /* seqx & seqy */
  char *dca_file_x;
  char *dca_file_y;
  char *frame;
  char *matrix;
  float gap_cost;
  float gap_open;
};

/* Parse a single option. */
error_t parse_opt (int key, char *arg, struct argp_state *state) {
  /* Get the input argument from argp_parse, which we
     know is a pointer to our arguments structure. */
  struct arguments *arguments = state->input;

  switch (key)
    {
    case 'f': arguments->frame = (char*)malloc(sizeof(char)*MAX_CHAR); strcpy(arguments->frame, arg); break;
    case 'x': arguments->dca_file_x = (char*)malloc(sizeof(char)*MAX_CHAR); strcpy(arguments->dca_file_x, arg); break;
    case 'y': arguments->dca_file_y = (char*)malloc(sizeof(char)*MAX_CHAR); strcpy(arguments->dca_file_y, arg); break;
    case 'g': arguments->gap_cost = atof(arg); break;
    case 'o': arguments->gap_open = atof(arg); break;
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
                   Bool gap, pair_p cur_pair, Bool last_pos) {
  int i, j, k, l;
  quartet_p* max_quartet = (quartet_p*)malloc(sizeof(quartet_p)*4);
  quartet_p prev_pos = NULL, tmp_kmer = NULL, cur_pos = NULL;
  float new_score;

  // Initialize the maximum quartets
  for (i = A; i <= G; ++i) {
    max_quartet[i] = (quartet_p)malloc(sizeof(quartet));
    max_quartet[i]->init = NULL;
  }

  // XXX: surely a better way to do this
  for (i = A; i <= G; ++i) {
    for (j = A; j <= G; ++j) {
      for (k = A; k <= G; ++k) {
        for (l = A; l <= G; ++l) {
          // save the kmer translated pair
          tmp_kmer = KMER_LIST[i][j][k][l];
          cur_pos = cur_matrix->quartet_entries[l];
          if (DCA == true)
            new_score = compute_score_dca(cur_pair, tmp_kmer->prime,
                                          dca_scores_x, dca_scores_y,
                                          cur_matrix->posx, cur_matrix->posy);
          else
            new_score = compute_score_sim(cur_pair, tmp_kmer->prime, MAT_TYPE);
          // if previous entry is not NULL (0, 0)
          if (prev_matrix) {
            // matching i —> previous l
            // prev_pos = quartet_p
            prev_pos = prev_matrix->quartet_entries[i];

            // if dca, then add the contribution from the previous pair
            if (DCA == true)
              new_score += compute_score_dca_pair(cur_pair, tmp_kmer->prime,
                                                  prev_pos->init, prev_pos->prime,
                                                  dca_scores_x, dca_scores_y,
                                                  cur_matrix, prev_matrix);

            // add gap cost and history
            if (gap == true) {
              new_score += GAP_PEN;
              /* if (prev_pos->gap == false && prev_pos->last_gap == false) */
              if (prev_pos->gap == false)
                new_score += GAP_OPN;
            }
            new_score += prev_pos->score;
          } else
            // no previous position
            prev_pos = NULL;

          // If score is better, save the quartet
          if (max_quartet[l]->init == NULL || max_quartet[l]->score < new_score) {
            // create a new quartet object
            max_quartet[l]->kmer = tmp_kmer->kmer;
            max_quartet[l]->prime = tmp_kmer->prime;
            max_quartet[l]->init = cur_pair;
            max_quartet[l]->score = new_score;
            max_quartet[l]->previous = (struct quartet*)prev_pos;
            max_quartet[l]->gap = gap;
            max_quartet[l]->last_gap = last_pos;
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
  int i, j, nuc;
  quartet_p* substitution;
  quartet_p* deletion;
  quartet_p* insertion;
  quartet_p max_quartet = NULL;
  quartet_p tmp_max;
  // array of possibilities
  pair_p subs_pair, delt_pair, insr_pair;
  Bool last_pos;
  result_quartet_p solution = malloc(sizeof(result_quartet));

  // init zeros
  for (i = 1; i <= matrix->len_x; ++i) {
    insr_pair = (pair_p)malloc(sizeof(pair));
    insr_pair->x = matrix->matrix_entries[i][0]->quartet_entries[A]->init->x;
    insr_pair->y = GAP;
    insertion = connect(matrix->matrix_entries[i-1][0], matrix->matrix_entries[i][0], true, insr_pair, true);
    for (nuc = A; nuc <= G; nuc++)
      copy_to_max(matrix->matrix_entries[i][0]->quartet_entries[nuc], insertion[nuc]);
    free(insertion);
  }

  for (j = 1; j <= matrix->len_y; ++j) {
    delt_pair = (pair_p)malloc(sizeof(pair));
    delt_pair->x = GAP;
    delt_pair->y = matrix->matrix_entries[0][j]->quartet_entries[A]->init->y;
    deletion = connect(matrix->matrix_entries[0][j-1], matrix->matrix_entries[0][j], true, delt_pair, true);
    for (nuc = A; nuc <= G; nuc++)
      copy_to_max(matrix->matrix_entries[0][j]->quartet_entries[nuc], deletion[nuc]);
    free(deletion);
  }

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
      last_pos = ((i == matrix->len_x || j == matrix->len_y) ? true : false);

      // Mutation, deletion or substitution move
      substitution = connect(matrix->matrix_entries[i - 1][j - 1], matrix->matrix_entries[i][j], false, subs_pair, last_pos);
      deletion = connect(matrix->matrix_entries[i][j - 1], matrix->matrix_entries[i][j], true, delt_pair, last_pos);
      insertion = connect(matrix->matrix_entries[i - 1][j], matrix->matrix_entries[i][j], true, insr_pair, last_pos);

      // for each nuc, find the best previous matrix entry
      // to do so, we allocate temporarily matrix entries for substitution, deletion and insertion
      for (nuc = A; nuc <= G; nuc++) {
        tmp_max = get_max_quartet(substitution[nuc], deletion[nuc], insertion[nuc]);

        // copy values values
        copy_to_max(matrix->matrix_entries[i][j]->quartet_entries[nuc], tmp_max);

        // save the best quartet for each nuc type
        if (max_quartet == NULL || max_quartet->score < matrix->matrix_entries[i][j]->quartet_entries[nuc]->score) {
          max_quartet = matrix->matrix_entries[i][j]->quartet_entries[nuc];
        }
        free(substitution[nuc]); free(deletion[nuc]); free(insertion[nuc]);
      }
      free(substitution); free(deletion); free(insertion);
      free(subs_pair); free(delt_pair); free(insr_pair);
    }
  }

  // get the best of the four quartet
  max_quartet = matrix->matrix_entries[matrix->len_x][matrix->len_y]->quartet_entries[A];
  for (nuc = A; nuc <= G; nuc++) {
    if (max_quartet->score < matrix->matrix_entries[matrix->len_x][matrix->len_y]->quartet_entries[nuc]->score) {
      max_quartet = matrix->matrix_entries[matrix->len_x][matrix->len_y]->quartet_entries[nuc];
    }
  }

  solution->result = max_quartet;
  solution->posx = i-1;
  solution->posy = j-1;
  return solution;
};

int main(int argc, char *argv[]) {
  struct arguments arguments;
  int* frames;
  char *seq_x, *seq_y;
  fasta_seq_p seq_x_el, seq_y_el;
  // default arguments
  arguments.frame = "-2";
  arguments.gap_cost = -3.0;
  arguments.gap_open = -16.0;
  arguments.dca_file_x = "-";
  arguments.dca_file_y = "-";
  arguments.matrix = "ident";

  // Parse the command line
  argp_parse(&argp, argc, argv, 0, NULL, &arguments);
  GAP_PEN = arguments.gap_cost;
  GAP_OPN = arguments.gap_open;
  MAT_TYPE = convert_matrix(arguments.matrix);
  frames = convert_frame(arguments.frame);

  // read sequences from fasta
  seq_x_el = read_fasta(arguments.args[0]);
  seq_y_el = read_fasta(arguments.args[1]);
  seq_x = seq_x_el->seq;
  seq_y = seq_y_el->seq;

  // Read DCA scores
  if (strcmp(arguments.dca_file_x, "-") != 0 && strcmp(arguments.dca_file_y, "-") != 0) {
    DCA = true;
    dca_scores_x = read_score(strlen(seq_x), arguments.dca_file_x, frames[0]);
    dca_scores_y = read_score(strlen(seq_y), arguments.dca_file_y, frames[1]);
  }

  if (frames[0] < 0)
    reverse(seq_x, 0, strlen(seq_x)-1);
  if (frames[1] < 0)
    reverse(seq_y, 0, strlen(seq_y)-1);

  printf("# GAP_PEN = %.1f\n# GAP_OPEN = %.1f\n# MAT = %s\n# FRAME = %s\n", arguments.gap_cost, arguments.gap_open, arguments.matrix, arguments.frame);
  printf("# NAMEX = %s\n# NAMEY = %s\n", seq_x_el->name, seq_y_el->name);
  if (DCA == true)
    printf("# DCA = true\n");
  else
    printf("# DCA = false\n");

  // Prepare the system
  matrix_qp matrix = init_matrix(seq_x, seq_y);
  result_quartet_p results;

  KMER_LIST = get_all_kmers(frames);

  results = overlapping_nc(matrix);
  if (frames[0] < 0)
    reverse(seq_x, 0, strlen(seq_x)-1);
  if (frames[1] < 0)
    reverse(seq_y, 0, strlen(seq_y)-1);

  printf("# SCORE = %.1f\n", results->result->score);
  backtrace(results, seq_x, seq_y, frames, true);

  free(seq_x); free(seq_y);
  free_kmer_list(KMER_LIST);
  if (DCA == true) {
    free_score(dca_scores_x, matrix->len_x);
    free_score(dca_scores_y, matrix->len_y);
  }
  free_matrix(matrix);
  free(results->result);
  free(results);
  return 0;
}

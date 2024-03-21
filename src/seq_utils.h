#ifndef SEQ_UTILS_

#include <stdio.h>
#include <string.h>

// Global variable
typedef enum {true, false} Bool;
typedef enum {A, C, T, G} Nucleotide;
typedef enum {ALA ,CYS ,THR ,GLU ,ASP ,PHE ,TRP ,ILE ,VAL ,LEU ,LYS ,
              MET ,ASN ,GLN ,SER ,ARG ,TYR ,HIS ,PRO ,GLY, STOP, GAP
} Aminoacid;

// save sequence and name
typedef struct {
  char* name;
  char* seq;
} fasta_seq;
typedef fasta_seq *fasta_seq_p;

// pair of residues
typedef struct {
  Aminoacid x;
  Aminoacid y;
} pair;
typedef pair *pair_p;

// Store sequences
typedef struct {
  char* name;
  char* seq;
} seq_str;
typedef seq_str *seq_str_p;

// one quartet
typedef struct {
  pair_p prime;
  pair_p init;
  float score;
  Nucleotide* kmer;
  Nucleotide* ckmer;
  struct quartet* previous;
  Bool gap;
  Bool last_gap;
} quartet;
typedef quartet *quartet_p;

// Save the best quartet
typedef struct {
  quartet_p result;
  int posx; int posy;
} result_quartet;
typedef result_quartet *result_quartet_p;

// 4 quartets, entry of the dynamic programming matrix
typedef struct {
  quartet_p* quartet_entries;
  Bool gap;
  int posx;
  int posy;
} matrix_entry;
typedef matrix_entry *matrix_entry_p;

// Matrix with 4 quartet entries
typedef struct {
  int len_x;
  int len_y;
  matrix_entry_p **matrix_entries;
} matrix_q;
typedef matrix_q *matrix_qp;

// Functions

// Init the list of 256 quartets
quartet_p**** get_all_kmers(int* frames);

// Init the dynamic programming matrix
matrix_qp init_matrix(char* seq_x, char* seq_y);

// free kmer list
void free_kmer_list(quartet_p**** kmer_list);

// free the matrix
void free_matrix(matrix_qp kmer_list);

// compute the similarity (blosum62) between two pairs of aa
float compute_score_sim(pair_p init, pair_p prime, int mat_type);

// compute the similarity score for local search
float compute_score_sim_local(pair_p init, pair_p prime, int mat_type, float thr_scr);

// get the quartet that maximizes that score
quartet_p get_max_quartet(quartet_p substitution, quartet_p deletion, quartet_p insertion);

// copy to the max quartet
quartet_p copy_to_max(quartet_p from_quartet, quartet_p to_quartet);

// read the matrix scores
float**** read_score(int len_seq, char* dca_score_file, int frame);

// free the dca scores
void free_score(float**** dca_score, int len_seq);

// compute the diagonal dca score
float compute_score_dca(pair_p init, pair_p prime, float**** mat_dca_x, float**** mat_dca_y, int pos_x, int pos_y);

// compute diagonal terms for connected quartet
float compute_score_dca_pair(pair_p cur_init, pair_p cur_prime, pair_p prev_init, pair_p prev_prime,
                             float**** mat_dca_x, float**** mat_dca_y, matrix_entry_p cur_mat,
                             matrix_entry_p prev_mat);
// reverse a string
void reverse(char *x, int begin, int end);

// reverse a string
int rand_int(int from_int, int to_int);

// random float
float rand_float(int min, int max);

// get a codon from a dna sequence and a given position
Aminoacid get_mut_codon(Nucleotide* dna, int frame, int posi, int icod, Nucleotide nnuc);

// return the ungapped id of a sequence
int* ungapped_seq_id(Aminoacid* seq, int len_seq);

// convert to enum
Nucleotide convert_nuc(char x);
Aminoacid convert_aa(char x);

// for parsing the frame str
int* convert_frame(char* frame);

// for parsing the matrix str
int convert_matrix(char* matrix_str);

// read a fasta file
fasta_seq_p read_fasta(char* fasta_file_str);

// progress infos
void printProgress (int step, int traj_len);

// backtrace to print sequences
void backtrace(result_quartet_p solution, char* seq_x_full, char* seq_y_full, int* frames, Bool glob);

#endif // SEQ_UTILS_

#ifndef DESIGN_

// pair of residues
typedef struct {
  Aminoacid* seq_x;
  Aminoacid* seq_y;
  Nucleotide* dna;

  // ungapped id sequence
  int* ugap_x;
  int* ugap_y;
  int len_x;
  int len_y;
  int len_dna;
  int frame_x;
  int frame_y;
  float nrj;
} over_conf;
typedef over_conf *over_conf_p;

// pair of residues
typedef struct {
  // position to change in the DNA sequence
  int pnuc;
  // corresponding codon in each overlapping frame
  // equivalent to the position in the aa sequence
  int posx;
  int posy;
  // proposed types
  Aminoacid aax;
  Aminoacid aay;
  // index in the codon X|Y, can be 0, 1, or 2
  int icodx;
  int icody;
  // the proposed nucleotide
  Nucleotide nuc;
} mutation;
typedef mutation *mutation_p;


/* over_conf_p initialize(char* seq_x,char* seq_y,char* dna,int frame_x,int frame_y); */
/* over_conf_p initialize(char* seq_x_str, char* seq_y_str, char* dna_str, int* frames) { */

// store the proposed mutation
/* Mutation prop_mutation(old_conf, tmp_mutation); */

// perform the mcmc for a given over pair
/* void mcmc_run(int nb_steps, over_conf old_conf); */

#endif // DESIGN_

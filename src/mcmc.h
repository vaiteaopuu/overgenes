#ifndef MCMC_

// Sequence configuration
typedef struct {
  Aminoacid* seq;
  int len;
  int count;
  double nrj;
} mc_conf;
typedef mc_conf *mc_conf_p;

// Proposed mutation in MCMC
typedef struct {
  // position to mutate
  int pos;
  // proposed types
  Aminoacid aa;
  // the proposed nucleotide
} mutation_mc;
typedef mutation_mc *mutation_mc_p;

#endif // MCMC_

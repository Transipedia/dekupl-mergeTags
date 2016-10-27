#include <stdio.h>  //fprintf
#include <stdlib.h> //free
#include <zlib.h>
#include <inttypes.h>
#include <math.h> // pow()
#include <stdint.h>

#include "kstring.h"
#include "kseq.h"
#include "khash.h"
#include "dna.h"

#define VERSION "0.0.1"

#define VIRGIN_KMER 0
#define ASSEMBLED_KMER 1

KSEQ_INIT(gzFile, gzread)

// Init k-mers hash
KHASH_MAP_INIT_INT64(kmers, uint32_t)

typedef khash_t(kmers) kmers_hash_t;

int load_kmers(const char* counts_file, size_t k_length, kmers_hash_t *h) {
  int nb_kmers = 0, dret = 0;
  gzFile fp;
	kstream_t *ks;
	kstring_t *str;
  khiter_t k;
  str = calloc(1, sizeof(kstring_t));

  fp = gzopen(counts_file, "r");
  if(!fp) { fprintf(stderr, "Failed to open %s\n", counts_file); exit(EXIT_FAILURE); }

  ks = ks_init(fp);
  while (ks_getuntil(ks, 0, str, &dret) >= 0) {
    //uint64_t hashed_kmer = hash_kmer(str->s, k_length);
    uint64_t hashed_kmer = dna_to_int(str->s, k_length);
    //fprintf(stderr, "%s => %" PRIu64 "\n", str->s, hashed_kmer);

    k = kh_put(kmers, h, hashed_kmer, &dret);
    kh_value(h, k) = VIRGIN_KMER;

    // skip the rest of the line
		if (dret != '\n') while ((dret = ks_getc(ks)) > 0 && dret != '\n');
    nb_kmers++;
  }

  ks_destroy(ks);
  gzclose(fp);
  free(str->s); free(str);
  return nb_kmers;
}

char *extend_kmer(uint64_t kmer, size_t k_length, kmers_hash_t *h, int reverse) {
  khiter_t k;
  kstring_t *assembly = calloc(1, sizeof(kstring_t));
  char *kmer_str = (char*)calloc(k_length + 1, 1);

  uint64_t neighbor_kmer, match_kmer, assembly_kmer;
  khiter_t match_k;
  int nb_match = 1;
  char nuc;

  if(reverse)
    assembly_kmer = int_revcomp(kmer,k_length);
  else
    assembly_kmer = kmer;

  int_to_dna(assembly_kmer, k_length, kmer_str);
  kputs(kmer_str,assembly);

  while (nb_match == 1) {
    nb_match = 0;

    for(int p = 0; p < NB_NUCLEOTIDES; p++) {
      neighbor_kmer = next_kmer(k_length, NUCLEOTIDES[p], assembly_kmer);
      if(reverse)
        k = kh_get(kmers, h, int_revcomp(neighbor_kmer,k_length));
      else
        k = kh_get(kmers, h, neighbor_kmer);
      if(k != kh_end(h) && kh_value(h,k) == VIRGIN_KMER) {
        nb_match++;
        nuc = NUCLEOTIDES[p];
        match_k = k;
        match_kmer = neighbor_kmer;
      }
      if(nb_match > 1)
        break;
    }

    if(nb_match == 1) {
      kh_value(h, match_k) = ASSEMBLED_KMER;
      kputc(nuc,assembly);
      assembly_kmer = match_kmer;
    }
  }

  char *return_value = ks_release(assembly);
  free(assembly);

  if(reverse) {
    char *reverse_assembly = malloc(strlen(return_value) + 1);
    revcomp(return_value, reverse_assembly, strlen(return_value));
    free(return_value);
    return_value = reverse_assembly;
  }

  return return_value;
}

void assemble_kmers(kmers_hash_t *h, size_t k_length) {
  khiter_t k;
  for(k = kh_begin(h); k != kh_end(h); ++k) {
    if(!kh_exist(h, k)) continue;
    if(kh_value(h,k) == VIRGIN_KMER) {
      kh_value(h, k) = ASSEMBLED_KMER;

      char *right_assembly  = extend_kmer(kh_key(h, k), k_length, h, 0);
      char *left_assembly   = extend_kmer(kh_key(h, k), k_length, h, 1);

      // Merge left and right assemblies
      char *assembly = malloc(strlen(left_assembly) + strlen(right_assembly) - k_length + 1);
      strcpy(assembly, left_assembly);
      strcat(assembly, &right_assembly[k_length]);

      //fprintf(stderr, "LEFT:  %s\n", left_assembly);
      //fprintf(stderr, "RIGHT: %s\n", right_assembly);
      fprintf(stdout, "%s\n", assembly);

      free(assembly);
      free(right_assembly);
      free(left_assembly);
    }
  }
}

int main(int argc, char *argv[])
{
	//fprintf(stderr, "[TOTO]\n");
  char *counts_file;
  int k_length = 31;

  if (optind == argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   mergeTags [options] <counts.tsv>\n\n");
		fprintf(stderr, "Options: -k INT    length of k-mers (max_value: 32) [%d]\n", k_length);
		fprintf(stderr, "\n");
		return 1;
	}

  counts_file = argv[optind++];
  khash_t(kmers) *h = kh_init(kmers);

  fprintf(stderr, "Loading k-mers into memory\n");
  int nb_kmers = load_kmers(counts_file,k_length,h);
  fprintf(stderr, "%d k-mers loaded\n", nb_kmers);
  fprintf(stderr, "Assembling k-mers\n");
  assemble_kmers(h,k_length);
	return 0;
}

/*
 * Quantify bycatch around sites targeted by enrichment.
 *
 * Copyright (c) 2015 Graham Gower <graham.gower@gmail.com>
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 */

/* Parts of this code derived from samools:

    Copyright (C) 2011 Broad Institute.
    Copyright (C) 2013, 2014 Genome Research Ltd.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <fcntl.h>
#include <stdint.h>
#include <errno.h>
#include <getopt.h>
#include <math.h>

#include <bam.h>

#include <zlib.h>
#include "kseq.h"
KSTREAM_INIT(int, read, 16384);

#include "klist.h"
#include "kbtree.h"

#define max(a, b) ((a)>(b)?(a):(b))


typedef struct {
	uint64_t id; // (tid << 32) | pos
	unsigned int window_cov; // sum of the site coverage,
				// for each site in the window around the locus
	unsigned int max_cov; // max coverage of all the sites
	unsigned int site_cov; // coverage at the site defining the window

	unsigned int n; // nonzero sites in the window
	double m; // cumulative mean
	double s; // cumulative variance*(n-1)
	double md; // median
} locus_t;

#define loci_cmp(a, b) (((b).id < (a).id) - ((a).id < (b).id))
KBTREE_INIT(loci, locus_t, loci_cmp);
#define pointer_free(x)
KLIST_INIT(loci, locus_t*, pointer_free);

struct bycatch {
	char *bed_file;
	char *bam_file;
	kbtree_t(loci) *hi; // bed loci
	bamFile bam_fp;
	bam_hdr_t* bam_hdr;

	unsigned int min_qual;
	unsigned int window_size;
};

/*
 * Load .bed loci into a B-tree.
 * Takes the chromosome name from the first column, matching it to one from
 * the bam header, and the position from the second column.
 */
static kbtree_t(loci) *
loadpos(const char *bed_fn, bam_hdr_t *bam_hdr)
{
	int fd;
	int dret;
	kstream_t *ks;
	kstring_t *str;
	kbtree_t(loci) *hash;

	fd = open(bed_fn, O_RDONLY);
	if (fd == -1) {
		fprintf(stderr, "open: %s: %s\n", bed_fn, strerror(errno));
		return NULL;
	}

	hash = kb_init(loci, KB_DEFAULT_SIZE);
	str = calloc(1, sizeof(kstring_t));
	ks = ks_init(fd);

	while (ks_getuntil(ks, KS_SEP_SPACE, str, &dret) >= 0) {
		int tid = bam_name2id(bam_hdr, str->s);
		if (tid >= 0 && dret != '\n') {
			if (ks_getuntil(ks, KS_SEP_SPACE, str, &dret) >= 0) {
				uint64_t x = (uint64_t)tid<<32 | atoi(str->s);
				locus_t l = {x, 0};
				kb_putp(loci, hash, &l);
			} else
				break;
		}
		if (dret != '\n') {
			// consume the rest of the line
			while ((dret = ks_getc(ks)) > 0 && dret != '\n')
				;
		}
		if (dret < 0)
			break;
	}
	ks_destroy(ks);
	close(fd);
	free(str->s);
	free(str);
	return hash;
}

/*
 * Get next alignment.
 */
static int
readaln(void *data, bam1_t *b)
{
	struct bycatch *byc = data;
	int ret;

	while (1) {
		ret = bam_read1(byc->bam_fp, b);
		if (ret < 0)
			break;

		if (b->core.flag & (BAM_FUNMAP | BAM_FQCFAIL | BAM_FDUP) ||
				b->core.qual < byc->min_qual)
			// skip these
			continue;

		break;
	}

	return ret;
}

static inline void
print_locus(locus_t *l, char **target_name)
{
	int tid = l->id >> 32;
	int pos = l->id & 0xffffffff;
	double shape = l->n ? l->m / l->n : 0.0;
//	double var = l->n ? l->s / l->n : 0.0;
//	double sd = sqrt(var);

	//printf("CHROM\tPOS\tarea\tnonzero sites\tsite depth\tmax depth\tmean\tshape\n");
	printf("%s\t%d\t%d\t%d\t%d\t%d\t%.3lf\t%.3lf\n",
			target_name[tid], pos, l->window_cov,
			l->n, l->site_cov, l->max_cov,
			l->m, shape);
}

/*
 * Record statistics based upon the depth of coverage for each site
 * within a window around a target site.
 */
int
bycatch(struct bycatch *byc)
{
	bam_plp_t iter;
	const bam_pileup1_t *plp;
	int i, tid, pos, n;

	klist_t(loci) *windows; // loci whose windows bound the current pos
	kliter_t(loci) *wi;
	locus_t *l, *u; // closest lower, upper targets to the pos
	unsigned int cov; // coverage for current pos

	byc->bam_fp = bam_open(byc->bam_file, "r");
	if (byc->bam_fp == NULL) {
		fprintf(stderr, "bam_open: %s: %s\n",
				byc->bam_file, strerror(errno));
		return -1;
	}

	byc->bam_hdr = bam_header_read(byc->bam_fp);
	byc->hi = loadpos(byc->bed_file, byc->bam_hdr);
	if (byc->hi == NULL)
		return -1;

	iter = bam_plp_init(readaln, byc);

	/*
	 * Iterate through the pileup.
	 */
	while ((plp = bam_plp_auto(iter, &tid, &pos, &n)) != 0) {

		/*
		 * Find the targets which have the current locus within
		 * their window, putting the target loci into the windows list.
		 */
		uint64_t x = ((uint64_t)tid)<<32 | max(1, pos - byc->window_size/2);
		locus_t lx = {x, 0};
		windows = kl_init(loci);
		kb_interval(loci, byc->hi, lx, &l, &u);

		while (u) {
			int tid2 = u->id >> 32;
			if (tid2 != tid)
				break;

			int pos2 = u->id & 0xffffffff;
			if (pos + byc->window_size/2 < pos2)
				break;

			*kl_pushp(loci, windows) = u;

			lx.id = u->id + 1;
			kb_interval(loci, byc->hi, lx, &l, &u);
		}

		if (windows->size == 0)
			// no target sites have this locus in their window
			goto next;

		cov = 0;
		for (i=0; i<n; i++) {
			const bam_pileup1_t *p = plp + i;

			if (p->is_del || p->is_refskip)
				continue;
			
			if (bam_get_qual(p->b)[p->qpos] < byc->min_qual)
				continue;

			cov++;
		}

		for (wi=kl_begin(windows); wi!=kl_end(windows); wi=kl_next(wi)) {
			l = kl_val(wi);
			l->window_cov += cov;
			if (cov) {
				l->n++;

				if (cov > l->max_cov)
					l->max_cov = cov;
				if (l->id == (((uint64_t)tid<<32) | pos))
					l->site_cov = cov;

				/*
				 * Maintain cumulative mean and variance.
				 * Knuth TAOCP vol 2, 3rd edition, page 232.
				 */
				if (l->n == 1) {
					l->m = cov;
					l->s = 0.0;
					l->md = 0.0;
				} else {
					double last_m = l->m;
					l->m += (cov - last_m)/l->n;
					l->s += (cov - l->m) * (cov - last_m);
					// median estimator - varies in accuracy
					l->md += copysign(l->m * 0.01, cov - l->md);
				}
			}
		}

next:
		kl_destroy(loci, windows);
	}

	bam_plp_destroy(iter);
	bam_close(byc->bam_fp);

	printf("CHROM\tPOS\tarea\tnonzero sites\ttarget depth\tmax depth\tmean\tshape\n");
#define traverse_f(p) print_locus((p), byc->bam_hdr->target_name)
	__kb_traverse(locus_t, byc->hi, traverse_f);

	bam_header_destroy(byc->bam_hdr);
	kb_destroy(loci, byc->hi);

	return 0;
}

void
usage(char *argv0)
{
	fprintf(stderr,
		"usage: %s [-q MINQUAL] [-w WINDOWSIZE] -l loci.bed in.bam\n",
		argv0);
}

int
main(int argc, char **argv)
{
	struct bycatch byc = {0,};
	int opt;

	// defaults
	byc.min_qual = 20;
	byc.window_size = 400;

	while ((opt = getopt(argc, argv, "l:q:w:")) != -1) {
		switch (opt) {
			case 'l':
				byc.bed_file = optarg;
				break;
			case 'q':
				errno = 0;
				byc.min_qual = strtoul(optarg, NULL, 0);
				if (errno || byc.min_qual > 255) {
					fprintf(stderr, "Error: expected MINQUAL in the range 0 to 255.\n");
					usage(argv[0]);
					return -1;
				}
				break;
			case 'w':
				errno = 0;
				byc.window_size = strtoul(optarg, NULL, 0);
				if (errno || byc.window_size > 1000000) {
					fprintf(stderr, "Error: expected WINDOWSIZE in the range 0 to 1000000.\n");
					usage(argv[0]);
					return -1;
				}
				break;
			default:
				usage(argv[0]);
				return -1;
		}
	}

	if (argc-optind != 1 || byc.bed_file == NULL) {
		usage(argv[0]);
		return -1;
	}

	byc.bam_file = argv[optind];

	return bycatch(&byc);
}

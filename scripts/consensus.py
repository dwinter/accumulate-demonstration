#!/usr/bin/env python

'''
consensus.py [input.vcf] > results.out

NOTE: This hard-codes the filtering steps used by the orginal authors
and assumes only three samples are included in the VCF file. Would need to be
re-written for new files
'''

import vcf
import sys

def minor_read_freq(s):
    if s.data.AD:
        if s.data.DP > 0:
            return(float(min(s.data.AD))/sum(s.data.AD))
    return(1.0)



def consensish(r):
    if r.num_called == 3:
        if 0 < r.num_hom_alt < 3:
            if sum( [s.data.DP > 3 for s in r.samples if s.data.DP]) == 3:
                return( all([minor_read_freq(s) < 0.1 for s in r.samples]))
    return(False)
            
def describe_mutant(r):
    if r.num_hom_alt == 1:
        mutant   = r.get_hom_alts()[0]
        ancestor = r.get_hom_refs()[0]
    else:        
        mutant   = r.get_hom_refs()[0]
        ancestor = r.get_hom_alts()[0]
    return( (r.CHROM, str(r.POS), mutant.sample, ancestor.gt_bases, mutant.gt_bases) )

def main():
    infile = sys.argv[1]
    recs = vcf.Reader(open(infile))
    consensus = (r for r in recs if consensish(r) )
    for r in consensus:
        sys.stdout.write( "\t".join(describe_mutant(r)) + "\n" )
    return(0)

if __name__ == "__main__":
    main()

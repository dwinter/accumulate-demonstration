#!/usr/bin/python
from Bio import SeqIO
from Bio.SeqUtils import GC
import sys

def GC_content(r):
    """calculates the GC content of a SeqRecord object (r) and multiples 
the percentage by the length of the object. returns the number of 
GC nucleotides per object. as each read is not identical in length, the 
GC content was calculated proportional to each object and added to 
create the GC content of the entire sequence file produced by Bio 
SeqIO.pars.
    """
    fraction=GC(r.seq)/100.0
    length=len(r.seq)
    GC_bases=fraction*length
    return(GC_bases)
      
    
def base_counter(ref):
    """counts the number of GC nucleotides in a SeqRecord object and the 
total number of nucleotides in a SeqRecord object
    """ 
    GC_bases=0
    bases=0
    for r in ref: 
        nbase = len(r.seq)
        GC_bases+=(GC(r.seq)/100)*nbase 
        bases+=nbase
    return(GC_bases, bases) 
    
 
def main():
    try:
        filename=sys.argv[1]
    except IndexError:
        print("Usage GC_content,py [reference genome (fasta)]")
        exit(1)
    ref=SeqIO.parse(filename,"fasta")
    GC_bases, bases=base_counter(ref)
    GC_percentage=GC_bases/bases
    AT_percentage=1-GC_percentage
    #print(GC_percentage) 
    print ("nfreqs={0:.3f} {1:.3f} {1:.3f} {0:.3f}".format(AT_percentage/2.0, GC_percentage/2.0))
   
  
if __name__ == "__main__":
    main()




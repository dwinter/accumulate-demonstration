#!/usr/bin/python
import sys
from Bio import SeqIO 

def main():
    """Represent a reference genome as a dictionary of chromosome names 
    and chromosome lengths for the purpose of using the dictionary to 
    create windows of non-overlapping genomic regions to run the software 
    accuMUlate in parallel. 
    Usage: 
        $python2 dictionary_converter.py [reference_genome.fasta]
    """  
    try: 
        file_input=sys.argv[1] 
    except IndexError: 
        print("Usage: dictionary_converter.py [file.fasta]")     
        sys.exit(1) 

    with open(sys.argv[1]) as file_input:
        records=SeqIO.parse(file_input, "fasta")
        fasta=[(x.id,len(x)) for x in records]
        if len(fasta)==0:
            raise RuntimeError("could not extract fasta records from file, check file format")
            sys.exit(1)
        for pair in fasta:
            print "{}\t{}".format(*pair)
        sys.exit(0)

    
if __name__ == "__main__":
    main()

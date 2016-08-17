#!/usr/bin/python
import sys
from collections import defaultdict

def parse_line(pairs):
    """Represent bam header key-value pairs as a dictionary

    Note: expects input to be a list header elements with the BAM 'key:value'
    format (e.g. ["VN:1.4", "GO:none","SO:coordinate"]). 
    """
    return( {k:v.strip() for k,v in [element.split(":", 1) for element in pairs]} )

def make_header_dict(handle):
    """Represent a complete BAM header as a dictoinary

    The returned dictionary contains keys for each unique header-type (i.e.
    thing that starts with an @ symbol) with each value being a list of
    dictionaries containing the key-value pairs defined in that header row.
    """
    res = defaultdict(list)
    for line in handle:
        elements = line.split("\t")
        name = elements[0].strip("@")
        res[name].append(parse_line(elements[1:]))
    return(res)

def print_usage():
    print ("""
    Usage, either of
        $ parse_rg.py [ancestral sample] [bam_header.txt]
        $ samtools view -H [alignment.bam] | parse_rg.py [ancestral sample] - 
    """)

def main():
    """ """
    try: 
        anc = sys.argv[1]
        fname = sys.argv[2]
    except IndexError:
        print_usage()
    if fname == "-":
        fhandle = sys.stdin
    else:
        fhandle = open(fname)
    d = make_header_dict(fhandle)
    samples =  [r["SM"] for r in d["RG"]]
    if anc in samples:        
        print("ancestor={}".format(anc))
    else:
        print("ERROR: Ancestral sample {} not found in header".format(anc))
        exit(1)
    for samp in set([s for s in samples if s != anc]):
        print( "sample-name={}".format(samp))
    exit(0)

if __name__ == "__main__":
    main()

#!/usr/bin/env python2

"""
Fetch all FASTQ files associated wth an SRA ID

usage:

    sra_fetch.py [sra-id]

"""
import urllib2
import sys
import os
import time
import re
import argparse
import subprocess 

import xml.etree.ElementTree as ET

is_num = re.compile("\\d+")
is_char = re.compile("[A-Z]+")
base_xml_url = "http://www.ebi.ac.uk/ena/data/view/{0}&display=xml"
base_fastq_url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{0}/{1}/{2}"
CHUNK = 16 * 1024


def extract_ranges(run_id_string):
    """Convert ranges of SRA run IDs into lists
    
    e.g. "SRXX121-SRXX123" -> ["SRXX121", "SRXX122", "SRXX123"]
    """
    un_ranged = run_id_string.split("-")
    if len(un_ranged) == 1:
        return(un_ranged)
    if len(un_ranged) > 2:
        raise ValueError (run_id_string)
    guff = is_char.findall(un_ranged[0])[0] + "{0}" 
    ids = [int(is_num.findall(id)[0]) for id in un_ranged]
    return( [guff.format(x) for x in range(ids[0], ids[1]+1)] )

def unnest(L):
    """ """
    return [item for sublist in L for item in sublist]

def fetch_run_ids(SRA_ID):
    """Find run ID associated with experiment ID
    
    Retuns a list of run ids as a boolean specifying whether the reads
    are paired-end.
    """
    rec = ET.parse(urllib2.urlopen(base_xml_url.format(SRA_ID)))
    ena_node = [r for r in rec.findall(".//XREF_LINK") if r[0].text == "ENA-RUN"][0]
    nested = [extract_ranges(x) for x in ena_node[1].text.split(",")]
    return unnest(nested), bool( rec.findall(".//PAIRED") )

def fetch_fastq(run_id, paired=True, overwrite=False):
    if paired:
        files = ["{0}_{1}.fastq.gz".format(run_id,i) for i in [1,2]]
    else:
        files = [run_id + ".fastq.gz"]
    for f in files:
        fname = "fastq/" + f
        if os.path.exists(fname):
            print("Already there..", overwrite)
            if not overwrite:
                continue
        url = base_fastq_url.format(run_id[:6], run_id, f)
        print "Downloading " + f + "..."
        response = urllib2.urlopen(url)
        with open(fname, "wb") as out:
            while True:
                next_chunk = response.read(CHUNK)
                if not next_chunk:
                    response.close()
                    time.sleep(6)#let's not DoS the ENA
                    break
                out.write(next_chunk)

#TODO: actually handle paried data
def make_run_alignment(exp_id, run_id, is_paired, ref_genome, nproc, overwrite=False):
    """ """
    if not overwrite:
        if os.path.exists("bam/{}.bam".format(run_id)):
            sys.stderr.write("bam already exists, skipping")
            return(0)
    RG = '"@RG\tID:{0}\tSM:{1}\tPL:illumina\tLB:{1}"'.format(run_id, exp_id)
    if is_paired:
        fastq_input = "fastq/{0}_1.fastq.gz fastq/{0}_2.fastq.gz".format(run_id)
    else:
        fastq_input = "fastq/{0}.fastq.gz".format(run_id) 
    cmd = "bwa mem -t {0} -M {1} {2} -R {3} | samtools sort -@ {0} -T /tmp/aln.sorted -o bam/{4}.bam".format(nproc, ref_genome, fastq_input, RG, run_id)
    print(cmd)
    ret_code = subprocess.call(cmd, shell=True)
    return(ret_code)

def parse_args():
    parser = argparse.ArgumentParser(prog="SRA_to_BAM.py")
    parser.add_argument("id",  help="SRA experiment id")
    parser.add_argument("ref", help="Reference genome")
    parser.add_argument("nproc", help="Number of processors")
    parser.add_argument("--overwrite", 
                        help="Overwrite file if it exists in fastq/", 
                        action="store_true")
    return(parser.parse_args())
    


def main():
    arg_vals = parse_args()
    ids, is_paired = fetch_run_ids(arg_vals.id)
    for run_id in ids:
        fetch_fastq(run_id, is_paired, arg_vals.overwrite)
        make_run_alignment(arg_vals.id, run_id, is_paired, arg_vals.ref, arg_vals.nproc, arg_vals.overwrite)
    return(0)



if __name__ == "__main__":
    main()

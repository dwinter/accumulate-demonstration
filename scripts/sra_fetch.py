#!/usr/bin/env python2

"""
Fetch all FASTQ files associated wth an SRA ID

usage:

    sra_fetch.py [sra-id]

"""
import xml.etree.ElementTree as ET
import urllib2
import sys
import time
import re

is_num = re.compile("\\d+")
is_char = re.compile("[A-Z]+")
base_xml_url = "http://www.ebi.ac.uk/ena/data/view/{0}&display=xml"
base_fastq_url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{0}/{1}/{2}"
CHUNK = 16 * 1024


def extract_ranges(run_id_string):
    """ """
    un_ranged = run_id_string.split("-")
    if len(un_ranged) == 1:
        return(un_ranged)
    if len(un_ranged) > 2:
        raise ValueError (run_id_string)
    guff = [is_char.findall(id)[0] + "{0}" for id in un_ranged]
    ids = [int(is_num.findall(id)[0]) for id in un_ranged]
    return( [guff[i].format(x) for i,x in enumerate(range(ids[0], ids[1]+1))] )

def unnest(L):
    """ """
    return [item for sublist in L for item in sublist]

def fetch_run_ids(SRA_ID):
    """ """
    rec = ET.parse(urllib2.urlopen(base_xml_url.format(SRA_ID)))
    ena_node = [r for r in rec.findall(".//XREF_LINK") if r[0].text == "ENA-RUN"][0]
    nested = [extract_ranges(x) for x in ena_node[1].text.split(",")]
    return unnest(nested), bool( rec.findall(".//PAIRED") )

def fetch_fastq(run_id, paired=True):
    if paired:
        files = ["{0}_{1}.fastq.gz".format(run_id,i) for i in [1,2]]
    else:
        files = [run_id + ".fastq.gz"]
    for f in files:
        url = base_fastq_url.format(run_id[:6], run_id, f)
        print "Downloading " + f + "..."
        response = urllib2.urlopen(url)
        with open(f, "wb") as out:
            while True:
                next_chunk = response.read(CHUNK)
                if not next_chunk:
                    response.close()
                    time.sleep(6)#let's not DoS the ENA
                    break
                out.write(next_chunk)

def main():
    SRA_ID = sys.argv[1]
    ids, is_paired = fetch_run_ids(SRA_ID)
    for sid in ids:
        fetch_fastq(sid, is_paired)

if __name__ == "__main__":
    main()

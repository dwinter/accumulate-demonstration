SHELL := /bin/bash
ref_genome := ref/Dicty.fasta
nproc := 20
basename := dicty


downloaded:
	mkdir -p bam
	mkdir -p fastq
	scripts/SRA_to_BAM.py SRX174012 $(ref_genome) $(nproc)
	scripts/SRA_to_BAM.py SRX174013 $(ref_genome) $(nproc)
	scripts/SRA_to_BAM.py SRX174014 $(ref_genome) $(nproc)
	touch downloaded

bam/$(basenmae).bam: downloaded
	$(MAKE) remove-intermediates #be sure some stale versions aren't getting in the way
	samtools merge -r bam/merge.bam  bam/*.bam
	# the PG tags in this merge include the RG tag used in the bwa mem command line
	# this is bad, because many parsers chose on the fact there are two IDs.
	# Here's a hack to remove the RG tags from the call
	samtools view -H bam/merge.bam | grep "^@PG" | cut -f 1-5 > PG.txt
	samtools view -H bam/merge.bam | grep -v "^@PG" | cat - PG.txt > header.txt
	echo -e "@RG\tID:ANCESTOR\tPL:illumina\tLB:XXXX\tSM:A0" >> header.txt
	samtools reheader -P header.txt bam/merge.bam > 
bam/$(basename).bam
	samtools index bam/$(basename).bam
	rm header.txt
	rm bam/merge.bam
	rm PG.txt

.PHONY: alignment
alignment:
	$(MAKE) bam/$(basename).bam


##Summary

stats/$(basename).idxstats: bam/realigned.bam
	samtools idxstats bam/realigned.bam > stats/$(basename).idxstats

stats/$(basename).metrics: bam/realigned.bam
	picard CollectWgsMetrics R=$(ref_genome) I=bam/realigned.bam 
O=stats/$(basename).metrics

stats/$(basename)_cov: bam/realigned.bam
	gatk -T DepthOfCoverage -R $(ref_genome) -I bam/realigned.bam -o 
stats/$(basename)_cov -nt $(nproc) --omitIntervalStatistics 


.PHONY: summary
summary:
	$(MAKE) stats/$(basename).idxstats
	$(MAKE) stats/$(basename).metrics

##Clean

bam/md.bam: bam/$(basename).bam
	mkdir -p stats
	mkdir -p logs 
	picard MarkDuplicates INPUT=bam/$(basename).bam 
OUTPUT=bam/md.bam  
METRICS_FILE=stats/dup_mark_metrics.txt ASSUME_SORTED=true 2> logs/dup_marking.log

indels.intervals: bam/md.bam
	samtools index bam/md.bam
	gatk -T RealignerTargetCreator -I bam/md.bam -R $(ref_genome) 
-nt $(nproc) -o indels.intervals 

bam/realigned.bam: bam/md.bam indels.intervals
	gatk -T IndelRealigner -I bam/md.bam -R $(ref_genome) 
-targetIntervals indels.intervals -o bam/realigned.bam
	samtools index bam/realigned.bam 

.PHONY: clean_alignment
clean_alignment:
	$(MAKE) bam/realigned.bam


##Analysis

vars/gatk_haploid_raw.vcf: bam/realigned.bam
	mkdir -vars 
	gatk -T HaplotypeCaller -R $(ref_genome) -I bam/realigned.bam 
-ploidy 1 -o 
vars/gatk_haploid_raw.vcf 2> logs/hc.log

results/consensus.out: vars/gatk_haploid_raw.vcf
	scripts/consensus.py vars/gatk_haploid_raw.vcf > results/consensus.out

##Parameters

parameters.ini: bam/realigned.bam 
	cp template_parameters.ini parameters.ini 
	samtools view -H bam/realigned.bam | python scripts/extract_samples.py A0 - >> parameters.ini
	python scripts/GC_content.py $(ref_genome) >> parameters.ini


.PHONY: parameters
parameters:
	$(MAKE) parameters.ini

##Results 
#results/raw_accuMUlate.out: parameters.ini 
#	./accuMUlate -c parameters.ini -b bam/realigned.bam -r ref/dicty_chromosomal.fasta >results/raw_accuMUlate.out


##Dictionary

dictionary_reference_genome.bed: $(ref_genome)
	python2 scripts/dictionary_converter.py $(ref_genome) | grep -v 
"^DDB0169550" > dictionary_reference_genome.bed

tmp/aa: dictionary_reference_genome.bed 
	mkdir -p tmp
	bedtools makewindows -g dictionary_reference_genome.bed -w 10000000 | split -l 1 - tmp/

.PHONY: dictionary
dictionary:
	$(MAKE) tmp 


##Parallel

results/raw_accuMUlate.out: tmp parameters.ini bam/realigned.bam 
	mkdir -p results
	rm -f results/unsorted.out 
	parallel -j $(nproc) ./accuMUlate -c parameters.ini -b 
bam/realigned.bam -r $(ref_genome) -i {} ::: tmp/* > 
results/unsorted.out
	sort -k1,1 -k2,2n results/unsorted.out >results/raw_accuMUlate.out

.PHONY: result
result:
	$(MAKE) results/raw_accuMUlate.out

.PHONY: analysis
analysis:
	$(MAKE) vars/gatk_haploid_raw.vcf
	$(MAKE) results/consensus.out

.PHONY: remove-intermediates
remove-intermediates:
	rm -f bam/dicty.bam
	rm -f bam/dicty.bam.bai
	rm -f bam/md.bam
	rm -f bam/md.bam.bai
	rm -f bam/merge.bam

.PHONY: clean
clean:
	rm -f indels.intervals
	rm -f vars/*
	rm -f stats/*
	rm -f logs/*
	$(MAKE) remove-intermediates 

.PHONY: all
all:
	$(MAKE) analysis
	$(MAKE) summary
	$(MAKE) remove-intermediates 

#bam/SRR534844.bam: fastq/SRR534844.fastq.gz ref/dicty_chromosomal.fasta
#	bwa mem -t 20 -M ref/dicty_chromosomal.fasta fastq/SRR534844.fastq.gz -R "@RG\tID:SRR534844\tSM:SRR534844\tPL:illumina\tLB:SRR534844" | samtools sort -T /tmp/aln.sorted -o bam/SRR534844.bam
#
#bam/SRR534465.bam: fastq/SRR534465.fastq.gz ref/dicty_chromosomal.fasta
#	bwa mem -t 20 -M ref/dicty_chromosomal.fasta fastq/SRR534465.fastq.gz -R "@RG\tID:SRR534465\tSM:SRR534465\tPL:illumina\tLB:SRR534465"	| samtools sort -T /tmp/aln.sorted -o bam/SRR534465.bam
#
#bam/SRR534464.bam: fastq/SRR534464.fastq.gz ref/dicty_chromosomal.fasta
#	bwa mem -t 20 -M ref/dicty_chromosomal.fasta fastq/SRR534464.fastq.gz -R "@RG\tID:SRR534464\tSM:SRR534464\tPL:illumina\tLB:SRR534464" | samtools sort -T /tmp/aln.sorted -o bam/SRR534464.bam
#

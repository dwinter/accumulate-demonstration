SHELL := /bin/bash
ref_genome := ref/Athal.fasta
nproc := 1
basename := athal

$(ref_genome).bwt:
	mkdir -p ref
	wget  https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas -O ref/temp.fasta
	awk '{print $1}' ref/temp.fasta  > $(ref_genome)
	rm ref/temp.fasta
	samtools faidx $(ref_genome)
	picard-tools CreateSequenceDictionary R=$(ref_genome) O=$(ref_genome:.fasta=.dict)
	bwa index -a bwtsw $(ref_genome)

.PHONY: reference_genome
reference_genome:
	$(MAKE) $(ref_genome).bwt

downloaded: $(ref_genome).bwt
	mkdir -p bam
	mkdir -p fastq
	scripts/SRA_to_BAM.py ERX386699 $(ref_genome) $(nproc) #59
	scripts/SRA_to_BAM.py ERX386701 $(ref_genome) $(nproc) #79
	scripts/SRA_to_BAM.py ERX386702 $(ref_genome) $(nproc) #89
	scripts/SRA_to_BAM.py ERX386703 $(ref_genome) $(nproc) #99
	scripts/SRA_to_BAM.py ERX386705 $(ref_genome) $(nproc) #119
	touch downloaded
#
bam/$(basename).bam: downloaded
	$(MAKE) remove-intermediates #be sure some stale versions aren't getting in the way
	samtools merge -r bam/merge.bam  bam/*.bam
	# the PG tags in this merge include the RG tag used in the bwa mem command line
	# this is bad, because many parsers chose on the fact there are two IDs.
	# Here's a hack to remove the RG tags from the call
	samtools view -H bam/merge.bam | grep "^@PG" | cut -f 1-5 > PG.txt
	samtools view -H bam/merge.bam | grep -v "^@PG" | cat - PG.txt > header.txt
	#add a fake ancestor
	echo -e "@RG\tID:ANCESTOR\tPL:illumina\tLB:XXXX\tSM:A0" >> header.txt 
	samtools reheader -P header.txt bam/merge.bam >bam/$(basename).bam
	samtools index bam/$(basename).bam
	rm header.txt
	rm bam/merge.bam
	rm PG.txt

.PHONY: alignment
alignment:
	$(MAKE) bam/$(basename).bam

#
##Summary
#
stats/$(basename).idxstats: bam/realigned.bai
	samtools idxstats bam/realigned.bam > stats/$(basename).idxstats

stats/$(basename).metrics: bam/realigned.bai
	picard-tools CollectWgsMetrics R=$(ref_genome) I=bam/realigned.bam O=stats/$(basename).metrics

stats/$(basename)_cov: bam/realigned.bai
	gatk -T DepthOfCoverage -R $(ref_genome) -I bam/realigned.bam -o stats/$(basename)_cov -nt $(nproc) --omitIntervalStatistics
#

.PHONY: summary
summary:
	$(MAKE) stats/$(basename).idxstats
	$(MAKE) stats/$(basename).metrics
	$(MAKE) stats/$(basename)_cov

#
###Clean the alignment up
#

bam/md.bam: bam/$(basename).bam
	mkdir -p stats
	mkdir -p logs
	picard-tools MarkDuplicates INPUT=bam/$(basename).bam OUTPUT=bam/md.bam  METRICS_FILE=stats/dup_mark_metrics.txt ASSUME_SORTED=true 2> logs/dup_marking.log
#
indels.intervals: bam/md.bam
	samtools index bam/md.bam
	gatk -T RealignerTargetCreator -I bam/md.bam -R $(ref_genome) -nt $(nproc) -o indels.intervals 
#
bam/realigned.bai: bam/md.bam indels.intervals
	gatk -T IndelRealigner -I bam/md.bam -R $(ref_genome) -targetIntervals indels.intervals -o bam/realigned.bam

.PHONY: clean_alignment
clean_alignment:
	$(MAKE) bam/realigned.bam

###Analysis
#
# Might as well do our usual trick of looking for variants with gatk
vars/gatk_haploid_raw.vcf: bam/realigned.bai
	mkdir -p vars
	gatk -T HaplotypeCaller -R $(ref_genome) -I bam/realigned.bam -ploidy 2 -o vars/gatk_haploid_raw.vcf 2> logs/hc.log

Athal.genome: 
	#only look at the nuclear genome, not the mt or plastid
	python2 scripts/dictionary_converter.py ref/Athal.fasta | egrep ^[1-5] > Athal.genome

windows/: Athal.genome
	mkdir -p windows
	bedtools makewindows -g Athal.genome -w 1000000  | split -l 1 - tmp/

params.ini: bam/realigned.bai
	cp params_template.ini params.ini
	samtools view -H bam/realigned.bam | python2 scripts/extract_samples.py A0 - >> params.ini
	python2 scripts/GC_content.py $(ref_genome) >> params.ini

	
results/accu_raw.out: bam/realigned.bai params.ini windows/
	mkdir -p results
	rm -f results/acc_raw_unsorted.out 	
	parallel -j $(nproc) accuMUlate -c params.ini -b bam/realigned.bam -x bam/realigned.bai -r $(ref_genome) -i {} -m30 '>>' results/acc_raw_unsorted.out  ::: tmp/* 
	sort -k1,1 -k2,2n results/acc_raw_unsorted.out  > results/accu_raw.out


random_intervals/: Athal.genome
	mkdir -p random_intervals/
	bedtools random -n 6000 -l 1000 -seed 123321 -g Athal.genome  | split -l 1 - random_intervals/

results/denom.out: bam/realigned.bai random_intervals/
	parallel -j $(nproc) denominate -b bam/realigned.bam -x bam/realigned.bai -r $(ref_genome) -i {} -m30 -c denom_params.ini  '>>' results/denom.out ::: random_intervals/*

athal_analysis.pdf: results/accu_raw.out results/denom.out athal_analysis.Rmd
	Rscript -e 'rmarkdown::render("athal_analysis.Rmd")'
	

.PHONY: analysis
analysis:
	$(MAKE) athal_analysis.pdf

.PHONY: remove-intermediates
remove-intermediates:
	rm -f bam/$(basename).bam
	rm -f bam/$(basename).bam.bai
	rm -f bam/md.bam
	rm -f bam/md.bam.bai
	rm -f bam/merge.bam
	rm -f header.txt
	rm -f PG.txt

.PHONY: all
all:
	$(MAKE) analysis
	$(MAKE) summary
	$(MAKE) remove-intermediates 


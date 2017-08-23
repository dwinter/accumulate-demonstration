is_het <- function(g) substr(g,1,1) != substr(g, 2,2)


#read in accumulate results files
read_accu <- function(file, ...){
    mu <- read.table(file, stringsAsFactors=FALSE, ...)
    names(mu) <- c("chrom", "from", "to", "ref", "mutant_sample", "mutation", 
                   "p_any", "p_one", "p_geno", "likelihood", "DP", 
                   "mutant_F", "mutant_R", "anc_in_mutant", "mutant_in_anc",
                   "MQ_diff", "insert_diff", "strand_pval", "pair_pval")
    mu
}


#read in denominate output
# If using a sampel of the genome 
# genome_size = the total number of base pairs used for calling mutations, 
# sample_size = number of base pairs used for denominator estimation

read_denom <- function(denom_file, params_file, genome_size=1, sample_size=1){
    scale <- genome_size/sample_size
    denom <- read.table(denom_file, stringsAsFactors=FALSE)
    denom_params <- read.table(params_file, 
                           sep="=", 
                           col.names=c("param", "value"), 
                           stringsAsFactors=FALSE)
    denom_by_base <- t(sapply(1:4, function(i) colSums(denom[,seq(i,ncol(denom), 4)]))) * scale
    colnames(denom_by_base) <- denom_params[denom_params$param == "sample-name","value"]
    rownames(denom_by_base) <- c("A", "C", "G", "T")
    denom_by_base
}

tidy_denom <- function(denom_table) melt(t(denom_table), value.name="n", varnames=c("sample", "base"))

summarize_denom <- function(denom_table){
    res <- data.frame(
        total = colSums(denom_table),
        AT = colSums(denom_table[c(1,4),]),
        GC = colSums(denom_table[2:3,])
     )
    res
}

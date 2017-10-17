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

read_denom <- function(denom_file, params_file, genome_size=1, interval_len = NULL){
    denom <- read.table(denom_file, stringsAsFactors=FALSE)
    if(interval_len){
        sample_size <- interval_len * nrow(denom)        
        scale <- genome_size/sample_size
    }else{
        scale <- 1
    }
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
    res_df <- data.frame(
        sample = colnames(denom_by_base),                     
        total = colSums(denom_table),
        AT = colSums(denom_table[c(1,4),]),
        GC = colSums(denom_table[2:3,])
     )
    colnames(res_df) <- NULL
    res_df    
}


mutation_specturm <- function(mu_table, denom_table, ngen){
    mu_by_type <- as.data.frame(table(mu_table$six_mus, dnn="mutation_type"))
    denom_summ <- summarize_denom(denom_table)
    denom_by_base <- rep(colSums(denom_summ[,2:3]), each=3)
    mu_by_type$rate <- mu_by_type$Freq/denom_by_base
    mu_by_type
}

.spectrum_denom <- function(denom_by_base, ngen){
    nsamp <- ncol(denom_by_base)
    snames <- colnames(denom_by_base)
    at_by_sample <- colSums(denom_by_base[c(1,4),] ) * ngen
    cg_by_sample <- colSums(denom_by_base[2:3,] ) * ngen
    res <- lapply(snames, 
                  function(x) rep(c(at_by_sample[x], cg_by_sample[x]),each=3))
    names(res) <- snames
    res
}

rate_table <- function(mu, denom_by_base, ngen, total=TRUE){
    snames <- colnames(denom_by_base)
    sample_spectrum <- aggregate(as.factor(six_mus) ~ mutant_sample, 
                                 data=mu, FUN=table, simplify=FALSE)
    spec_list <- structure(sample_spectrum[[2]], .Names=sample_spectrum$mutant_sample)
    D <- .spectrum_denom(denom_by_base, ngen)
    res <- as.data.frame(sapply(snames, function(x) spec_list[[x]]/D[[x]]))
    if(total){
        total_denom <- rep(rowSums(sapply(D, "[", c(1,4))), each=3)
        total_n <- rowSums(sapply(spec_list, "+"))
        res$total <- total_n/total_denom
    }
    data.frame(t(res))
}


filtering_table <- function(filtered, unfiltered, denom_by_base, ngen){
    denom_by_sample <- colSums(denom_by_base) * ngen
    sample_names <- names(denom_by_sample)
    rate_by_sample <- round(mu_by_sample/denom_by_sample,12)[sample_names]
    res <- data.frame(
        sample = c(sample_names, "All hets", "All"),
        ngen = c(ngen, rep(sum(ngen),2)),
        n.initial = c(table(mu$mutant_sample)[sample_names],sum(mu$desc_het), nrow(mu)),
        n.filter = c(mu_by_sample[sample_names], sum(snm$desc_het), nrow(snm)),
        rate = c(rate_by_sample, "-", round(nrow(snm) / sum(denom_by_sample), 12))
    )
    res
}

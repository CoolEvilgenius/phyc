# Output files
set_output <- function(out_dir, files) {
    return(sapply(1:length(files), function(x) paste(out_dir, files[x], sep='')))
}

write_df <- function (df, file_name) {
    write.table(df, file = file_name, quote = FALSE, sep = "\t", 
                row.names = FALSE, col.names = FALSE)
}

filter_results <- function(df, stats, MAF) {
    names(df)[2] <- "name"
    df$freq <- stats$freq
    df <- df[df$freq >= MAF,]
    
    return(df)
}

#
#   Plots
#
manhattan <- function(df, bc, genome_size) {
    df$x <- sapply(1:nrow(df), 
                   function(x) {as.integer(substring(df$name[x], 2)) })
    
    p <- ggplot(df, aes(x=x, y=logp)) +
        xlab("Position (Mb) along N315 chromosome") + 
        ylab(expression(-log[10](P))) +
        geom_point() +
        scale_x_continuous(breaks = round(seq(0, genome_size, by = 250000),1), 
                           labels = round(seq(0, genome_size, 
                                              by = 250000)/1000000,1)) +
        geom_hline(yintercept=bc, linetype="dashed") +
        theme_bw() +
        theme(text = element_text(family="Arial", size=20))

    return(p)
}

qqplot <- function(df, bc) {
    df$x <- sort(-1*log10(1:nrow(df) / nrow(df)))
    
    p <- ggplot(df, aes(x=x, y=logp)) +
        xlab(expression(Expected~~-log[10](P))) + 
        ylab(expression(Observed~~-log[10](P))) +
        geom_point() +
        geom_line() +
        geom_abline(intercept = 0, slope = 1) +
        geom_hline(yintercept=bc, linetype="dashed") +
        theme_bw() +
        theme(text = element_text(family="Arial", size=16))

    return(p)     
}


#
#   ROADTRIPS related
#
write_roadtrips <- function(snps, mic, threshold, rt_dir) {
    
    phenotype <- mic_to_phenotype(mic, threshold)
    pedigree <- phenotype_to_pedigree(phenotype)
    genotype <- t(snps[,2:ncol(snps)])
    snpnames <- data.frame(snpnames=colnames(snps)[2:ncol(snps)])
    prevalence <- data.frame(
        prevalence=length(phenotype$status[phenotype$status == 2])/
            nrow(phenotype))
    
    # write files
    output <- set_output(rt_dir, c("phenofile", "pedinfo", "genofile",
                                   "prevalence","snpnames"))
    system(paste("mkdir -p ", rt_dir, sep=''))
    write_df(phenotype, output[1])
    write_df(pedigree, output[2])
    write_df(genotype, output[3])
    write_df(prevalence, output[4])
    write_df(snpnames, output[5])
}

mic_to_phenotype <- function (mic, threshold) {
    # PAP-AUC not tested on 120338, but is known VISA
    if (is.na(mic[mic$strain == "120338", ]$mic)) {
        mic[mic$strain == "120338", ]$mic <- threshold
    }
    
    # satus -> 1:supsceptible, 2:resistant
    df <- data.frame(family_id=1:nrow(mic),
                     individual_id=rep(1, nrow(mic)), 
                     status=ifelse(mic$mic >= threshold, 2, 1))

    return(df)
}

phenotype_to_pedigree <- function (phenotype) {
    
    df <- data.frame(family_id=phenotype$family_id,
                     individual_id1=phenotype$individual_id,
                     individual_id2=phenotype$individual_id,
                     pedigree=rep(0,nrow(phenotype)))
    return(df)
}

run_roadtrips <- function(rt_dir, rt_path = "ROADTRIPS") {
    system(paste("cd ", rt_dir, " && ",rt_path," -n snpnames > roadtrips.out", sep=''))
    df <- read.table(paste(rt_dir, "ROADTRIPStest.pvalues", sep=''), header=TRUE)
    df$logp <- -1 * log10(df$RM)
    
    return(df)
}


#
# QROADTRIPS
#
write_qtrips <- function(snps, mic, qrt_dir) {
    phenotype <- qtrips_phenotype(mic)
    genotype <- qtrips_genotype(t(snps))
    
    
    # write files
    output <- set_output(qrt_dir, c("phenofile", "genofile"))
    system(paste("mkdir -p ", qrt_dir, sep=''))
    write_df(phenotype, output[1])
    write_df(genotype, output[2])
}

qtrips_phenotype <- function (mic) {
    # satus -> 1:supsceptible, 2:resistant
    df <- data.frame(family_id=1:nrow(mic),
                     individual_id=rep(1, nrow(mic)), 
                     father_id=rep(0, nrow(mic)), 
                     mother_id=rep(0, nrow(mic)), 
                     sex_id=rep(0, nrow(mic)),
                     status=mic$mic)
    return(df)
}

qtrips_genotype <- function(snps) {
    rows <- nrow(snps)
    tmp <- rep(0,nrow(snps)-1)
    df <- data.frame(chr=tmp, 
                     snps=rownames(snps)[2:nrow(snps)], 
                     dist=tmp,
                     pos=tmp)
    
    for (i in 1:ncol(snps)) {
        tmp = as.integer(snps[2:nrow(snps),i])+1
        df[paste(snps[1,i], "_f", sep="")] <- tmp
        df[paste(snps[1,i], "_m", sep="")] <- tmp
    }
    
    return(df)
}

run_qtrips <- function(qrt_dir) {
    system(paste("cd ", qrt_dir," && qroadtrips -p phenofile -g genofile",
                 " > roadtrips.out", sep=""))
    
    df <- read.table(paste(qrt_dir, "QROADTRIPStest.pvalues", sep=''), 
                     header=TRUE)
    names(df)[3] <- "P"
    df2 <- read.table(paste(qrt_dir, "QROADTRIPStest.testvalues", sep=''), 
                     header=TRUE)
    df$QR <- df2$QR
    df$logp <- -1 * log10(df$P)
    
    return(df)
}


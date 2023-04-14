library(TwoSampleMR)
library(ieugwasr)
library(tidyverse)
library(LDlinkR)
library(MRPRESSO)
library(argparser)
library(glue)


###### get arguments from command line #############
p <- arg_parser("Perform mendelian randomization analysis using TwoSampleMR")
p <- add_argument(p, "--exp", help = "exposure gwas summary")
p <- add_argument(p, "--otc", help = "outcome gwas summary")
p <- add_argument(p, "--nexp", help = "sample size of exposure")
p <- add_argument(p, "--notc", help = "sample size of outcome")
p <- add_argument(p, "--prefix", help = "output prefix for results", default = "TMR_results")
p <- add_argument(p, "--expname", help = "exposure name", default = "exposure")
p <- add_argument(p, "--otcname", help = "outcome name", default = "outcome")
p <- add_argument(p, "--proxy", help = "use proxy snps", flag = TRUE)
p <- add_argument(p, "--online", help = "search proxy online", flag =TRUE)
p <- add_argument(p, "--pop", help = "specify LD population <default = EUR>", default = "EUR")
p <- add_argument(p, "--rsq", help = "specify LD r-square <default = 0.8>", default = 0.8)
p <- add_argument(p, "--bfile", help = "plink format file to use LD reference panel", default = "")
p <- add_argument(p, "--plink", help = "path to plink", default = "")
argv <- parse_args(p)

######## function for get proxy #########
get_proxy <- function(exp_dat, otc_dat, r2 = 0.8, pop = "EUR", online = TRUE, bin_plink, ref_panel, outfile = "proxy") {
    exp.snp = exp_dat$SNP
    otc.snp = otc_dat$SNP
    otc_dat = otc_dat %>% mutate(
        target_snp.outcome = NA,
        proxy_snp.outcome = NA,
        target_a1.outcome = NA,
        target_a2.outcome = NA,
        proxy_a1.outcome = NA,
        proxy_a2.outcome = NA,
        proxy.outcome = NA
        )
    miss.expsnps = setdiff(exp.snp, intersect(exp.snp, otc.snp))
    remain.otcsnps = setdiff(otc.snp, intersect(exp.snp, otc.snp))
    for (miss.snp in miss.expsnps) {
        if (online) {
            proxies.df = LDproxy(snp = miss.snp, pop = pop, token = "6fb632e022ef")
            top_proxies.df = proxies.df %>% filter(R2 >= r2, Distance > 0, RS_Number %in% remain.otcsnps) %>% arrange(-R2) %>% head(1)  
        } else {
            system(glue("{bin_plink} -bfile {ref_panel} --r2 --ld-snp {miss.snp} --ld-window 1000000 --ld-window-kb 1000 --ld-window-r2 {r2} --out {outfile}_{miss.snp}_proxy"))
            proxies.df = read_table(glue("{outfile}_{miss.snp}_proxy.ld"))
            top_proxies.df = proxies.df %>% filter(SNP_A != SNP_B, SNP_B %in% remain.otcsnps) %>% arrange(-R2) %>% head(1) %>% rename(RS_Number = SNP_B)
        }
        if (nrow(top_proxies.df) > 0) {
            cat("Found proxies of ", miss.snp,"\n", sep = "")
            rsid = top_proxies.df$RS_Number
            effect_allele.exposure = (exp_dat %>% filter(SNP == miss.snp) %>% head(1))$effect_allele.exposure
            other_allele.exposure = (exp_dat %>% filter(SNP == miss.snp) %>% head(1))$other_allele.exposure
            if (online) {
                proxy.allele1 = (str_replace_all((top_proxies.df %>% head(1))$Alleles, "[()]", "") %>% str_split("/"))[[1]][1]
                proxy.allele2 = (str_replace_all((top_proxies.df %>% head(1))$Alleles, "[()]", "") %>% str_split("/"))[[1]][2]
                target.allele1 = (((top_proxies.df %>% head(1))$Correlated_Alleles  %>% str_split(","))[[1]][1] %>% str_split("="))[[1]][1]
                target.allele2 = (((top_proxies.df %>% head(1))$Correlated_Alleles  %>% str_split(","))[[1]][1] %>% str_split("="))[[1]][2]
            } else {
                allele_info <- system(glue('{bin_plink} -bfile {ref_panel} --ld {miss.snp} {rsid} | grep \'In phase\' | awk \'{{print substr($NF,1,1)"\\n"substr($NF,4,1)"\\n"substr($NF,2,1)"\\n"substr($NF,5,1)}}\''), intern = TRUE)
                target.allele1 = allele_info[1]
                target.allele2 = allele_info[2]
                proxy.allele1 = allele_info[3]
                proxy.allele2 = allele_info[4]
            }

            otc_dat[otc_dat$SNP == rsid,]$target_snp.outcome = miss.snp
            otc_dat[otc_dat$SNP == rsid,]$proxy_snp.outcome = rsid
            otc_dat[otc_dat$SNP == rsid,]$target_a1.outcome = target.allele1
            otc_dat[otc_dat$SNP == rsid,]$target_a2.outcome = target.allele2
            otc_dat[otc_dat$SNP == rsid,]$proxy_a1.outcome = proxy.allele1
            otc_dat[otc_dat$SNP == rsid,]$proxy_a2.outcome = proxy.allele2
            otc_dat[otc_dat$SNP == rsid,]$proxy.outcome = TRUE
            otc_dat[otc_dat$SNP == rsid,]$effect_allele.outcome = effect_allele.exposure
            otc_dat[otc_dat$SNP == rsid,]$other_allele.outcome = other_allele.exposure
            otc_dat[otc_dat$SNP == rsid,]$SNP = miss.snp
        } else {
            cat("Not found proxies of ", miss.snp, "\n", sep = "")
        }
    }
    return(otc_dat %>% filter(SNP %in% exp.snp))
}
########## perform MR #####################
exp_gwas <- read_tsv(argv$exp)
exp_dat <- format_data(exp_gwas %>% filter(p_value <= 5e-8) %>% mutate(phenotype = argv$expname, id = argv$expname),
    snp_col = "variant_id",
    beta_col = "beta",
    se_col = "standard_error",
    chr_col = "chromosome",
    pos_col = "base_pair_location",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency",
    pval_col = "p_value",
    phenotype_col = "phenotype",
    id_col = "id",
    type = "exposure")

# exp_dat <- clump_data(exp_dat)
rsid_pval_exp <- exp_dat %>% select(SNP, pval.exposure) %>% rename(rsid = SNP, pval = pval.exposure)
rsid_pval_exp.clump <- ld_clump_local(rsid_pval_exp , clump_kb = 10000, clump_r2 = 0.001, clump_p = 1, plink_bin = argv$plink, bfile = argv$bfile)
exp_dat <- exp_dat %>% filter(SNP %in% rsid_pval_exp.clump$rsid)
cat("Finished clump\n")
otc_dat <- read_outcome_data(
    argv$otc,
    sep = "\t", 
    snp_col = "variant_id",
    beta_col = "beta",
    se_col = "standard_error",
    chr_col = "chromosome",
    pos_col = "base_pair_location",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency",
    pval_col = "p_value"
    )
if (argv$proxy) {
    cat("Searching proxy snps...\n")
    otc_dat <- get_proxy(
        exp_dat, 
        otc_dat, 
        r2 = as.numeric(argv$rsq), 
        online = argv$online, 
        pop = argv$pop, 
        bin_plink = argv$plink,
        ref_panel = argv$bfile,
        outfile = argv$prefix
        )
}
otc_dat <- otc_dat %>% mutate(outcome = argv$otcname, id.outcome = argv$otcname, originalname.outcome = argv$otcname)
cat("Harmonising data......\n")
exp_otc_dat <- harmonise_data(exp_dat, otc_dat)
exp_otc_dat <- exp_otc_dat %>% mutate(samplesize.exposure = as.numeric(argv$nexp), samplesize.outcome = as.numeric(argv$notc))
cat("Performing MR...\n")
exp_otc_res <- mr(exp_otc_dat, method_list = c("mr_wald_ratio", "mr_ivw", "mr_egger_regression", "mr_weighted_median"))
out.df <- data.frame(
    ivw_or_wald.b = (exp_otc_res %>% filter(method %in% c("Wald ratio", "Inverse variance weighted")))$b, 
    nsnp = (exp_otc_res %>% filter(method %in% c("Wald ratio", "Inverse variance weighted")))$nsnp, 
    ivw_or_wald.se = (exp_otc_res %>% filter(method %in% c("Wald ratio", "Inverse variance weighted")))$se, 
    ivw_or_wald.pval = (exp_otc_res %>% filter(method %in% c("Wald ratio", "Inverse variance weighted")))$pval)
cat("Finished MR\n")
####### sensitive test ##############
cat("Performing sensitive MR...\n")
# 0. Same direction for all method 
same_direction <- all(exp_otc_res$b > 0) | all(exp_otc_res$b < 0)
out.df$same_direction <- same_direction
# 1. Heterogeneity test
if (nrow(exp_otc_dat) >= 2) {
    exp_otc_mrHet <- mr_heterogeneity(exp_otc_dat, method_list = "mr_ivw")
    q.het <- exp_otc_mrHet$Q_pval
    out.df$q.het <- q.het
} else {
    out.df$q.het <- NA
}
# 2. Horizontal pleiotropy test
if (nrow(exp_otc_dat) >= 3) {
    out.df$eggReg.b <- (exp_otc_res %>% filter(method == "MR Egger"))$b
    out.df$eggReg.se <- (exp_otc_res %>% filter(method == "MR Egger"))$se
    out.df$eggReg.pval <- (exp_otc_res %>% filter(method == "MR Egger"))$pval
    out.df$weightMed.b <- (exp_otc_res %>% filter(method == "MR Egger"))$b
    out.df$weightMed.se <- (exp_otc_res %>% filter(method == "MR Egger"))$se
    out.df$weightMed.pval <- (exp_otc_res %>% filter(method == "MR Egger"))$pval
    exp_otc_mrEggerIntercept <- mr_pleiotropy_test(exp_otc_dat)
    out.df$p.eggerIntercept <- exp_otc_mrEggerIntercept$pval
    

} else {
    out.df$p.eggerIntercept <- NA
    out.df$eggReg.b <- NA
    out.df$eggReg.se <- NA
    out.df$eggReg.pval <- NA 
    out.df$weightMed.b <- NA 
    out.df$weightMed.se <- NA
    out.df$weightMed.pval <- NA
}
# 3. global pleiotropy test
if (nrow(exp_otc_dat) >= 4) {
    exp_otc_mrPresso <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                    OUTLIERtest = FALSE, DISTORTIONtest = FALSE, data = exp_otc_dat, SignifThreshold = 0.05)
    out.df$p.globalPleiotropy <- exp_otc_mrPresso$`MR-PRESSO results`$`Global Test`$Pvalue
} else {
    out.df$p.globalPleiotropy <- NA
}
# 4. Reverse causation test
exp_otc_Steiger <- directionality_test(exp_otc_dat)
(exp_otc_Steiger$snp_r2.exposure < exp_otc_Steiger$snp_r2.outcome) & exp_otc_Steiger$steiger_pval < 0.05
out.df$snp_r2.exposure <- exp_otc_Steiger$snp_r2.exposure
out.df$snp_r2.outcome <- exp_otc_Steiger$snp_r2.outcome
out.df$steiger_pval <- exp_otc_Steiger$steiger_pval
out.df$reverse_caution <- (exp_otc_Steiger$snp_r2.exposure < exp_otc_Steiger$snp_r2.outcome) & exp_otc_Steiger$steiger_pval < 0.05
cat("Finished sensitive test\n")
if (out.df$nsnp == 1) {
    out.df$pass_alltests <- ifelse((out.df$ivw_or_wald.pval < 0.05) & (!out.df$reverse_caution), TRUE, FALSE)
} else if (out.df$nsnp == 2) {
    out.df$pass_alltests <- ifelse((out.df$ivw_or_wald.pval < 0.05) & (out.df$q.het > 0.05) & (!out.df$reverse_caution), TRUE, FALSE)
} else if (out.df$nsnp == 3) {
    out.df$pass_alltests <- ifelse((out.df$ivw_or_wald.pval < 0.05) & (out.df$q.het > 0.05) & (out.df$p.eggerIntercept > 0.05) &(out.df$same_direction) & (!out.df$reverse_caution), TRUE, FALSE)
} else {
    out.df$pass_alltests <- ifelse((out.df$ivw_or_wald.pval < 0.05) & (out.df$q.het > 0.05) & (out.df$p.eggerIntercept > 0.05) &(out.df$same_direction) & (!out.df$reverse_caution) & (out.df$p.globalPleiotropy > 0.05), TRUE, FALSE)
}
# write to file
cat("Write to file... \n")
out.tsv.path <- paste(argv$prefix, ".tsv", sep = "")
write_tsv(out.df, out.tsv.path)
cat("Successfully write results to ", out.tsv.path, "\n", sep="")
setwd("~/phdProject/MR")
library(gwasrapidd)
library(tidyverse)
library(readxl)
library(MRInstruments)
library(epigraphdb)
library(TwoSampleMR)
library(MRPRESSO)
load(file = "~/phdProject/MR/frailtyMR.RData")
mr_filter <- function(exp_otc_dat, exp_otc_res) {
    p.bonf <- p.adjust((exp_otc_res %>% filter(method %in% c("Wald ratio", "Inverse variance weighted")))$pval, method = "bonferroni")
    p.bh <- p.adjust((exp_otc_res %>% filter(method %in% c("Wald ratio", "Inverse variance weighted")))$pval, method = "BH")
    exp_otc_wald.id <- paste((exp_otc_res %>% filter(method == "Wald ratio"))$id.exposure, (exp_otc_res %>% filter(method == "Wald ratio"))$id.outcome)
    exp_otc_ple <- mr_pleiotropy_test(exp_otc_dat)
    if (length(exp_otc_ple > 1)) {
        exp_otc_ple.id <- paste((exp_otc_ple %>% filter(pval > 0.05))$id.exposure, (exp_otc_ple %>% filter(pval > 0.05))$id.outcome)
    } else {
        exp_otc_ple.id = NA
    }
    exp_otc_het <- mr_heterogeneity(exp_otc_dat)
    if (length(exp_otc_het > 1)) {
        exp_otc_het.id <- paste((exp_otc_het %>% filter(Q_pval > 0.05))$id.exposure, (exp_otc_het %>% filter(Q_pval > 0.05))$id.outcome)
    } else {
        exp_otc_het.id = NA
    }
    exp_otc_select.id <- union(exp_otc_wald.id, intersect(exp_otc_ple.id, exp_otc_het.id))
    exp_otc_select.id <- exp_otc_select.id[!is.na(exp_otc_select.id)]
    exp_otc_res.bh <- exp_otc_res %>% mutate(idcol = paste(id.exposure, id.outcome)) %>% filter(method %in% c("Wald ratio", "Inverse variance weighted")) %>% mutate(p.bh = p.bh) %>% filter(p.bh < 0.05, idcol %in% exp_otc_select.id)
    exp_otc_res.bonf <- exp_otc_res %>% mutate(idcol = paste(id.exposure, id.outcome)) %>% filter(method %in% c("Wald ratio", "Inverse variance weighted")) %>% mutate(p.bonf = p.bonf) %>% filter(p.bonf < 0.05, idcol %in% exp_otc_select.id)
    all_res <- list(exp_otc_res.bh = exp_otc_res.bh, exp_otc_res.bonf = exp_otc_res.bonf, exp_otc_ple = exp_otc_ple, exp_otc_het = exp_otc_het)
    return(all_res)
}
head(risk_fi_res %>% mutate(idcol = paste(id.exposure, id.outcome)))
############ FI risk factors ############
fi_risk_factors <- read_lines("~/phdProject/MR/fi_risk_factors.txt")

fi_risk_exposures <- query_epigraphdb(
  route="/mr",
  params=list(exposure_trait=NULL, outcome_trait=fi_risk_factors[1], pval_threshold=1e-08),
  mode="table"
)
for (id in res[2:length(res)]) {
  temp.res <-  query_epigraphdb(
  route="/mr",
  params=list(exposure_trait=NULL, outcome_trait=fi_risk_factors[1], pval_threshold=1e-08),
  mode="table"
  )
  fi_risk_exposures <- rbind(fi_risk_exposures, temp.res)
}
fi_risk.id <- unique(fi_risk_exposures$exposure.id)

##########   MR.0 risk factors - FI  ##########
fi.dfname <- filter(ao, grepl("[Ff]railty.*", trait))$id
risk_fi_id <- filter(ao, id %in% fi_risk.id)$id
risk_exp_dat <- extract_instruments(risk_fi_id)
risk_exp_dat <- clump_data(risk_exp_dat, clump_r2 = 0.1)
fi_otu_dat0 <- extract_outcome_data(snps = risk_exp_dat$SNP, outcomes =  fi.dfname)
risk_fi_dat <- harmonise_data(risk_exp_dat, fi_otu_dat0)
risk_fi_res <- TwoSampleMR::mr(risk_fi_dat)
risk_fi.ntests <- nrow(risk_fi_res %>% distinct(id.exposure))
risk_fi.p_threshold <- 0.05 / risk_fi.ntests
risk_fi_res %>% filter(id.exposure %in% (risk_fi_res %>% filter(method %in% c("Wald ratio", "Inverse variance weighted"), pval <= risk_fi.p_threshold))$id.exposure)
risk_fi_res %>% filter(id.exposure %in% (risk_fi_res %>% filter(method %in% c("Wald ratio", "Inverse variance weighted"), pval <= 0.05))$id.exposure) %>% distinct(exposure)

##########   MR.1 FI - metabolites  ###########
data(gwas_catalog)
head(gwas_catalog)
# metabo_studies <- get_studies(pubmed_id = "35050183")
# save(metabo_studies, file = "~/phdProject/MR/metabo_studies.RData")
# load(file = "~/phdProject/MR/available_outcomes.RData")
# ao <- available_outcomes()
# save(ao, file = "~/phdProject/MR/available_outcomes.RData")
# metabo_associations <- get_associations(study_id = metabo_studies@studies$study_id)
# save(metabo_associations, file = "~/phdProject/MR/metabo_associations.RData")
# load(file = "~/phdProject/MR/metabo_studies.RData")
# load(file = "~/phdProject/MR/metabo_associations.RData")
# load(file = "~/phdProject/MR/available_outcomes.RData")
# metabo_dat <- metabo_associations@associations %>%
#     mutate(
#         SNP = metabo_associations@risk_alleles$variant_id,
#         beta = beta_number * case_when(
#             beta_direction == "decrease" ~ -1,
#             beta_direction == "increase" ~ 1
#         ),
#         effect_allele = metabo_associations@risk_alleles$risk_allele,
#         association_id = as.character(association_id)
#     ) %>%
#     select(SNP, beta, standard_error, effect_allele, pvalue, association_id)
# mtb_exp_dat <- format_data(
#     metabo_dat,
#     snp_col = "SNP",
#     beta_col = "beta",
#     se_col = "standard_error",
#     other_allele_col = "", 
#     effect_allele_col = "effect_allele",
#     pval_col = "pvalue",
#     phenotype_col = "association_id",
#     type = "exposure"
# )
# mtb_gwas <- subset(gwas_catalog, PubmedID == 24816252, Phenotype_info = "Blood metabolite levels")
# mtb_gwas <- mtb_gwas %>% rename(Phenotype = "Phenotype_info", Phenotype_info = "Phenotype") %>% 
            # mutate(Phenotype = str_replace_all(Phenotype , "(^\\()|(\\)$)", ""))
# mtb_exp_dat <- format_data(mtb_gwas, type = "exposure") %>% mutate(exposure = str_replace(exposure, " \\(.*", ""))
# mtb_otc_dat <- format_data(mtb_gwas, type = "outcome") %>% mutate(outcome = str_replace(outcome, " \\(.*", ""))
# mtb_exp_dat <- clump_data(mtb_exp_dat)
mtb_id.list <- filter(ao, pmid == "24816252")$id
mtb_exp_dat <- extract_instruments(mtb_id.list)
mtb_exp_dat <- clump_data(mtb_exp_dat, clump_r2 = 0.1)
fi_otc_dat <- extract_outcome_data(snps=mtb_exp_dat$SNP, outcomes =  fi.dfname)
fi_exp_dat <- extract_instruments(fi.dfname)
mtb_otc_dat <- extract_outcome_data(snps=fi_exp_dat$SNP, outcomes = mtb_id.list)


######## mtb to fi ######### 
mtb_exp_dat <- clump_data(mtb_exp_dat)
mtb_fi_dat <- harmonise_data(mtb_exp_dat, fi_otc_dat)
mtb_fi_res <- TwoSampleMR::mr(mtb_fi_dat)
####### fi to mtb ##########
fi_exp_dat <- clump_data(fi_exp_dat)
fi_mtb_dat <- harmonise_data(fi_exp_dat, mtb_otc_dat)
fi_mtb_res <- TwoSampleMR::mr(fi_mtb_dat)



############## MR.2 FI - disease ##############
########### disease to FI ###################
disease_studies <- c(
    "ieu-a-44",
    "ieu-a-1058",
    "ieu-a-12",
    "ieu-a-996",
    "ieu-a-1054",
    "finn-b-K11_IBD",
    "ieu-a-1025",
    "ieu-a-1112",
    "ieu-a-833",
    "ieu-a-815",
    "finn-b-E4_DM1",
    "ieu-a-973",
    "ieu-a-1169",
    "ieu-a-1170",
    "ieu-a-1171",
    "ieu-a-975",
    "ieu-a-1126",
    "ieu-a-1013",
    "ieu-a-965",
    "ieu-a-966",
    "ieu-a-816",
    "ieu-a-822",
    "ieu-a-967",
    "ieu-a-1082",
    "ieu-b-4965",
    "ieu-a-1109",
    "ieu-a-7",
    "ieu-a-814",
    "ieu-a-1108",
    "ieu-a-1110",
    "ieu-a-1111",
    "ieu-a-23",
    "ieu-a-1102",
    "ieu-a-1081",
    "ieu-a-1183",
    "ieu-b-5067",
    "ieu-a-1086",
    "ieu-a-1186",
    "ieu-a-1185",
    "ieu-a-990",
    "ieu-a-1029",
    "ieu-a-1188",
    "ieu-a-1189",
    "ieu-a-812",
    "ieu-a-803",
    "ieu-a-22",
    "finn-b-I9_HYPTENS")
dis_exp_dat <- extract_instruments(disease_studies)
dis_exp_dat <- clump_data(dis_exp_dat)
fi_otc_dat2 <- extract_outcome_data(snps=dis_exp_dat$SNP, outcomes =  fcisi.dfname)
dis_fi_dat <- harmonise_data(dis_exp_dat, fi_otc_dat2)
dis_fi_res <- TwoSampleMR::mr(dis_fi_dat)
dis_fi_res %>% filter(pval <= 0.05)

########### FI to disease ##################
dis_otc_dat <- extract_outcome_data(snps=fi_exp_dat$SNP, outcomes = disease_studies)
fi_dis_dat <- harmonise_data(fi_exp_dat, dis_otc_dat)
fi_dis_res <- TwoSampleMR::mr(fi_dis_dat)
fi_dis_res.ntests <- nrow(fi_dis_res %>% distinct(id.outcome))
fi_dis.p_threshold <- 0.05 / fi_dis_res.ntests

fi_dis_res %>% filter(id.outcome %in% (fi_dis_res %>% filter(method %in% c("Wald ratio", "Inverse variance weighted"), pval <= fi_dis.p_threshold))$id.outcome) %>% filter(method %in% c("Wald ratio", "Inverse variance weighted"))
fi_dis_res %>% filter(id.outcome %in% (fi_dis_res %>% filter(method %in% c("Wald ratio", "Inverse variance weighted"), pval <= 0.05))$id.outcome) %>% filter(method %in% c("Wald ratio", "Inverse variance weighted"))
write_tsv(
fi_dis_res %>% filter(id.outcome %in% (fi_dis_res %>% filter(method %in% c("Wald ratio", "Inverse variance weighted"), pval <= 0.05))$id.outcome) %>% filter(method %in% c("Wald ratio", "Inverse variance weighted")),
file = "./fi_dis_res.tsv"
)
fi_dis_pleio <- mr_pleiotropy_test(fi_dis_dat)
fi_dis_heterogeneity <- mr_heterogeneity(fi_dis_dat)
fi_dis_leave1 <- mr_leaveoneout(fi_dis_dat)
########## MR.3 FI - metagenome #################
mt_studies <- get_studies(pubmed_id = "35115689")
mt_associations <-get_associations(study_id = mt_studies@studies$study_id)
mt_variants <- get_variants(study_id = mt_studies@studies$study_id)
mt_traits <- get_traits(study_id = mt_studies@studies$study_id)


########## MR.4 FI - proteome ###################
prt_id.list <- filter(ao, pmid == "29875488")$id
prt_exp_dat <- extract_instruments(prt_id.list)
prt_exp_dat <- clump_data(prt_exp_dat, clump_r2 = 0.10)
prt_exp_dat %>% mutate(exposure = str_replace_all(exposure, " \\|\\|.*", "")) %>% filter(exposure %in% pqtl_info$`Target fullname`, SNP %in% pqtl_info$`Sentinel variant*`)
# prt_exp_dat.cisqtl <- prt_exp_dat %>% filter(SNP %>% 
# save(prt_exp_dat, file = "~/phdProject/MR/prt_exp_dat.RData")
fi_otc_dat3 <- extract_outcome_data(snps=prt_exp_dat$SNP, outcomes =  fi.dfname)
prt_fi_dat <- harmonise_data(prt_exp_dat, fi_otc_dat3)
prt_fi_res <- TwoSampleMR::mr(prt_fi_dat)
prt_fi_res.ntests <- nrow(prt_fi_res %>% distinct(id.exposure))
prt_fi.p_threshold <- 0.05 / prt_fi_res.ntests

prt_fi_res %>% filter(id.exposure %in% (prt_fi_res %>% filter(method %in% c("Wald ratio", "Inverse variance weighted"), pval <= prt_fi.p_threshold))$id.exposure)
prt_fi_res %>% filter(method %in% c("Wald ratio", "Inverse variance weighted") , pval <= prt_fi.p_threshold)
############ cis  to fi #####################
pqtl_info <- read_tsv("~/phdProject/MR/pqtl_summary.txt")
pcis_info <- pqtl_info %>% filter(`cis/ trans` == "cis")
pcis_exp_dat <- format_data(
    pcis_info,
    snp_col = "Sentinel variant*",
    beta_col = "beta",
    se_col = "SE",
    effect_allele_col = "Effect Allele (EA)",
    other_allele_col = "Other Allele (OA)",
    eaf_col = "EAF",
    pval_col = "p",
    type = "exposure",
    phenotype_col = "Target fullname",
    id_col = "Target",
    gene_col = "Mapped gene"
)
pcis_exp_dat <- clump_data(pcis_exp_dat, clump_r2 = 0.1)
fi_otc_dat6 <- extract_outcome_data(snps=pcis_exp_dat$SNP, outcomes =  fi.dfname)
pcis_fi_dat <- harmonise_data(pcis_exp_dat, fi_otc_dat6)
pcis_fi_res <- mr(pcis_fi_dat)
# pcis_fi_eggertest <- mr_pleiotropy_test(pcis_fi_dat)
# pcis_fi_eggertest %>% filter(pval >= 0.05)
pci_fi_hetero <- mr_heterogeneity(pcis_fi_dat)
pcis_fi_res.ntests <- nrow(pcis_fi_res %>% distinct(id.exposure))
pcis_fi.p_threshold <- 0.05 / pcis_fi_res.ntests
pcis_fi.p.bh <- p.adjust((pcis_fi_res %>% filter(method %in% c("Wald ratio", "Inverse variance weighted")))$pval, method = "BH")
pcis_fi_signif <- pcis_fi_res %>% filter(pval <= pcis_fi.p_threshold, method %in% c("Wald ratio", "Inverse variance weighted"))
pcis_fi_signif.bh <- (pcis_fi_res %>% filter(method %in% c("Wald ratio", "Inverse variance weighted")))[pcis_fi.p.bh<=0.05,] 
pcis_fi_res_filter <- mr_filter(pcis_fi_dat, pcis_fi_res)
############### fi to cis ###############
pcis_id.list <- (filter(ao, trait %in% pcis_info$`Target fullname`))$id
pcis_otc_dat <- extract_outcome_data(snps = fi_exp_dat$SNP, outcomes = pcis_id.list)
############# cis-pqtl to risk ##########################
pcis_exp_dat2 <- pcis_exp_dat %>% filter(id.exposure %in% pcis_fi_signif$id.exposure)
riskpap_otc_dat <- extract_outcome_data(snps=pcis_exp_dat2$SNP, outcomes = c(fi_riskpap$id, disease_studies))
cis_riskpap_dat <- harmonise_data(pcis_exp_dat2, riskpap_otc_dat)
cis_risk_res <- mr(cis_riskpap_dat)
cis_risk_res %>% filter(pval <= 0.05)
############ risk to cis-pqtl ###################

############ cis + trans #####################
pqtl_info <- read_tsv("~/phdProject/MR/pqtl_summary.txt")
pqtl_exp_dat <- format_data(
    pqtl_info,
    snp_col = "Sentinel variant*",
    beta_col = "beta",
    se_col = "SE",
    effect_allele_col = "Effect Allele (EA)",
    other_allele_col = "Other Allele (OA)",
    eaf_col = "EAF",
    pval_col = "p",
    type = "exposure",
    phenotype_col = "Target fullname",
    id_col = "Target",
    gene_col = "Mapped gene"
)
pqtl_exp_dat <- clump_data(pqtl_exp_dat, clump_r2 = 0.1)
fi_otc_dat7 <- extract_outcome_data(snps=pqtl_exp_dat$SNP, outcomes =  fi.dfname)
pqtl_fi_dat <- harmonise_data(pqtl_exp_dat, fi_otc_dat7)
pqtl_fi_res <- mr(pqtl_fi_dat)

pqtl_fi_res.ntests <- nrow(pqtl_fi_res %>% distinct(id.exposure))
pqtl_fi.p_threshold <- 0.05 / pqtl_fi_res.ntests
pqtl_fi_res %>% filter(pval <= pqtl_fi.p_threshold, method %in% c("Wald ratio", "Inverse variance weighted"))
########### metagenome - FI #################
# mgwas.summary <- read_tsv("/home/ywwang/phdProject/MR/mgwas-pubmedId_32572223.tsv")
library(gwasrapidd)
mgwas_associations <- get_associations(pubmed_id = "32572223")
mgwas_variants <- get_variants(pubmed_id = "32572223")
mtg_dat <- mgwas_associations@associations %>%
    mutate(
        SNP = mgwas_associations@risk_alleles$variant_id,
        beta = beta_number * case_when(
            beta_direction == "decrease" ~ -1,
            beta_direction == "increase" ~ 1
        ),
        effect_allele = mgwas_associations@risk_alleles$risk_allele,
        eaf = mgwas_associations@risk_alleles$risk_frequency,
        pvalue_description = str_replace_all(pvalue_description, "[()]", "")
    ) %>%
    select(SNP, beta, standard_error, effect_allele, eaf, pvalue, pvalue_description)
mtg_exp_dat <- format_data(
    mtg_dat,
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    eaf_col = "eaf",
    pval_col = "pvalue",
    type = "exposure",
    phenotype_col = "pvalue_description"
)
mtg_exp_dat <- mtg_exp_dat %>% mutate(id.exposure = exposure)
# mtg_exp_dat <- clump_data(mtg_exp_dat, clump_r2 = 0.1)
fi_otc_dat4 <- extract_outcome_data(snps=mtg_exp_dat$SNP, outcomes =  fi.dfname)
mtg_fi_dat <- harmonise_data(mtg_exp_dat, fi_otc_dat4)
mtg_fi_res <- mr(mtg_fi_dat)
mtg_fi_res.ntests <- nrow(mtg_fi_res %>% distinct(id.exposure))
mtg_fi.p_threshold <- 0.05 / mtg_fi_res.ntests
mtg_fi_res %>% filter(id.exposure %in% (mtg_fi_res %>% filter(method %in% c("Wald ratio", "Inverse variance weighted"), pval <= mtg_fi.p_threshold))$id.exposure)
mtg_fi_res %>% filter(id.exposure %in% (mtg_fi_res %>% filter(method %in% c("Wald ratio", "Inverse variance weighted"), pval <= 0.05))$id.exposure)
write_tsv(
mtg_fi_res %>% filter(id.exposure %in% (mtg_fi_res %>% filter(method %in% c("Wald ratio", "Inverse variance weighted"), pval <= 0.05))$id.exposure) %>% filter(method %in% c("Wald ratio", "Inverse variance weighted")),
file = "mtg_fi_res.tsv"
)
########### metagenome - FI risk #################
ris_otc_dat <- extract_outcome_data(snps=mtg_exp_dat$SNP, outcomes = fi_risk.id)
mtg_ris_dat <- harmonise_data(mtg_exp_dat, ris_otc_dat)
mtg_ris_res <- mr(mtg_ris_dat)

########### Parabacteroides - disease  #######################
know_disease_id <- c(
    "finn-b-L12_ALOPECAREATA",
    "finn-b-E4_POCS",
    "ebi-a-GCST008036",
    "finn-b-M13_PSORIARTH",
    "finn-b-M13_ANKYLOSPON",
    "finn-b-K11_IBD",
    "ieu-b-18",
    "finn-b-NAFLD",
    "finn-b-O15_PREG_DM",
    "ebi-a-GCST005195",
    "finn-b-E4_OBESITY"
)
# alop.id <- filter(ao, trait == "Alopecia areata")$id
parabac.exp.dat <- mtg_exp_dat %>% filter(grepl("Parabacteroides", id.exposure))
knowd.otc.dat <- extract_outcome_data(snps=parabac.exp.dat$SNP, outcomes = know_disease_id)
parabac_knowd_dat <- harmonise_data(parabac.exp.dat, knowd.otc.dat)
parabac_knowd_mr <- mr(parabac_knowd_dat)

###### metabolome - fi_items ###########
fi_items_id <- (ao %>% filter(trait %in% fi_risk_factors, population == "European"))$id
fi_items_otc_dat <- extract_outcome_data(snps=mtb_exp_dat$SNP, outcomes = fi_items_id)
mtg_fi_items_dat <- harmonise_data(mtg_exp_dat, fi_items_otc_dat)
mtg_fi_items.res <- mr(mtg_fi_items_dat)
mtg_fi_items.res %>% filter(pval <= 0.05)
mtg_fi_items.res.ntests <- nrow(mtg_fi_items.res %>% distinct(id.outcome))


############## fi_risk_paper to frailty ###############
jobs_studies <- ao %>% filter(grepl("Job SOC coding", trait))
riskpap <- read_csv("/home/ywwang/phdProject/MR/fi_risk_frompaper.txt", col_names = c("risk", "id"))
riskpap_exp_dat <- extract_instruments(c(riskpap$id, jobs_studies$id))
riskpap_exp_dat  <- clump_data(riskpap_exp_dat)
fi_otc_dat5 <- extract_outcome_data(snps=riskpap_exp_dat$SNP, outcomes =  fi.dfname)
riskpap_fi_dat <- harmonise_data(riskpap_exp_dat, fi_otc_dat5)
riskpap_fi_res <- TwoSampleMR::mr(riskpap_fi_dat, method_list = c(subset(mr_method_list(), use_by_default)$obj,  "mr_ivw_mre", "mr_ivw_fe"))

############## frailty to fi_risk_paper ###############
############# proteom to fi_risk_casual #################
risk_casual.id <- (riskpap_fi_res %>% filter(pval <= risk_fi.p_threshold, method %in% c("Wald ratio", "Inverse variance weighted")))$id.exposure
risk_casual_otc_dat <- extract_outcome_data(snps=prt_exp_dat$SNP, outcomes =  risk_casual.id)
risk_casual_prt <- harmonise_data(prt_exp_dat, risk_casual_otc_dat)
prt_risk_casual_res <- mr(risk_casual_prt)
prt_risk_casual_res.ntests <- nrow(prt_risk_casual_res   %>% distinct(id.exposure))
prt_risk_casual.p_threshold <- 0.05 / prt_risk_casual_res.ntests
prt_risk_casual_res %>% filter(pval <= 5e-5)
save.image(file = "~/phdProject/MR/frailtyMR.RData")




########## fi to phenome ############
disphenome <- read_lines("ukb_disphenome.txt")
t1 <- ao %>% filter(trait %in% disphenome, population == "European", sex == "Males and Females") %>% arrange(-nsnp) %>% distinct(trait, .keep_all = T)
t2 <- ao %>% filter(grepl("Disease.*", category), population == "European", sex == "Males and Females") %>% arrange(-nsnp) %>% distinct(trait, .keep_all = T)
t3 <- ao %>% filter(grepl("disease.*", trait), population == "European", sex == "Males and Females") %>% arrange(-nsnp) %>% distinct(trait, .keep_all = T)
t4 <- ao %>% filter(grepl("Diagnoses.*", trait), population == "European", sex == "Males and Females") %>% arrange(-nsnp) %>% distinct(trait, .keep_all = T)
t5 <- ao %>% filter(grepl("cancer.*", trait), population == "European", sex == "Males and Females") %>% arrange(-nsnp) %>% distinct(trait, .keep_all = T)
t6 <- ao %>% filter(id %in% c(disease_studies))
t7 <- rbind(t1, t2, t3, t4, t5, t6) %>% arrange(-nsnp) %>% distinct(trait, .keep_all = T) %>% arrange(trait)
disphenome_otc_dat <- extract_outcome_data(snps=fi_exp_dat$SNP, outcomes = t7$id)
fi_disphenome_dat <- harmonise_data(fi_exp_dat, disphenome_otc_dat)
fi_disphenome_res <- mr(fi_disphenome_dat)
fi_disphenome_res_filter <- mr_filter(fi_disphenome_dat, fi_disphenome_res)
################## phenome to fi ###############################
disphenome_exp_dat <- extract_instruments(t7$id)
disphenome_exp_dat <- clump_data(disphenome_exp_dat)
disphenome_fi_otc <- extract_outcome_data(snps = disphenome_exp_dat$SNP, outcomes = fi.dfname)
disphenome_fi_dat <- harmonise_data(disphenome_exp_dat, disphenome_fi_otc)
disphenome_fi_res <- mr(disphenome_fi_dat)
disphenome_fi_res_filter <- mr_filter(disphenome_fi_dat, disphenome_fi_res)
############## selected cis-qtl to disphenome
############# DQA2 and sRAGE to phenome ######################
twocis_exp_dat <-pcis_exp_dat %>% filter(id.exposure %in% c("DQA2", "sRAGE"))
twocis_dis_otc <- extract_outcome_data(snps = twocis_exp_dat$SNP, outcomes = t7$id)
twocis_dis_dat <- harmonise_data(twocis_exp_dat, twocis_dis_otc)
twocis_dis_res <- mr(twocis_dis_dat)
test <- mr_heterogeneity(twocis_dis_dat)
twocis_dis_res_filter <- mr_filter(twocis_dis_dat, twocis_dis_res)
############# phenome to DQA2 and sRAGE ######################
twocis.id <- filter(ao, pmid == "29875488", trait %in% c("HLA class II histocompatibility antigen, DQ alpha 2 chain", "Advanced glycosylation end product-specific receptor, soluble"))$id
# pqtl_otc_dat <- extract_outcome_data(snps = disphenome_exp_dat$SNP, outcomes = prt_id.list)
twocis_otc_dat <- extract_outcome_data(snps = disphenome_exp_dat$SNP, outcomes = twocis.id)
dis_twocis_dat <- harmonise_data(disphenome_exp_dat, twocis_otc_dat)
dis_twocis_res <- mr(dis_twocis_dat)
dis_twocis_res_filter <- mr_filter(dis_twocis_dat, dis_twocis_res)
################ fi - COVID 19 #######################
# covid_severe.gwas <- read_tsv("/home/ywwang/phdProject/MR/COVID19_HGI_A2_ALL_eur_leave23andme_20220403_GRCh37.tsv.gz")
# covid_severe_filter.gwas <- covid_severe.gwas %>% filter(all_inv_var_meta_p < 5e-8) 
# save(covid_severe.gwas, file = "/home/ywwang/phdProject/MR/covid_severe_gwas.RData")
# load("/home/ywwang/phdProject/MR/covid_severe_gwas.RData")
# covid_hosp.gwas <- read_tsv("/home/ywwang/phdProject/MR/COVID19_HGI_B2_ALL_eur_leave23andme_20220403_GRCh37.tsv.gz")
# covid_hosp_filter.gwas <- covid_hosp.gwas %>% filter(all_inv_var_meta_p < 5e-8)
# save(covid_hosp.gwas, file = "/home/ywwang/phdProject/MR/covid_hosp.gwas.RData") 
# load("/home/ywwang/phdProject/MR/covid_hosp.gwas.RData")

# covid_severe.exp <- format_data(
#     covid_severe_filter.gwas %>% mutate(phenotype = "covid_severe"),
#     snp_col = "rsid",
#     beta_col = "all_inv_var_meta_beta",
#     se_col = "all_inv_var_meta_sebeta",
#     effect_allele_col = "ALT",
#     other_allele_col = "REF",
#     eaf_col = "all_meta_AF",
#     pval_col = "all_inv_var_meta_p",
#     type = "exposure",
#     phenotype_col = "phenotype"
# )
# covid_severe_id <- "ebi-a-GCST011075"
# covid_hosp_id <- "ebi-a-GCST011081"
covid.id <- c(covid_severe_id, covid_hosp_id)
covid.id.list <- (ao %>% filter(grepl("COVID-19", trait), population == "European"))$id
covid_exp_dat <- extract_instruments(covid.id.list)
covid_exp_dat <- clump_data(covid_exp_dat)
covid_fi_otc <- extract_outcome_data(snps = covid_exp_dat$SNP, outcomes = fi.dfname)
covid_fi_dat <- harmonise_data(covid_exp_dat, covid_fi_otc)
covid_fi_res <- mr(covid_fi_dat)
covid_fi_res %>% filter(pval <= 0.05)
fi_covid_otc <- extract_outcome_data(snps = fi_exp_dat$SNP, outcomes = covid.id)
fi_covid_dat <- harmonise_data(fi_exp_dat, fi_covid_otc)
fi_covid_res <- mr(fi_covid_dat)
# covid_severe.exp <- clump_data(covid_severe.exp)
# covid_severe_fi_otc <- extract_outcome_data(snps = covid_severe.exp$SNP, outcomes =  fi.dfname)
# covid_severe_fi_dat <- harmonise_data(covid_severe.exp, covid_severe_fi_otc)
# covid_severe_fi_res <- TwoSampleMR::mr(covid_severe_fi_dat)

# covid_severe.otc <- format_data(
#     covid_severe_filter.gwas %>% mutate(phenotype = "covid_severe"),
#     snp_col = "rsid",
#     beta_col = "all_inv_var_meta_beta",
#     se_col = "all_inv_var_meta_sebeta",
#     effect_allele_col = "ALT",
#     other_allele_col = "REF",
#     eaf_col = "all_meta_AF",
#     pval_col = "all_inv_var_meta_p",
#     type = "outcome",
#     phenotype_col = "phenotype"
# )
fi_covid_otc <- extract_outcome_data(snps = fi_exp_dat$SNP, outcomes = covid.id.list)
fi_covid_dat <- harmonise_data(fi_exp_dat, fi_covid_otc)
fi_covid_res <- mr(fi_covid_dat)
fi_covid_res %>% filter(pval <= 0.05)
save(fi_covid_res, file = "./fi_covid_res.RData")
write_lines(fi_exp_dat$SNP, file = "fi_exp_snp.txt")

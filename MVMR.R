library(dplyr)
library(TwoSampleMR)
library(MendelianRandomization)
library(MVMR)
library(xlsx)
library(data.table)

###TwosampleMR包
##import data
setwd("E:/reproductive/rawdata")
a<-fread("GCST90029036_buildGRCh37.tsv",header=T)
b<-fread("repro_MENOPAUSE_AGE.sumstats.gz",header=T)
c<-fread("GCST90000045_buildGRCh37.tsv.gz",header=T)
d<-fread("GCST90000048_buildGRCh37.tsv.gz",header=T)
names(a)<-c("SNP","CHR","BP","effect_allele","other_allele","REF","eaf","beta","se","pval","samplesize","INFO")
names(b)[4:11]<-c("other_allele","effect_allele","REF","eaf","beta","se","pval","samplesize")
names(c)[1:6]<-c("SNP","effect_allele","Other_allele","beta","se","pval")
names(d)[1:6]<-c("SNP","effect_allele","Other_allele","beta","se","pval")
exposure_dat<-mv_extract_exposures_local(c("menarche_raw.csv","menopause_raw.csv","AFS_raw.csv","afb_raw.csv"),
                                         sep = ",",snp_col = "SNP",
                                         beta_col = "beta",
                                         se_col = "se",
                                         eaf_col = "eaf",
                                         effect_allele_col = "effect_allele",
                                         other_allele_col = "other_allele",
                                         pval_col = "pval",samplesize_col = "samplesize",min_pval = 1e-200,
                                         log_pval = FALSE,
                                         pval_threshold = 5e-08,
                                         clump_r2 = 0.001,
                                         clump_kb = 10000,
                                         harmonise_strictness = 2)
write.csv(exposure_dat,"mv_exposure_dat.csv")
setwd("E:/Rstudy/MRGWASdata/ADdata")
mvmr_ad<-fread("AD_sumstats_Jansenetal_2019sept.txt",header = T)
names(mvmr_ad)[4:14]<-c("effect allele","other allele","SNP","Z","p","Nsum","Neff","dir","eaf","beta","se")
exp_out_merged<-merge(exposure_dat,mvmr_ad,by.x = "SNP",by.y = "SNP")
write.csv(exp_out_merged,file = "outcome_mvmr_merge.csv")
outcome_dat<-read_outcome_data(snps =exposure_dat$SNP,filename = "outcome_mvmr_merge.csv",sep = ",",snp_col = "SNP",beta_col = "beta",se_col = "se",effect_allele_col = "effect allele",other_allele_col = "other allele",pval_col = "p")
mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)

##data standardization
MRInputObject <- mr_mvinput(bx = mvdat$exposure_beta,
                            bxse= mvdat$exposure_se,
                            by = mvdat$outcome_beta,
                            byse = mvdat$outcome_se,
                            snps = rownames(mvdat$exposure_beta),
                            exposure = c("menar","meno","AFS","AFB"),outcome = "AD")

mv_end_dat1<- format_mvmr(mvdat$exposure_beta, mvdat$outcome_beta, mvdat$exposure_se, mvdat$outcome_se,rownames(mvdat$exposure_beta))
mvmr_res <- mvmr(mv_end_dat1, 0, 1)
save(mvmr_res, file="mvmr_res_dpwandcsi_final.Rdata")

##MR-IVW, MR-Egger
mr_mvivw <- mr_mvivw(MRInputObject)
save(mr_mvivw, file="mr_mvivw_dpwandcsi_final.Rdata")
mr_mvegger <- mr_mvegger(MRInputObject, orientate = 1)
save(mr_mvegger, file="mr_mvegger_dpwandcsi_final.Rdata")
mvmr_results_menar <- c(exp(mr_mvivw$Estimate[1]), exp(mr_mvivw$CILower[1]), exp(mr_mvivw$CIUpper[1]), mr_mvivw$Pvalue[1], exp(mr_mvegger$Estimate[1]), exp(mr_mvegger$CILower.Est[1]), exp(mr_mvegger$CIUpper.Est[1]), mr_mvegger$Pvalue.Est[1])
mvmr_results_menar <- c(format(mvmr_results_menar, scientific=F))
names(mvmr_results_menar) <- c("IVW_OR", "IVW_CIL", "IVW_CIU", "IVW_P", "Egger_OR", "Egger_CIL", "Egger_CIU", "Egger_P")

mvmr_results_meno <- c(exp(mr_mvivw$Estimate[2]), exp(mr_mvivw$CILower[2]), exp(mr_mvivw$CIUpper[2]), mr_mvivw$Pvalue[2], exp(mr_mvegger$Estimate[2]), exp(mr_mvegger$CILower.Est[2]), exp(mr_mvegger$CIUpper.Est[2]), mr_mvegger$Pvalue.Est[2])
mvmr_results_meno <- c(format(mvmr_results_meno, scientific=F))
names(mvmr_results_meno) <- c("IVW_OR", "IVW_CIL", "IVW_CIU", "IVW_P", "Egger_OR", "Egger_CIL", "Egger_CIU", "Egger_P")

mvmr_results_AFS <- c(exp(mr_mvivw$Estimate[3]), exp(mr_mvivw$CILower[3]), exp(mr_mvivw$CIUpper[3]), mr_mvivw$Pvalue[3], exp(mr_mvegger$Estimate[3]), exp(mr_mvegger$CILower.Est[3]), exp(mr_mvegger$CIUpper.Est[3]), mr_mvegger$Pvalue.Est[3])
mvmr_results_AFS <- c(format(mvmr_results_AFS, scientific=F))
names(mvmr_results_AFS) <- c("IVW_OR", "IVW_CIL", "IVW_CIU", "IVW_P", "Egger_OR", "Egger_CIL", "Egger_CIU", "Egger_P")

mvmr_results_AFB <- c(exp(mr_mvivw$Estimate[4]), exp(mr_mvivw$CILower[4]), exp(mr_mvivw$CIUpper[4]), mr_mvivw$Pvalue[4], exp(mr_mvegger$Estimate[4]), exp(mr_mvegger$CILower.Est[4]), exp(mr_mvegger$CIUpper.Est[4]), mr_mvegger$Pvalue.Est[4])
mvmr_results_AFB <- c(format(mvmr_results_AFB, scientific=F))
names(mvmr_results_AFB) <- c("IVW_OR", "IVW_CIL", "IVW_CIU", "IVW_P", "Egger_OR", "Egger_CIL", "Egger_CIU", "Egger_P")

mvmr_results <- cbind(mvmr_results_menar, mvmr_results_meno,mvmr_results_AFS,mvmr_results_AFB)  
write.csv(mvmr_results, "mvmr_results_ivwandegger.csv") 
mvmregger_results <- c(mr_mvegger$Intercept[1], mr_mvegger$StdError.Int[1], mr_mvegger$CILower.Int[1], mr_mvegger$CIUpper.Int[1], mr_mvegger$Pvalue.Int[1])
names(mvmregger_results) <- c("intercept", "stderror.int", "cilower.int", "ciupper.int", "pvalue.int")
write.csv(mvmr_results, "mvmr_results_egger.csv")

##Instrument strength 
strength_mvmr <- strength_mvmr(mv_end_dat1, gencov=0)
pleiotropy_mvmr <- pleiotropy_mvmr(mv_end_dat1, gencov=0)

##MR-LASSO
lasso1<-mr_mvlasso(MRInputObject,
                   +                    orientate = 1,
                   +                    distribution = "normal", alpha = 0.05,
                   +                    lambda = numeric(0))
mvmr_lasso_estimate<-cbind(lasso1@Exposure,lasso1@Estimate,lasso1@StdError,lasso1@CILower,lasso1@CIUpper,lasso1@Pvalue)
names(mvmr_lasso_estimate)<-c("Exposure","Estimate","StdError","CILower","CIUpper","Pvalue")
write.csv(mvmr_lasso_estimate,"mvmr_lasso_estimate.csv")
working.dir <- '~/Desktop/BBL/data/joy/BBL/'
setwd(working.dir)
#genomics <- read.csv(paste0(working.dir,"/studies/pnc/n9498_dataFreeze/genetics/pnc_dem_dx_cnb_prs_psy_20170522.csv"))
###############################
######sample construction######
###############################
load(paste0(working.dir,'projects/prsConnectivity/sample.RData')) #sample constructed (n=1015)
#sample_genomics <- merge(genomics,sample_merge_qa, by='bblid')
il_ea_prs <- read.csv('~/Desktop/illumina_ea_pgriskscore.csv') # PRS for euro ancestry
#volData <- read.csv('~/Desktop/volumeDataForGeneStudy_20171030.csv') #vol data from Chad for scanid matching
#sample_geno_volData <- merge(sample_genomics,volData, by = 'scanid')
#sample_geno_volData_il_ea <- merge(sample_geno_volData,il_ea_prs, by.x = 'bblid.y', by.y = 'pnc_id')
#colnames(sample_geno_volData_il_ea)[which(colnames(sample_geno_volData_il_ea) == "bblid.y")] <- 'pnc_id'
#colnames(sample_geno_volData_il_ea)[which(colnames(sample_geno_volData_il_ea) == "bblid.x")] <- 'bblid'
n1601_ids <- read.csv('~/Desktop/n1601_dbgap_newids.csv')
colnames(n1601_ids) <- c('pnc_id','scanid')
n1601_prs<-merge(il_ea_prs,n1601_ids)
n1601_prs_sample <- merge(n1601_prs,sample_merge_qa)

###############################
##########Get networks#########
###############################
n1601_prs_sample_adj <- get_adj(n1601_prs_sample,'power')
power_com <- read.csv(paste0(working.dir,'studies/pnc/n1601_dataFreeze/neuroimaging/rest/restNetwork_power/Consensus264.csv'))
power_com_nouk <- subset(power_com, Community > 0)

########################################
##########all original networks#########
########################################
prs_adj_com <- lapply(n1601_prs_sample_adj,function(adj) get_community_weight(adj,power_com))
power_labels <- rev(c("CRB","DAT","VAT","SBC","SAL","FPT","VIS","MEM",
                               "DMN","AUD","COP","SMT-M","SMT-H"))
power_com_vec <- get_com_mat_lables(power_labels)

prs_adj_com_plots <- lapply(prs_adj_com,function(com_adj) plot_community_adj(com_adj,power_labels,0.001))
ave_adj_com_plot<-plot_community_adj(apply(simplify2array(prs_adj_com), 1:2, mean),power_labels,0.0001)

com_wt_flat<-t(sapply(prs_adj_com, function(com_adj) com_adj[upper.tri(com_adj,diag = TRUE)]))
com_wt_flat_data <- data.frame(age = n1601_prs_sample$ageAtScan1,sex = as.factor(n1601_prs_sample$sex),
                               motion = n1601_prs_sample$restRelMeanRMSMotion, scz_prs = n1601_prs_sample$scz_prs, com_wt_flat)

com_wt_gam <- lapply(5:95,function(x) gam(com_wt_flat_data[,x] ~ s(age) + sex + motion + scz_prs , data = com_wt_flat_data, method="REML"))
com_wt_gam_ptable <- sapply(com_wt_gam, function(gam) summary(gam)$p.table )
scz_prs_pval <- com_wt_gam_ptable[dim(com_wt_gam_ptable)[1],]

scz_prs_results<-data.frame(community = power_com_vec, pval_raw = scz_prs_pval, pval_fdr = p.adjust(scz_prs_pval,method = 'fdr'))
subset(scz_prs_results, pval_raw < 0.05)
com_sig <- which(scz_prs_results$pval_raw < 0.05)
cov = y ~ s(age) + sex + motion + scz_prs
scz_prs_gam_plots<-lapply(as.list(com_sig),function(sig_com) plot_gam(x= 'scz_prs', y = com_wt_flat[,sig_com], model =  cov,
                                                                      data = com_wt_flat_data, pval_results = scz_prs_results[sig_com,]))
plot_grid(plotlist = scz_prs_gam_plots,nrow = 3)

########################################
################less  networks##########
########################################
power_com_merge <-power_com
power_com_merge$Community <- 0
new_SMT <- which(power_com$System == 'Sensory/somatomotor Hand' | power_com$System == 'Sensory/somatomotor Mouth' | power_com$System == 'Auditory' | power_com$System == 'Visual')
new_DMN <- which(power_com$System == 'Default mode' | power_com$System == 'Memory retrieval?')
new_ATT <- which(power_com$System == 'Ventral attention' | power_com$System == 'Dorsal attention')
new_nonCX <- which(power_com$System == 'Subcortical' | power_com$System == 'Cerebellar')
new_COP <- which(power_com$System == 'Cingulo-opercular Task Control' | power_com$System == 'Salience')
new_AUD_VIS <- which(power_com$System == "Auditory" | power_com$System == "Visual")

power_com_merge[new_SMT,'Community'] <- 1
power_com_merge[c(new_ATT,new_COP),'Community'] <- 2
#power_com_merge[new_AUD_VIS,'Community'] <- 3
power_com_merge[new_DMN,'Community'] <- 3
#power_com_merge[new_ATT,'Community'] <- 4
power_com_merge[which(power_com$System == "Fronto-parietal Task Control"),'Community'] <- 4
power_com_merge[new_nonCX,'Community'] <- 0

prs_adj_com <- lapply(n1601_prs_sample_adj,function(adj) get_community_weight(adj,power_com_merge))
power_labels <- rev(c("FPT","DMN","COP/ATT","SMT"))
power_com_vec <- get_com_mat_lables(power_labels)

prs_adj_com_plots <- lapply(prs_adj_com,function(com_adj) plot_community_adj(com_adj,power_labels,0.0001))
ave_adj_com_plot<-plot_community_adj(apply(simplify2array(prs_adj_com), 1:2, mean),power_labels,0.0001)
com_wt_flat<-t(sapply(prs_adj_com, function(com_adj) com_adj[upper.tri(com_adj,diag = TRUE)]))
com_wt_flat_data <- data.frame(age = n1601_prs_sample$ageAtScan1,sex = as.factor(n1601_prs_sample$sex),
                               motion = n1601_prs_sample$restRelMeanRMSMotion, scz_prs = n1601_prs_sample$scz_prs, com_wt_flat)


com_wt_gam <- lapply(5:dim(com_wt_flat_data)[2],function(x) gam(com_wt_flat_data[,x] ~ s(age) + sex + motion + scz_prs , data = com_wt_flat_data, method="REML"))
com_wt_gam_ptable <- sapply(com_wt_gam, function(gam) summary(gam)$p.table )

scz_prs_pval <- com_wt_gam_ptable[dim(com_wt_gam_ptable)[1],]

com_sig<-which(p.adjust(scz_prs_pval,method = 'fdr') <0.05)
power_com_vec[com_sig]

scz_prs_results<-data.frame(community = power_com_vec, pval_raw = scz_prs_pval, pval_fdr = p.adjust(scz_prs_pval,method = 'fdr'))
cov = y ~ s(age) + sex + motion + scz_prs
scz_prs_gam_plots<-lapply(as.list(com_sig),function(sig_com) plot_gam(x= 'scz_prs', y = com_wt_flat[,sig_com], model =  cov,
                                                   data = com_wt_flat_data, pval_results = scz_prs_results[sig_com,]))
plot_grid(plotlist = scz_prs_gam_plots,nrow = 3)
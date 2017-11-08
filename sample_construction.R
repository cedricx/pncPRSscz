working.dir <- '~/Desktop/BBL/data/joy/BBL/'
setwd(working.dir)
#genomics <- read.csv(paste0(working.dir,"/studies/pnc/n9498_dataFreeze/genetics/pnc_dem_dx_cnb_prs_psy_20170522.csv"))
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

n_sample <- dim(sample_qa)[1]
sample_net<-array(NA, c(264, 264, n_sample))
for (i in 1:n_sample){
  scanid <- sample_qa$scanid[i]
  netpath<- paste("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/rest/restNetwork_264PowerPNC/264PowerPNCNetworks/",scanid,"_264PowerPNC_network.txt",sep="")
  sample_net[,,i] <- as.matrix(read.table(netpath))
  print(paste(i,"."," copying ",scanid,"_","Power",sep=""))
}
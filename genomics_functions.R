get_power_path <- function(i,n,scanid){
  netpath<- paste0(working.dir,"studies/pnc/n1601_dataFreeze/neuroimaging/rest/restNetwork_264PowerPNC/264PowerPNCNetworks/",scanid,"_264PowerPNC_network.txt")
  print(paste0(i,'/',n,'-',"Copying ",scanid,"_","Power"))
  return(netpath)
}

get_adj <- function(sample,modality){
  scanid_list <- as.list(sample$scanid)
  nsample <- length(scanid_list)
  if (modality == 'power') {
    sample_adj<-lapply(seq_along(scanid_list), function(i) as.matrix(read.table(get_power_path(i,nsample,scanid_list[[i]]))))
  }
  return(sample_adj)
}

get_community_weight <- function(adj, com_assign) {
  n_com <- max(com_assign$Community)
  adj_com <- matrix(NA, nrow = n_com, ncol = n_com)
  for (com_i in 1:n_com) {
    nodes_in_com_i <- which(com_assign$Community == com_i)
    for (com_j in 1:n_com) {
      nodes_in_com_j <- which(com_assign$Community == com_j)
      adj_com_temp <- adj[nodes_in_com_i,nodes_in_com_j]
      if (com_i == com_j) {
      adj_com[com_i,com_j] <-mean(adj_com_temp[upper.tri(adj_com_temp)])}
      else {
        adj_com[com_i,com_j] <-mean(adj_com_temp)  
      }
    }
  }
  return(adj_com)
}

plot_community_adj <- function(com_adj, com_labels, zero){
  keycol = c(rev(brewer.pal(7,"Reds")),'white',brewer.pal(7,"Blues"))
  half_range <- seq(1,zero,length.out = 8)
  full_range <- c(half_range,rev(-half_range))
  com_adj_plot <-levelplot(com_adj, col.regions = rev(keycol), at = full_range,
                           xlab="",ylab = "",
                           scales=list(x=list(at = 1:13,labels=com_labels,rot=90, tck = 0),
                                       y=list(at = 1:13,labels=com_labels, tck = 0)))
}

get_com_mat_lables <- function(com_labels) {
  com_name_mat<- matrix(NA,nrow = length(com_labels), ncol = length(com_labels))
  for (com_i in com_labels) {
    for (com_j in com_labels) {
      i = which(com_labels == com_i)
      j = which(com_labels == com_j)
      com_name_mat[i,j] <- paste0(com_i,'-',com_j)
    }
  }
  com_mat_vec<-com_name_mat[upper.tri(com_name_mat,diag = T)]
  return(com_mat_vec)
}


plot_gam <- function(x,y, model, data, pval_results){
  data = cbind(data, y)
  gam_model <- gam(eval(model), data = data, method = "REML")
  plotdata <- visreg(gam_model,x,type = "conditional",scale = "linear", plot = FALSE)
  smooths <- data.frame(Variable = plotdata$meta$x, 
                        x=plotdata$fit[[plotdata$meta$x]], 
                        smooth=plotdata$fit$visregFit, 
                        lower=plotdata$fit$visregLwr, 
                        upper=plotdata$fit$visregUpr)
  predicts <- data.frame(Variable = x, 
                         x=plotdata$res[,x],
                         y=plotdata$res$visregRes)
  
  p_fdr_text <- paste0("p_fdr=", round(pval_results$pval_fdr,3))
  p_raw_text <- paste0("p_raw=", round(pval_results$pval_raw,4))
  
  gam_plot<-ggplot() +
    geom_point(data = predicts, aes(x, y), colour = "grey" ) +
    geom_line(data = smooths, aes(x = x, y = smooth)) +
    geom_line(data = smooths, aes(x = x, y=lower), linetype="dashed") + 
    geom_line(data = smooths, aes(x = x, y=upper), linetype="dashed") +
    annotate("text",x = -Inf, y = Inf, hjust = -0.1,vjust = 1,label = p_fdr_text ,size = 5, colour = "black",fontface ="italic" ) +
    annotate("text",x = -Inf, y = Inf, hjust = -0.1,vjust = 3,label = p_raw_text ,size = 5, colour = "black",fontface ="italic" ) +
    theme(legend.position = "none") +
    ggtitle(paste0(toString(model[2]),toString(model[1]),toString(model[3]))) +
    labs(x = x, y = pval_results$community) 
  return(gam_plot)
}







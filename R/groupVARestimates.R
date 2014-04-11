groupVARestimates <- function(object)
{
  mats <- extractParFit(object$Model$Results)
  
  Vars <- rownames(mats$latGraphs[[1]]$est)
  
  if (length(mats$latGraphs)==1) names <- '' else  names <- sapply(seq_along(mats$latGraphs), function(i)paste0(colnames(object$ClusterInfo$clusterDF),': ', object$ClusterInfo$clusterDF[i,], collapse = " - "))
  
  Table <- do.call(rbind,lapply(seq_along(mats$latGraphs), function(i){
    x <- mats$latGraphs[[i]]
    
    data.frame(
      Cluster = names[i],
      From = rep(Vars, times = length(Vars)),
      To = rep(Vars, each = length(Vars)),
      est = round(c(x$est),3),
      se = round(c(x$se),3),
      z = round(c(x$z),3),
      p.value = c(x$pvalue)
    )}))
  
  Table$p.value.Bonf <- pmin(Table[['p.value']] * nrow(Table), 1)
  
  Table$significance <- sapply(sapply(Table$p.value.Bonf,function(x){
    x[is.na(x)] <- 1
    sum(x < c(0.001,0.01,0.05))+1
  }),switch,'','*','**','***')
  
  Table$p.value.Bonf <- round(Table$p.value.Bonf,3)
  Table$p.value <- round(Table$p.value,3)
  
  #     Table[is.na(Table)] <- ''
  return(Table)
}
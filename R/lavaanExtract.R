lavextractNetwork <- function(fit, net, Ni)
{
  Res <- list()
  pars <- parameterEstimates(fit)
  if (!is.null(pars$group)) pars <- pars[pars$group==1,]
  if (missing(net)) net <- sort(unique(as.numeric(gsub("(^N)|(\\_.*$)","",pars$lhs[pars$op=="=~"]))))
  if (missing(Ni)) Ni <- sqrt(sum(grepl("^N1\\_",pars$lhs) & grepl("=~",pars$op)))
  
  for (i in net)
  {
    Res[[i]] <- matrix( pars[grepl(paste0("N",i,"_"),pars$lhs) & pars$op == "=~", "est"], Ni, Ni , byrow = TRUE)
    rownames(Res[[i]]) <- colnames(Res[[i]]) <- lavNames(fit,'ov.nox')
  }
  
  return(Res)
}

lavextractNu <- function(fit, nNet)
{
  pars <- parameterEstimates(fit)
  
  if (missing(nNet)) nNet <- max(unique(as.numeric(gsub("(^N)|(\\_.*$)","",pars$lhs[pars$op=="=~"]))))
  
  pars <- pars[pars$op == '~',]
  pars <- pars[ !duplicated(pars$label), ]
  
  Nu <- matrix(pars$est, , nNet, byrow = TRUE)
  
  return(Nu)
}


extractParFit <- function(fit, net, Ni)
{
  FullRes <- list()
  
  Res <- list()
  pars <- parameterEstimates(fit)
  if (!is.null(pars$group)) pars <- pars[pars$group==1,]
  if (missing(net)) net <- sort(unique(as.numeric(gsub("(^N)|(\\_.*$)","",pars$lhs[pars$op=="=~"]))))
  if (missing(Ni)) Ni <- sqrt(sum(grepl("^N1\\_",pars$lhs) & grepl("=~",pars$op)))
  
  for (i in net)
  {
    Res[[i]] <- list()
    
    Res[[i]]$est <- matrix( pars[grepl(paste0("N",i,"_"),pars$lhs) & pars$op == "=~", "est"], Ni, Ni , byrow = TRUE)
    Res[[i]]$se <- matrix( pars[grepl(paste0("N",i,"_"),pars$lhs) & pars$op == "=~", "se"], Ni, Ni , byrow = TRUE)
    Res[[i]]$z <- matrix( pars[grepl(paste0("N",i,"_"),pars$lhs) & pars$op == "=~", "z"], Ni, Ni , byrow = TRUE)
    Res[[i]]$pvalue <- matrix( pars[grepl(paste0("N",i,"_"),pars$lhs) & pars$op == "=~", "pvalue"], Ni, Ni , byrow = TRUE)
    
    for (j in 1:length(Res[[i]])) rownames(Res[[i]][[j]]) <- colnames(Res[[i]][[j]]) <- lavNames(fit,'ov.nox')
  }
  
  FullRes$latGraphs <- Res

  nNet <- max(net)
  pars <- parameterEstimates(fit)
  
  if (missing(nNet)) nNet <- max(unique(as.numeric(gsub("(^N)|(\\_.*$)","",pars$lhs[pars$op=="=~"]))))
  
  pars <- pars[pars$op == '~',]
  pars <- pars[ !duplicated(pars$label), ]
  
  Nu <- list(
    est = matrix(pars$est, , nNet, byrow = TRUE),
    se = matrix(pars$est, , nNet, byrow = TRUE),
    z = matrix(pars$est, , nNet, byrow = TRUE),
    pvalue = matrix(pars$est, , nNet, byrow = TRUE))
  
  FullRes$Nu <- Nu
  
  return(FullRes)
}

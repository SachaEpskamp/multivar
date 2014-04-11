suVARglm <-
function(x,family,parallel=FALSE,fit = c("chi2","aic","bic"),k=2,alpha=0.05)
{
  fit <- fit[1]
  
  if (missing(x)) stop("'x' must be assigned")
  x <- as.matrix(x)
  
  Ni <- ncol(x)
  Nt <- nrow(x)  
  
  if (missing(family)) 
  {
    if (identical(c(0,1),sort(unique(c(x))))) family <- rep("binomial",Ni) else family <- rep("gaussian",Ni)
  }
  if (length(family)==1)
  {
    family <- list(family)
    if (Ni > 1) for (i in 2:Ni) family[[i]] <- family[[1]]
  }
  if (length(family)!=Ni) stop("Length of family is not equal to number of variables.")
  
  # Stepup per node:
  if (isTRUE(parallel)) 
  {
    library("parallel")
    Res <- mclapply(seq_len(Ni),function(i)suVARglminner(x,i,family,fit,k,alpha), mc.cores=getOption("mc.cores", detectCores()))
  } else {
    Res <- lapply(seq_len(Ni),function(i)suVARglminner(x,i,family,fit,k,alpha))
  }  
  
  Out <- list(
    adjacency = as.matrix(do.call(cbind,lapply(Res,'[[','edges'))),
    graph = as.matrix(do.call(cbind,lapply(Res,'[[','estimates'))),
    history = lapply(Res,'[[','history')
  )
  
  return(Out)
}

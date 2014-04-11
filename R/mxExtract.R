mxextractNetwork <- function(fit, net, Ni)
{
  Ni <- length(fit@submodels[[1]]@manifestVars)/2
  nNet <- nrow(fit@submodels[[1]]@matrices$A@values) / Ni - 2
  
  Res <- list()
  for (i in 1:nNet)
  {
    Res[[i]] <- t(fit@submodels[[1]]@matrices$A@values[(Ni+1):(Ni*2), (Ni*2+1 + (i-1)*Ni) : (Ni*2 + i*Ni)  ])
  }
  return(Res)
}

mxextractNu <- function(fit, nNet)
{
  Ni <- sapply(seq_along(fit@submodels),   function(i) length(fit@submodels[[i]]@manifestVars)/2)
  nNet <- sapply(seq_along(fit@submodels),   function(i) nrow(fit@submodels[[i]]@matrices$A@values) / Ni[i] - 2)
  
  Nu <- lapply(seq_along(fit@submodels), function(i)  fit@submodels[[i]]@matrices$A@values[ seq(2*Ni[i] + 1, (2 + nNet[i]) * Ni[i], by = Ni[i]) ,1] )
  Nu <- do.call(rbind, Nu)
  
  return(Nu)
}

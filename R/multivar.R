# This function returns an output list of class 'multivar' with following ekements:
# latGraphs: latent networks
# Nu: The nu parameter
# persGraphs: personal networks
# Results: Results object

# Run multi-VAR analysis:
multivar <- function( 
  Data, # List of data per subject
  Ng = 1, # Number of graphs to extract
  pkg = 'lavaan', # Package to be used.
  ... # Arguments sent to GenerateMixVAR
)
{
#   # Add variance if missing:
#   for (i in seq_along(Data))
#   {
#     for (j in 1:ncol(Data[[i]]))
#     {
#      if (sd(Data[[i]][,j])==0) 
#       {
#        Data[[i]][,j] <- Data[[i]][,j] + rnorm(nrow(Data[[i]]),0,1e-6)
#       }
#     }
#   }
  
  # Generate model:
  mod <- GenerateMixVAR(Data, Ng, pkg, ...)
  
  # Output list:
  Res <- list(
    latGraphs = NULL,
    Nu = NULL,
    persGraphs = NULL,
    Results = NULL
    )
  
  # Run analysis:
  if (grepl("lavaan",pkg,ignore.case=TRUE))
  {  
#     Res$Results <- lavaan(mod$Model, sample.cov = lapply(mod$Data, cov,use = "pairwise.complete.obs"), sample.nobs=lapply(mod$Data, nrow), meanstructure=FALSE, 
#                           missing = "fiml")
    
    
    for (i in seq_along(mod$Data)) 
    {
      mod$Data[[i]] <- as.data.frame(mod$Data[[i]])
      mod$Data[[i]]$GroupID <- i
    }
    Res$Results <- lavaan(mod$Model, data = do.call(rbind,mod$Data), group = "GroupID", meanstructure=FALSE,
                          estimator = "GLS")
    
    Res$latGraphs <- lavextractNetwork(Res$Results)
    Res$Nu <- lavextractNu(Res$Results)
    Res$persGraphs <- lapply(1:nrow(Res$Nu),function(p)Reduce("+", mapply(Nu=Res$Nu[p,],G=Res$latGraphs,FUN=function(Nu,G)Nu*G,SIMPLIFY=FALSE)))
  }
  
  # Run analysis:
  if (grepl("mx",pkg,ignore.case=TRUE))
  {
    Res$Results <- mxRun(mod$Model)
    
    Res$latGraphs <- mxextractNetwork(Res$Results)
    Res$Nu <- mxextractNu(Res$Results)
    Res$persGraphs <- lapply(1:nrow(Res$Nu),function(p)Reduce("+", mapply(Nu=Res$Nu[p,],G=Res$latGraphs,FUN=function(Nu,G)Nu*G,SIMPLIFY=FALSE)))
  }

  class(Res) <- 'multivar'
  return(Res)
}
### User friendly function for estimating group different VAR networks:
# Data must be in long format

groupVAR <- function(
  data, # dataset, in long format
  vars, # vector of variables to use in analysis
  id.var, # subject id variable, as character
  group.vars, # Group variable or vector of group variables, if missing null model (one network) is fitted instead
  time.var, # Measurment variable, if missing data is assumed to be ordered
  method = c("LVCA","MultiLevel"),
  ... # Send to multivar
)
{

  ### add missing:
  stopifnot(!missing(vars))
  stopifnot(!missing(data))
  stopifnot(!missing(id.var))
  
  data <- as.data.frame(data)
  
  if (missing(group.vars) || is.null(group.vars))
  {
    data[['GROUPVAR']] <- 1
    group.vars <- 'GROUPVAR'
  }
  
  if (missing(time.var))
  {
    data[['TIMEVAR']] <- ave(data[[id.var]],data[[id.var]],FUN=seq_along)
    time.var <- 'TIMEVAR'
  }
  
  ## Group membership:
  GroupDF <- data[!duplicated(data[[id.var]]),group.vars,drop=FALSE]
  unGroups <- unique(GroupDF)
  
  Np <- length(unique(data[[id.var]]))
  
  membership <- integer(Np)
  for (i in seq_along(membership))
  {
    for (j in seq_len(nrow(unGroups)))
    {
      if (all(GroupDF[i,] == unGroups[j,]))
      {
        membership[i] <- j
        break
      }
    }
  }
  
  data <- data[order(data[[id.var]],data[[time.var]]),]
  
  # Check data:
  if (any(tapply(data[[time.var]], data[[id.var]], function(x) any(diff(sort(x)) > 1)))) stop("Skip in measurments found. This is not yet supported.")
  if (any(tapply(data[[id.var]], data[[id.var]], function(x) any(diff(sort(x)) > 1)))) stop("Skip in id found. This is not yet supported.")
  
  if (method[[1]] == "LVCA")
  {
    ## Split data per group:
    splData <- split(data[,vars,drop=FALSE], data[[id.var]])
    Np <- length(splData)

    # Run NULL model (1 network)
    Res0 <- multivar(splData, 1 ,pkg = 'lavaan', nu.groups = rep(1, Np), ...)
    
    # Run alternative model:
    ResA <- multivar(splData, max(membership), pkg = 'lavaan', nu.groups = membership, ...)
    
    # Store results:
    
    RetVal <- list(Model = ResA, Null = Res0, ClusterInfo = list(membership = membership, clusterDF = unGroups))
    class(RetVal) <- c('groupVAR','list')
    
    return(RetVal)
    
  } else if (method[[1]] == "MultiLevel")
  {
    
    stop("Not yet supported")
    
    stopifnot(require("lme4"))
    
    # Add membership to data:
    data[['GroupMembership']] <- as.factor(membership[data[[id.var]]])
    
    # Create exogenous and endogenous dataset (exo has last measure removed, endo first):
    exoDat <- data[as.logical(ave(data[[time.var]],data[[id.var]],FUN=function(x) x!=max(x) )) , vars]
    endoDat <- data[as.logical(ave(data[[time.var]],data[[id.var]],FUN=function(x) x!=min(x) )) , c(vars,id.var,'GroupMembership')]
    
    exoVars <- paste0("Exo_",vars)
    names(exoDat) <- exoVars
    
    fullData <- cbind(exoDat,endoDat)
    
    lmerRes <- list()
    for (v in seq_along(vars))
    {
      formula <- as.formula(paste0(vars[v],' ~ ',paste0('(',vars,'| GroupMembership )', collapse = '+')))
      lmerRes[[v]] <- lmer(formula, fullData)  
    }
    
    
  
  } else stop("Method not supported")

}



#### METHODS ####
plot.groupVAR <- function(object, panel = TRUE, values = c('est','sig'), sigbonf = FALSE, ...)
{
  stopifnot(require("qgraph"))
  Ng <- max(object$ClusterInfo$membership)
  if (panel)
  {
    layout(t(seq_len(Ng)))
  }
  
  Graph <- qgraph(object$Null$latGraphs[[1]], weighted = TRUE, DoNotPlot = TRUE, ...)
  
  Res <- list()
  if (values[[1]]=='est')
  {
    for (g in seq_len(Ng))
    {
      Res[[g]] <- qgraph(object$Model$latGraphs[[g]], weighted = TRUE, maximum = max(sapply(object$Model$latGraphs,function(x)max(abs(x)))), layout = Graph$layout, ...)
      if (Ng > 1) title(paste0(colnames(object$ClusterInfo$clusterDF),': ', object$ClusterInfo$clusterDF[g,], collapse = " - "), adj = 0, line = 3)
    }    
  } else if (values[[1]]=='sig')
  {

    fit <- extractParFit(object$Model$Results)
    
    for (g in seq_len(Ng))
    {
      Res[[g]] <- qgraph(fit$latGraphs[[g]]$pvalue + 1e-14, weighted = TRUE, layout = Graph$layout, mode = 'sig', bonf = sigbonf, ...)
      if (Ng > 1) title(paste0(colnames(object$ClusterInfo$clusterDF),': ', object$ClusterInfo$clusterDF[g,], collapse = " - "), adj = 0, line = 3)
    }  
    
  } else stop("Wrong 'graph' argument.")

  
  invisible(Res)
}


print.groupVAR <- function(object)
{
  Ng <- max(object$ClusterInfo$membership)
  
  if (Ng == 1)
  {
    fitTable <- data.frame(SingleNetwork = fitMeasures(object$Model$Result) )
  } else {
    fitTable <- cbind(SingleNetwork = round(fitMeasures(object$Null$Result),3),
                     MultipleNetworks = round(fitMeasures(object$Model$Result),3))
  
    # ANOVA #
    SingleNetwork <- object$Null$Result
    MultipleNetworks <- object$Model$Result
    
    aov <- anova(SingleNetwork, MultipleNetworks)
    
    cat('Difference test single network - multiple network models:\n\n')
    print(aov)
  }
  
  cat('\n\nFit measures:\n')
  print(fitTable)  
}


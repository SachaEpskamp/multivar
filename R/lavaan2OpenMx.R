lavaan2OpenMx <- function( lav, Data, model = "OpenMx Model")
{
  stopifnot(require("lavaan"))
  stopifnot(require("OpenMx"))
  
  # Check lav model:
  # Fix eq.id to zero (not sure why)
  lav$eq.id <- 0
  lav <- as.data.frame(lavaan:::lav2check(lav))
  
  # Add constraints:
  lav$lbound <- -Inf
  lav$ubound <- Inf
  
  if (any(lav$op == ">"))
  {
    for (l in which(lav$op==">"))
    {
      lav$lbound[lav$label == lav$lhs[l]] <- as.numeric(as.character(lav$rhs[l])) 
    }
  }
  
  if (any(lav$op == "<"))
  {
    for (l in which(lav$op=="<"))
    {
      lav$lbound[lav$label == lav$lhs[l]] <- as.numeric(as.character(lav$rhs[l])) 
    }
  }
  
  ## Set NA for empty labels:
  lav$label[lav$label == ''] <- NA
  
  ### RUN PER MODEL ###
  GroupModels <- list()
  if (is.null(names(Data))) names(Data) <- paste0("g",1:length(Data))
  GroupNames <- names(Data)
  
  NGroup <- max(lav$group)
  
  for (g in 1:NGroup)
  {

    # Extract labels:
    Labels <- sort(unique(c(lav[lav$group == g, ]$lhs[lav[lav$group == g, ]$op%in%c('=~','~','~~')],lav[lav$group == g, ]$rhs[lav[lav$group == g, ]$op%in%c('=~','~','~~')])))
    latLabs <- Labels[Labels %in% lav[lav$group == g, ]$lhs[lav[lav$group == g, ]$op == '=~']]
    manLabs <- Labels[! Labels %in% latLabs ]
    
    # Extract data:
    Cov <- cov(Data[[g]])
    Nobs <- nrow(Data[[g]])
    
    ### Generate model:
    ModelList <-  list(GroupNames[[g]], 
                       
                       # Data
                       type="RAM",
                       mxData(
                         observed=Cov, 
                         type="cov",
                         numObs = Nobs
                       ),
                       manifestVars = manLabs,
                       latentVars = latLabs)
    
    
    # Add factor loadings:
    if (any(lav[lav$group == g, ]$op=='=~'))
    {
      ModelList[[length(ModelList) + 1]] <- mxPath(
        from = lav[lav$group == g, ]$lhs[lav[lav$group == g, ]$op == '=~'],
        to = lav[lav$group == g, ]$rhs[lav[lav$group == g, ]$op == '=~'],
        arrows = 1,
        free = lav[lav$group == g, ]$free[lav[lav$group == g, ]$op == '=~']!=0,
        values = lav[lav$group == g, ]$ustart[lav[lav$group == g, ]$op == '=~'],
        labels = lav[lav$group == g, ]$label[lav[lav$group == g, ]$op == '=~'],
        lbound = lav[lav$group == g, ]$lbound[lav[lav$group == g, ]$op == '=~'],
        ubound = lav[lav$group == g, ]$ubound[lav[lav$group == g, ]$op == '=~']
      ) 
    }
    
    # Add regressions:
    if (any(lav[lav$group == g, ]$op=='~'))
    {
      ModelList[[length(ModelList) + 1]] <- mxPath(
        from = lav[lav$group == g, ]$rhs[lav[lav$group == g, ]$op == '~'],
        to = lav[lav$group == g, ]$lhs[lav[lav$group == g, ]$op == '~'],
        arrows = 1,
        free = lav[lav$group == g, ]$free[lav[lav$group == g, ]$op == '~']!=0,
        values = lav[lav$group == g, ]$ustart[lav[lav$group == g, ]$op == '~'],
        labels = lav[lav$group == g, ]$label[lav[lav$group == g, ]$op == '~'],
        lbound = lav[lav$group == g, ]$lbound[lav[lav$group == g, ]$op == '~'],
        ubound = lav[lav$group == g, ]$ubound[lav[lav$group == g, ]$op == '~']
      ) 
    }
    
    # Add variances
    if (any(lav[lav$group == g, ]$op=='~~'))
    {
      ModelList[[length(ModelList) + 1]] <- mxPath(
        from = lav[lav$group == g, ]$lhs[lav[lav$group == g, ]$op == '~~'],
        to = lav[lav$group == g, ]$rhs[lav[lav$group == g, ]$op == '~~'],
        arrows = 2,
        free = lav[lav$group == g, ]$free[lav[lav$group == g, ]$op == '~~']!=0,
        values = lav[lav$group == g, ]$ustart[lav[lav$group == g, ]$op == '~~'],
        labels = lav[lav$group == g, ]$label[lav[lav$group == g, ]$op == '~~'],
        lbound = lav[lav$group == g, ]$lbound[lav[lav$group == g, ]$op == '~~'],
        ubound = lav[lav$group == g, ]$ubound[lav[lav$group == g, ]$op == '~~']
      ) 
    }
    
    # Create model and return:
    GroupModels[[g]] <- do.call(mxModel, ModelList)
  }

  Objective <- eval(parse(text=paste( "mxAlgebra(",paste0(GroupNames,".objective",collapse = " + "),", name='obj')")))
  
  FullModel <- do.call(mxModel,c(
    list(model),
    GroupModels,
    list(Objective,mxAlgebraObjective("obj"))))
  
  return(FullModel)
}
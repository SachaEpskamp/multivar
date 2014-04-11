
# nu.groups: factor indicating group membership, or NA indicating free to estimate.

GenerateMixVAR <- function(
  Data, # List of dataframes to be analyse
  nNet = 1, # Number of latent networks.
  package = 'lavaan', # Package to be used
  nu, # Nu setup. NA if free, or value. If missing will be fixed automatically
  auto.fix.nu = missing(nu), # fix nu, only if nu is missing.
  auto.fix.gamma = FALSE, 
  bounds.nu = c(0, Inf), 
  adjacency)
{
  stopifnot(require("lavaan"))
  stopifnot(require("plyr"))
  options(stringsAsFactors = FALSE)
  
  Ni <- ncol(Data[[1]])
  if (missing(adjacency))
  {
    adjacency <- matrix(1,Ni,Ni)
  }
  adjacency[is.na(adjacency)] <- 0
  if (!all(unique(c(adjacency)) %in% c(0,1))) stop("Adjacency must contain only 1 (free) and 0.")
  
  # Data augmentation function:
  datAug <- function(x) 
  {
    if (is.null(colnames(x)))
    {
      Names <- paste0("y",1:ncol(x))
    } else Names <- colnames(x)
    
    Names <- c(paste0("x",Names), Names)
    
    x <- cbind(x[-nrow(x),], x[-1,])
    colnames(x) <- Names
    
    return(x)
  }
  
  # Prepare Data:
  Np <- length(Data)
  fullData <- lapply(Data, datAug)
  
  # Fix group membership:
#   if (missing(nu.groups))
#   {
#     nu.groups <- rep(NA, Np)
#   } else
#   {
#     nu.groups <- as.numeric(as.factor(nu.groups))
#   }

  
  # Variable names:
  varNames <- unique(gsub("^x","",c(unlist(sapply(fullData,colnames)))))
  
  # Network names:
  netNames <- paste0("N",1:nNet)
  

  ### NU SETUP ####
  if (missing(nu))
  {
    # k-means for starting values:
    VARcoef <- lapply(Data,function(x)c(VARglm(x, adjacency = adjacency)$graph))
    
    VARcoef <- do.call(rbind,VARcoef)
    if (nNet > Np) stop("Number of networks larger than number of persons.")
    if (nNet == 0) stop("0 graphs not supported.")
    if (nNet==1)
    {
      km <- list( 
        centers = matrix(colMeans(VARcoef),1,),
        cluster = rep(1,Np))
    } else if (nNet == Np)
    {
      km <- list( 
        centers = VARcoef,
        cluster = 1:Np)    
    } else km <- kmeans(VARcoef,nNet)
    kmeanGraphs <- alply(km$centers,1,matrix,nrow=Ni,ncol=Ni)
    kmClosest <- apply(km$centers,1,function(x) which.min(colSums(abs(t(VARcoef) - x))))
    
    nu <- matrix(NA,Np, nNet) 
    
    if (auto.fix.nu)
    {
      for (n in 1:nNet)
      {
        nu[ kmClosest[n], n] <- 1
      }
    }
  } else if (auto.fix.nu) warning("'auto.fix.nu' is not used if 'nu' is assigned.")



  # Dummy parameter table:
  dumPar <- list(
    id = integer(0),
    lhs = character(0),
    op = character(0),
    rhs = character(0),
    user = integer(0),
    group = integer(0),
    free = integer(0),
    ustart = numeric(0),
    exo = numeric(0),
    label = character(0),
    eq.id = numeric(0),
    unco = numeric(0)
  )
  
  # Full parameter table:
  parTab <- as.data.frame(dumPar,)
  
  
  # Helper functions:
  'addTo<-' <- function(x,value) 
  {
    x[(length(x)+1):(length(x)+length(value))] <- value
    return(x)
  }
  
  ParList2DF <- function(x)
  {
    l <- length(x$lhs)
    for (i in seq_along(x)) if (length(x[[i]])==0) x[[i]] <- NA
    return(as.data.frame(x))
  }
  
  
  ### Add nu weighting parameter ###
  # free over individuals, equal per network.
  ParTemp <- dumPar
  
  ParTemp$exo <- 0
  
  for (p in 1:Np)
  {
    
    for (n in 1:nNet)
    {
#       if (is.na(nu.groups[p]))
#       {
        addTo(ParTemp$user) <- rep(1,Ni)
        addTo(ParTemp$lhs) <- paste0(netNames[n],'_',varNames)
        addTo(ParTemp$op) <- rep('~',Ni)
        addTo(ParTemp$rhs) <- paste0('x',varNames)
        addTo(ParTemp$group) <- rep(p,Ni)
        addTo(ParTemp$label) <- rep(paste0("Nu_p",p,'_',n), Ni)
        
        if (!is.na(nu[p,n]))
          #         if (auto.fix.nu && kmClosest[n]==p)
        {
          addTo(ParTemp$ustart) <- rep(nu[p,n], Ni)
          addTo(ParTemp$free) <- rep(0,Ni)        
        } else 
        {
#           addTo(ParTemp$ustart) <- rep(0.1 + 0.8 * (km$cluster[p]==n), Ni)
          addTo(ParTemp$ustart) <- rep(0.5, Ni)
          addTo(ParTemp$free) <- rep(1,Ni)
          
          ## Constraints:
          addTo(ParTemp$user) <- 0
          addTo(ParTemp$lhs) <- paste0("Nu_p",p,'_',n)
          addTo(ParTemp$op) <- '>'
          addTo(ParTemp$rhs) <- bounds.nu[1]
          addTo(ParTemp$group) <- 0
          addTo(ParTemp$label) <- ''
          addTo(ParTemp$free) <- 0
          addTo(ParTemp$ustart) <- NA
          
          ## Constraints:
          addTo(ParTemp$user) <- 0
          addTo(ParTemp$lhs) <- paste0("Nu_p",p,'_',n)
          addTo(ParTemp$op) <- '<'
          addTo(ParTemp$rhs) <- bounds.nu[2]
          addTo(ParTemp$group) <- 0
          addTo(ParTemp$label) <- ''
          addTo(ParTemp$free) <- 0
          addTo(ParTemp$ustart) <- NA
        }
#       } else {
#         
#         addTo(ParTemp$user) <- rep(1, Ni)
#         addTo(ParTemp$lhs) <- paste0(netNames[n],'_',varNames)
#         addTo(ParTemp$op) <- rep('~',Ni)
#         addTo(ParTemp$rhs) <- paste0('x',varNames)
#         addTo(ParTemp$group) <- rep(p,Ni)
#         addTo(ParTemp$label) <- rep(paste0("Nu_p",p,'_',n), Ni)
#         
#         addTo(ParTemp$ustart) <- rep(1*(nu.groups[p] == n),Ni)
#         addTo(ParTemp$free) <- rep(0,Ni)        
#       }
    }
  }
  
  # Add (incomplete) to par table:
  parTab <- rbind(parTab, ParList2DF(ParTemp))
  
  
  ### Network parameters ###
  # Fixed over individual
  # Starting value set to kmeans center
  # Constrain highest positive edge
  ParTemp <- dumPar
  
  ParTemp$user <- 1
  ParTemp$exo <- 0
  
  for (p in 1:Np)
  {
    
    for (n in 1:nNet)
    {
#       maxN <- which(kmeanGraphs[[n]] == max(kmeanGraphs[[n]]),arr.ind=TRUE)
      for (i in 1:Ni)
      {
#         if ( auto.fix.gamma && i == maxN[2])
#         {
#           addTo(ParTemp$free) <- ifelse(1:Ni != maxN[1] & adjacency[i,]==1, 1, 0)
#         } else addTo(ParTemp$free) <- ifelse(adjacency[i,]==1,1,0)
#         
        addTo(ParTemp$free) <- ifelse(adjacency[i,]==1,1,0)
        
        addTo(ParTemp$lhs) <- rep(paste0(netNames[n],'_',varNames[i]), Ni)
        addTo(ParTemp$op) <- rep('=~',Ni)
        addTo(ParTemp$rhs) <- varNames
        addTo(ParTemp$group) <- rep(p,Ni)
        addTo(ParTemp$label) <- paste0(netNames[n],'_',varNames[i],'_',varNames)
#         addTo(ParTemp$ustart) <- ifelse(adjacency[i,]==1,kmeanGraphs[[n]][i,],0)
        addTo(ParTemp$ustart) <- ifelse(adjacency[i,]==1,0.5,0)
        
        
      }
    }
  }
  
  # Add (incomplete) to par table:
  parTab <- rbind(parTab, ParList2DF(ParTemp))
  
  # Residuals:
  ParTemp <- dumPar
  
  ParTemp$user <- 1
  ParTemp$free <- 1
  ParTemp$exo <- 0
  
  for (p in 1:Np)
  {
    addTo(ParTemp$lhs) <- varNames
    addTo(ParTemp$op) <- rep('~~',Ni)
    addTo(ParTemp$rhs) <- varNames
    addTo(ParTemp$group) <- rep(p,Ni)
    addTo(ParTemp$label) <- paste0('Residual_p',p,'_',1:Ni)
  }
  
  # Add (incomplete) to par table:
  parTab <- rbind(parTab, ParList2DF(ParTemp))
  
  ## Ad exogenous variables with a little trick:
  exoComb <- combn(paste0('x',varNames), 2)
  for (i in 1:Np)
  {
    ParTemp <- dumPar
    ParTemp$lhs <- exoComb[2,]
    ParTemp$rhs <- exoComb[1,]
    ParTemp$op <- '~~'
    ParTemp$user <- 0
    ParTemp$group <- i
    ParTemp$free <- 0
    ParTemp$ustart <- cov(fullData[[i]][,seq_along(varNames)])[lower.tri(cov(fullData[[i]][,seq_along(varNames)]))]
    ParTemp$exo <- 1
    ParTemp$label <- ""
    ParTemp$eq.id <- 0
    ParTemp$unco <- 0
    
    # Add (incomplete) to par table:
    parTab <- rbind(parTab, ParList2DF(ParTemp))
    
    ParTemp <- dumPar
    ParTemp$lhs <- paste0('x',varNames)
    ParTemp$rhs <- paste0('x',varNames)
    ParTemp$op <- '~~'
    ParTemp$user <- 0
    ParTemp$group <- i
    ParTemp$free <- 0
    ParTemp$ustart <- diag(cov(fullData[[i]][,seq_along(varNames)]))
    ParTemp$exo <- 1
    ParTemp$label <- ""
    ParTemp$eq.id <- 0
    ParTemp$unco <- 0
    
    # Add (incomplete) to par table:
    parTab <- rbind(parTab, ParList2DF(ParTemp))
  }
  #   
  #   
  #   exo <- lavaanify(paste("FOO ~ ",paste('x',varNames,collapse = '+')))
  #   for (i in 1:Np)
  #   {
  #     exo$group <- i
  #     parTab <- rbind(parTab,exo[exo$exo==1,])
  #   }
  
  ### Fix other entries ###
  parTab$id <- 1:nrow(parTab)
  unLabs <- unique(parTab$label[parTab$label!='' & parTab$free!=0])
  parTab$free[parTab$free!=0] <- match(parTab$label[parTab$free!=0],unLabs)
  parTab$free[is.na(parTab$free)] <- 0  
  parTab$free[parTab$free!=0] <- as.numeric(as.factor(parTab$free[parTab$free!=0]))
  
  parTab$eq.id <- ave(parTab$id, parTab$free, FUN = function(x) rep(x[1],length(x)))
  parTab$eq.id[parTab$free==0] <- 0
  
  parTab$unco <- cumsum(parTab$free!=0)
  parTab$unco[parTab$free==0] <- 0
  
  parTab$user <- as.integer(parTab$user)
  parTab$group <- as.integer(parTab$group)
  parTab$free <- as.integer(parTab$free)
  parTab$exo <- as.integer(parTab$exo)
  parTab$eq.id <- as.integer(parTab$eq.id)
  parTab$unco <- as.integer(parTab$unco)
  
  # Remove fixed to zero pars:
#   parTab <- parTab[!(parTab$free==0 & parTab$ustart == 0),]
  
  ### Run analysis:
  samCov <- lapply(fullData, cov)
  names(samCov) <- names(fullData)
  obs <- lapply(fullData, nrow)
  names(obs) <- names(fullData)
  
  #   parTab <- parTab[order(parTab$group),]
  
  Results <- list(Data = fullData)
  
  # mplus:
  if (grepl("mplus",package,ignore.case=TRUE))
  {
    Results$Model <- lavaan:::lav2mplus(parTab, group.label=NULL)
    return(Results)
  } else if (grepl("lavaan",package,ignore.case=TRUE))
  {
    Results$Model <- parTab
    return(Results)
  } else if (grepl("mx",package,ignore.case=TRUE))
  {
    Results$Model <- lavaan2OpenMx(parTab, fullData)
    return(Results)    
  } else stop("package not yet supported.")
}

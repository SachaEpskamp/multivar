
simVARnetsRandom <- function( 
  n, # Number of graphs
  k, # Number of items
  range = c(0.4,0.8), # Uniform range of weights
  edgeProb = 0.2, # Probability of edge
  negProb = 0.5, # Probability edge is negative
  checkEigen = TRUE, # Check if all eigenvalues are below 1
  checkEigenMax = 1000 # Maximum number of repeats before breaking generation of graphs.
)
{
  Graphs <- lapply(1:n, function(x) sample(0:1,k^2,TRUE, prob = c(1 - edgeProb, edgeProb)) *  matrix( sample(c(-1,1), k, TRUE)*runif(k, range[1], range[2]), k, k))
  
  # Redo till eigen are correct:
  count <- 1
  if (checkEigen && any(!sapply(Graphs,function(x)all(abs(eigen(x)$values) < 1))))
  {
    repeat
    {
      repGraphs <- !sapply(Graphs,function(x)all(abs(eigen(x)$values) < 1))
      Graphs[repGraphs] <- lapply(repGraphs, function(x) sample(0:1,k^2,TRUE, prob = c(1 - edgeProb, edgeProb)) *  matrix( sample(c(-1,1), k, TRUE)*runif(k, range[1], range[2]), k, k))
      if (all(sapply(Graphs,function(x)all(abs(eigen(x)$values) < 1)))) 
      {
        break
      }
      count <- count + 1
      if (count > checkEigenMax) stop("Could not construct a graph with eigenvalues < 1")
    }
  }
  
  return(Graphs)
}

constructPersGraphs <- function(
  Np, # Number of persons
  Nt, # Number of time points
  latGraphs, # List contaning latent graphs
  ErrorSD = 0.1, # Error sd on the personal graphs.
  Nu, # Nu parameter (can be missing to simulate random, using below properties)
  NuMinorRange = c(0, 0.3), # Range or minor graphs
  NuMajorRange = c(0.4, 0.8), # Range of major graphs
  checkEigen = TRUE, # Check if all eigenvalues are below 1 (only when Nu is missing)
  checkEigenMax = 1000 # Maximum number of repeats before breaking generation of graphs.
)
{
  stopifnot(!missing(Np))
  stopifnot(!missing(Nt))
  stopifnot(!missing(latGraphs))
  
  # Number of graphs:
  Ng <- length(latGraphs)
  
  # Number of items:
  if (length(unique(sapply(latGraphs,nrow))) > 1) stop("Latent networks are not of same size.")
  Ni <- nrow(latGraphs[[1]])
  
  # If Nu is missing, simulate it:
  NuEst <- missing(Nu)
  if (NuEst)
  {
    MajorGraph <- sample(1:Ng,Np,TRUE)
    
    Nu <- matrix(runif(Np*Ng,NuMinorRange[1],NuMinorRange[2]),Np,Ng)
    for (i in 1:Np) 
    {
      Nu[i,MajorGraph[i]] <- runif(1,NuMajorRange[1],NuMajorRange[2])
      Nu[i,] <- Nu[i,]
    }  
  }
  
  # Personal graphs:
  PersGraphs <- lapply(1:(Np), function(p) Reduce("+", mapply(Nu=Nu[p,],G=latGraphs,FUN=function(Nu,G)Nu*G,SIMPLIFY=FALSE)))
  for (i in seq_along(PersGraphs)) PersGraphs[[i]] <- PersGraphs[[i]] + matrix(rnorm(Ni^2,0,ErrorSD),Ni,Ni)
  
  if (any(!sapply(PersGraphs,function(x)all(abs(eigen(x)$values) < 1))))
  {
    if (!NuEst)
    {
      warning("Graphs with eigenvalues > 1 detected. VAR might not be stable")
    } else if (checkEigen)
    {
      count <- 1
      repeat
      {
        repGraphs <- !sapply(PersGraphs,function(x)all(abs(eigen(x)$values) < 1))
        Nu[repGraphs,] <- runif(sum(repGraphs),NuMinorRange[1],NuMinorRange[2])
        for (i in which(repGraphs)) 
        {
          Nu[i,MajorGraph[i]] <- runif(1,NuMajorRange[1],NuMajorRange[2])
        }
        PersGraphs <- lapply(1:(Np), function(p) Reduce("+", mapply(Nu=Nu[p,],G=latGraphs,FUN=function(Nu,G)Nu*G,SIMPLIFY=FALSE)))
        if (all(sapply(PersGraphs,function(x)all(abs(eigen(x)$values) < 1)))) 
        {
          break
        }
        count <- count + 1
        if (count > checkEigenMax) stop("Could not construct a graph with eigenvalues < 1")
      }
    }
  }
  
  return(list(
    persGraphs = PersGraphs,
    latGraphs = latGraphs,
    Nu = Nu))
}



# Sampling function:
sampleVARs <- function(
  Graphs, # Either a single graph or list of graphs
  Nt, # Either a  single number or a vector per person
  init, # Initial state, can be omitted or matrix with row per person), defaults to zeros
  errorSD = 1 # standard deviation of error
)
{
  if (is.matrix(Graphs)) Graphs <- list(Graphs)
  Np <- length(Graphs)
  Nt <- rep(Nt,length=Np)
  
  if (any(sapply(Graphs,function(x) length(unique(dim(x))) > 1 ))) stop ("non-square graph detected.")
  
  Ni <- sapply(Graphs, ncol)
  
  if (missing(init)) 
  {
    init <- matrix(0,Np,Ni)
  }
  if (is.vector(init))
  {
    init <- do.call(rbind, lapply(1:Np, function(x) init ))
  }
  
  Res <- list()
  for (p in 1:Np)
  {
    Res[[p]] <- matrix(,Nt,Ni)
    Res[[p]][1,] <- init[p,]
    for (t in 2:(Nt[p]))
    {
      Res[[p]][t,] <- t(Graphs[[p]]) %*% Res[[p]][t-1,] + rnorm(Ni[p], 0, errorSD)
    }    
  }
  
  if (length(Res)==1) Res <- Res[[1]]
  return(Res)
}


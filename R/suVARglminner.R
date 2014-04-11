suVARglminner <-
function(x,var,family,fit,k,alpha)
{
  message(paste("Testing incoming edges on node",var))
  
  family <- family[[var]]
  
  Ni <- ncol(x)
  inEdges <- inEdgesEsts <- rep(0,Ni)
  if (fit=="chi2") 
  {
    curfit <- -1 * logLik(glm(x[-1,var] ~ NULL))
  } else if (fit=="aic")
  {
    curfit <- AIC(glm(x[-1,var] ~ NULL),k=k)
  } else if (fit=="bic")
  {
    curfit <- BIC(glm(x[-1,var] ~ NULL))
  } else stop("'fit' must be 'chi2', 'aic' or 'bic'")
  
  addEdge <- function(x,i) {
    x[i] <- 1
    x
  }
  
  EdgesHist <- list()
  EdgesHist[[1]] <- inEdges
  c <- 2
  
  repeat
  {
    propEdges <- which(inEdges==0)
    glmRes <- lapply(propEdges,function(i)glm(x[-1,var] ~ x[-nrow(x),as.logical(addEdge(inEdges,i))],family=family))
    
    if (fit=="chi2") 
    {
      propfit <- -1 * sapply(glmRes,logLik)
      if (!any(pchisq(2*curfit - 2*propfit, 1, lower.tail=FALSE) < alpha)) break
    } else if (fit=="aic")
    {
      propfit <- sapply(glmRes,AIC)
      if (all(propfit >= curfit)) break
    } else if (fit=="bic")
    {
      propfit <- sapply(glmRes,BIC)
      if (all(propfit >= curfit)) break
    } 
    
    inEdges[propEdges[which.min(propfit)]] <- 1
    EdgesHist[[c]] <- inEdges
    c <- c + 1
    inEdgesEsts[as.logical(inEdges)] <- coef(glmRes[[which.min(propfit)]])[-1]
    curfit <- min(propfit)
  }
  return(list(edges = inEdges, estimates = inEdgesEsts, fit = curfit, history = EdgesHist))
}

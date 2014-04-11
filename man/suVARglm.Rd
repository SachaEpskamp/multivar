\name{suVARglm}
\alias{suVARglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Step-up estimation of lag-1 VAR network.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
suVARglm(x, family, parallel = FALSE, fit = c("chi2", "aic", "bic"), k = 2, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{family}{
%%     ~~Describe \code{family} here~~
}
  \item{parallel}{
%%     ~~Describe \code{parallel} here~~
}
  \item{fit}{
%%     ~~Describe \code{fit} here~~
}
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, family, parallel = FALSE, fit = c("chi2", "aic", 
    "bic"), k = 2, alpha = 0.05) 
{
    fit <- fit[1]
    if (missing(x)) 
        stop("'x' must be assigned")
    x <- as.matrix(x)
    Ni <- ncol(x)
    Nt <- nrow(x)
    if (missing(family)) {
        if (identical(c(0, 1), sort(unique(c(x))))) 
            family <- rep("binomial", Ni)
        else family <- rep("gaussian", Ni)
    }
    if (length(family) == 1) {
        family <- list(family)
        if (Ni > 1) 
            for (i in 2:Ni) family[[i]] <- family[[1]]
    }
    if (length(family) != Ni) 
        stop("Length of family is not equal to number of variables.")
    if (isTRUE(parallel)) {
        library("parallel")
        Res <- mclapply(seq_len(Ni), function(i) suVARglminner(x, 
            i, family, fit, k, alpha), mc.cores = getOption("mc.cores", 
            detectCores()))
    }
    else {
        Res <- lapply(seq_len(Ni), function(i) suVARglminner(x, 
            i, family, fit, k, alpha))
    }
    Out <- list(adjacency = as.matrix(do.call(cbind, lapply(Res, 
        "[[", "edges"))), graph = as.matrix(do.call(cbind, lapply(Res, 
        "[[", "estimates"))), history = lapply(Res, "[[", "history"))
    return(Out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

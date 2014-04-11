\name{VARglm}
\alias{VARglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit lag-1 VAR network using glm
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VARglm(x, family, vars, adjacency, icfun = BIC, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{family}{
%%     ~~Describe \code{family} here~~
}
  \item{vars}{
%%     ~~Describe \code{vars} here~~
}
  \item{adjacency}{
%%     ~~Describe \code{adjacency} here~~
}
  \item{icfun}{
%%     ~~Describe \code{icfun} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (x, family, vars, adjacency, icfun = BIC, ...) 
{
    if (missing(x)) 
        stop("'x' must be assigned")
    x <- as.matrix(x)
    Ni <- ncol(x)
    Nt <- nrow(x)
    if (missing(vars)) 
        vars <- 1:Ni
    No <- length(vars)
    if (missing(adjacency)) 
        adjacency <- matrix(1, Ni, No)
    if (is.vector(adjacency)) 
        adjacency <- as.matrix(adjacency)
    if (!is.matrix(adjacency) && ncol(adjacency) != No && nrow(adjacency) != 
        Ni) 
        stop("'adjacency' must be square matrix with a row for each predictor and column for each outcome variable.")
    if (missing(family)) {
        if (identical(c(0, 1), sort(unique(c(x))))) 
            family <- rep("binomial", No)
        else family <- rep("gaussian", No)
    }
    if (length(family) == 1) {
        family <- list(family)
        if (No > 1) 
            for (i in 2:No) family[[i]] <- family[[1]]
    }
    if (length(family) != No) 
        stop("Length of family is not equal to number of outcome variables.")
    Out <- list()
    Out$graph <- matrix(0, Ni, No)
    Out$IC <- 0
    for (i in 1:No) {
        if (is.function(family[[i]])) 
            fam <- family[[i]]
        else fam <- get(family[[i]])
        if (any(as.logical(adjacency[, i]))) {
            Res <- glm(x[-1, vars[i]] ~ x[-nrow(x), as.logical(adjacency[, 
                i])], family = fam)
        }
        else {
            Res <- glm(x[-1, vars[i]] ~ NULL, family = fam)
        }
        Out$graph[as.logical(adjacency[, i]), i] <- coef(Res)[-1]
        Out$IC <- Out$IC + icfun(Res, ...)
    }
    return(Out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

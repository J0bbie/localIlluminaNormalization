\name{filterFun}
\alias{filterFun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filteres probes that have a low expression value, thresholds are given as parameters.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
filterFun(x.lumiRaw, x.lumiNorm, filter.Th, filter.dp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.lumiRaw}{
%%     ~~Describe \code{x.lumiRaw} here~~
}
  \item{x.lumiNorm}{
%%     ~~Describe \code{x.lumiNorm} here~~
}
  \item{filter.Th}{
%%     ~~Describe \code{filter.Th} here~~
}
  \item{filter.dp}{
%%     ~~Describe \code{filter.dp} here~~
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
function (x.lumiRaw, x.lumiNorm, filter.Th, filter.dp) 
{
    dataMatrix = exprs(x.lumiNorm)
    presentCount <- detectionCall(x.lumiRaw, Th = filter.Th)
    filtered.normData <- dataMatrix[presentCount > filter.dp, 
        ]
    rmProbes = length(featureNames(normData)) - nrow(filtered.normData)
    cat("\nCreating a filtered data table to speed up the processing and reduce false positives by removing the unexpressed probes genes:\n", 
        "Settings: Rows which do not have p-value < ", filter.Th, 
        " in at least ", filter.dp + 1, " column will be excluded.\n", 
        " - Normalized data table contains: ", length(featureNames(normData)), 
        " probe rows;\n", " - Filtered table contains: ", nrow(filtered.normData), 
        " probe rows;\n", " - Removed ", rmProbes, " unexpressed probe rows from original normalized data file.\n", 
        sep = "")
    return(filtered.normData)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

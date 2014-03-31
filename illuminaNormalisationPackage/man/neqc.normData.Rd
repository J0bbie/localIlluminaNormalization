\name{neqc.normData}
\alias{neqc.normData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a lumiBatch R object containing the normalized expressions using neqc method. Uses the raw expressions and normalizes based on user/script input.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
neqc.normData(x.lumi, controlData)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.lumi}{
%%     ~~Describe \code{x.lumi} here~~
}
  \item{controlData}{
%%     ~~Describe \code{controlData} here~~
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
function (x.lumi, controlData) 
{
    print("Normalizing data using the 'limma::neqc function' using log2 and quantile normalization")
    controlData.neg <- controlData[grep("negative", tolower(controlData$controlType)), 
        ]
    controlData.neg <- controlData.neg[, 3:length(colnames(controlData.neg))]
    expData <- exprs(x.lumi)
    totalData <- rbind(expData, controlData.neg)
    status <- c(rep("regular", nrow(expData)), rep("negative", 
        nrow(controlData.neg)))
    normData.neqc <- limma::neqc(totalData, status)
    normData <- x.lumi
    exprs(normData) <- normData.neqc[1:nrow(normData.neqc), ]
    return(normData)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

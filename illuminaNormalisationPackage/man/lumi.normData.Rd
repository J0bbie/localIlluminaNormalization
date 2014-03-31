\name{lumi.normData}
\alias{lumi.normData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a lumiBatch R object containing the normalized expressions using lumi method. Uses the raw expressions and normalizes based on user/script input.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lumi.normData(x.lumi, bg.correct, bgCorrect.m, variance.stabilize, variance.m, normalize, normalization.m, normDataQC)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.lumi}{
%%     ~~Describe \code{x.lumi} here~~
}
  \item{bg.correct}{
%%     ~~Describe \code{bg.correct} here~~
}
  \item{bgCorrect.m}{
%%     ~~Describe \code{bgCorrect.m} here~~
}
  \item{variance.stabilize}{
%%     ~~Describe \code{variance.stabilize} here~~
}
  \item{variance.m}{
%%     ~~Describe \code{variance.m} here~~
}
  \item{normalize}{
%%     ~~Describe \code{normalize} here~~
}
  \item{normalization.m}{
%%     ~~Describe \code{normalization.m} here~~
}
  \item{normDataQC}{
%%     ~~Describe \code{normDataQC} here~~
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
function (x.lumi, bg.correct, bgCorrect.m, variance.stabilize, 
    variance.m, normalize, normalization.m, normDataQC) 
{
    normData = lumiExpresso(x.lumi, bg.correct = bg.correct, 
        bgcorrect.param = list(method = bgCorrect.m), variance.stabilize = variance.stabilize, 
        varianceStabilize.param = list(method = variance.m), 
        normalize = normalize, normalize.param = list(method = normalization.m), 
        QC.evaluation = normDataQC)
    return(normData)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{density.plot}
\alias{density.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Makes a density plot of the expressions.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
density.plot(x.lumi, ns, col = NULL, maxArray = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.lumi}{
%%     ~~Describe \code{x.lumi} here~~
}
  \item{ns}{
%%     ~~Describe \code{ns} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{maxArray}{
%%     ~~Describe \code{maxArray} here~~
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
function (x.lumi, ns, col = NULL, maxArray = NULL) 
{
    nrPages = ceiling(length(sampleNames(x.lumi))/maxArray)
    plotsPerPage = ceiling(length(sampleNames(x.lumi))/nrPages)
    pdf.file = paste(ns, "density.pdf", sep = "_")
    pdf(paper = "a4r", file = pdf.file, onefile = T, width = 18, 
        height = 12)
    for (p in 1:nrPages) {
        start = (p - 1) * plotsPerPage + 1
        end = min(p * plotsPerPage + (p > 1), length(sampleNames(x.lumi)))
        plot(x.lumi[, start:end], col = col[start:end], what = "density")
    }
    dev.off()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

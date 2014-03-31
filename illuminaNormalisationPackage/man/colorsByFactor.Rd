\name{colorsByFactor}
\alias{colorsByFactor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gets colors for the plots.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
colorsByFactor(experimentFactor)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{experimentFactor}{
%%     ~~Describe \code{experimentFactor} here~~
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
function (experimentFactor) 
{
    if (class(experimentFactor) != "factor") 
        stop("Parameter 'experimentFactor' must be of class 'factor'")
    if (length(levels(experimentFactor)) == 1) {
        plotColors <- rainbow(length(experimentFactor), s = 0.8, 
            v = 0.7)
        legendColors <- "white"
    }
    else {
        tab.tmp <- table(experimentFactor)
        colors.light <- rainbow(length(levels(experimentFactor)), 
            s = 1 - sapply(tab.tmp, min, 5) * 0.1)
        colors.dark <- rainbow(length(levels(experimentFactor)), 
            v = 1 - sapply(tab.tmp, min, 5) * 0.14)
        plotColors <- NULL
        legendColors <- NULL
        for (l in 1:length(levels(experimentFactor))) {
            colorFun <- colorRampPalette(c(colors.light[l], colors.dark[l]))
            tmpColors <- colorFun(tab.tmp[l])
            plotColors[experimentFactor == levels(experimentFactor)[l]] <- tmpColors
            legendColors[l] <- tmpColors[ceiling(length(tmpColors)/2)]
        }
    }
    return(list(plotColors = plotColors, legendColors = legendColors))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

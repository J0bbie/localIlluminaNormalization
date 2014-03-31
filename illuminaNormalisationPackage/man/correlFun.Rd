\name{correlFun}
\alias{correlFun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Makes a correlation plots of the raw/normed expressions.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
correlFun(x.lumi, normalized = FALSE, clusterOption1 = clusterOption1, clusterOption2 = clusterOption2, experimentFactor = NULL, legendColors = NULL, WIDTH = 1000, HEIGHT = 1414, POINTSIZE = 24, MAXARRAY = 41, normalization.m = NULL, fileName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.lumi}{
%%     ~~Describe \code{x.lumi} here~~
}
  \item{normalized}{
%%     ~~Describe \code{normalized} here~~
}
  \item{clusterOption1}{
%%     ~~Describe \code{clusterOption1} here~~
}
  \item{clusterOption2}{
%%     ~~Describe \code{clusterOption2} here~~
}
  \item{experimentFactor}{
%%     ~~Describe \code{experimentFactor} here~~
}
  \item{legendColors}{
%%     ~~Describe \code{legendColors} here~~
}
  \item{WIDTH}{
%%     ~~Describe \code{WIDTH} here~~
}
  \item{HEIGHT}{
%%     ~~Describe \code{HEIGHT} here~~
}
  \item{POINTSIZE}{
%%     ~~Describe \code{POINTSIZE} here~~
}
  \item{MAXARRAY}{
%%     ~~Describe \code{MAXARRAY} here~~
}
  \item{normalization.m}{
%%     ~~Describe \code{normalization.m} here~~
}
  \item{fileName}{
%%     ~~Describe \code{fileName} here~~
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
function (x.lumi, normalized = FALSE, clusterOption1 = clusterOption1, 
    clusterOption2 = clusterOption2, experimentFactor = NULL, 
    legendColors = NULL, WIDTH = 1000, HEIGHT = 1414, POINTSIZE = 24, 
    MAXARRAY = 41, normalization.m = NULL, fileName) 
{
    if (is.null(experimentFactor)) 
        stop("the 'exerimentFactor' parameter is required")
    if (is.null(legendColors)) 
        stop("the 'legendColors' parameter is required")
    if (normalized) {
        Type = "NORM"
        text1 <- paste("Array correlation plot\nafter", normalization.m, 
            "normalization")
    }
    else {
        Type <- "RAW"
        text1 <- "Raw data correlation plot"
    }
    if (length(sampleNames(x.lumi)) < 2) {
        warning("Only one array in dataset, no correlation plot made")
    }
    else {
        png(file = paste(fileName, Type, "dataArrayCorrelation.png", 
            sep = "_"), width = WIDTH, height = HEIGHT, pointsize = POINTSIZE)
        if (length(sampleNames(x.lumi)) < MAXARRAY) {
            par(oma = c(17, 0, 0, 0), cex.axis = 0.7, cex.main = 0.8)
        }
        else {
            par(oma = c(17, 0, 0, 0), srt = 90, las = 2, cex.axis = 0.5, 
                cex.main = 0.8)
        }
        COpt1 <- "pearson"
        if (tolower(clusterOption1) == "spearman") 
            COpt1 <- "spearman"
        crp <- cor(exprs(x.lumi), use = "complete.obs", method = COpt1)
        text1 <- paste(text1, "\ncorrelation method:", COpt1, 
            "\ncluster method:", clusterOption2)
        switch(tolower(clusterOption1), pearson = {
            my.dist <- function(x) cor.dist(x, abs = FALSE)
        }, spearman = {
            my.dist <- function(x) spearman.dist(x, abs = FALSE)
        }, euclidean = {
            my.dist <- function(x) euc(x)
        })
        my.hclust <- function(d) hclust(d, method = clusterOption2)
        sideColors <- legendColors[as.numeric(experimentFactor)]
        heatmap.2(crp, distfun = my.dist, hclustfun = my.hclust, 
            trace = "none", symm = TRUE, density.info = "density", 
            main = text1, dendrogram = "row", ColSideColors = sideColors, 
            ncol = 3)
        x <- dev.off()
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

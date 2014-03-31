\name{clusterFun}
\alias{clusterFun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Makes a clustering of the raw/normed expression data based on the groups defined in the description file.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
clusterFun(x.lumi, normalized = FALSE, experimentFactor = NULL, clusterOption1 = clusterOption1, clusterOption2 = clusterOption2, plotColors = NULL, legendColors = NULL, plotSymbols = NULL, legendSymbols = NULL, WIDTH = 1000, HEIGHT = 1414, POINTSIZE = 24, MAXARRAY = 41, normalization.m = NULL, fileName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x.lumi}{
%%     ~~Describe \code{x.lumi} here~~
}
  \item{normalized}{
%%     ~~Describe \code{normalized} here~~
}
  \item{experimentFactor}{
%%     ~~Describe \code{experimentFactor} here~~
}
  \item{clusterOption1}{
%%     ~~Describe \code{clusterOption1} here~~
}
  \item{clusterOption2}{
%%     ~~Describe \code{clusterOption2} here~~
}
  \item{plotColors}{
%%     ~~Describe \code{plotColors} here~~
}
  \item{legendColors}{
%%     ~~Describe \code{legendColors} here~~
}
  \item{plotSymbols}{
%%     ~~Describe \code{plotSymbols} here~~
}
  \item{legendSymbols}{
%%     ~~Describe \code{legendSymbols} here~~
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
function (x.lumi, normalized = FALSE, experimentFactor = NULL, 
    clusterOption1 = clusterOption1, clusterOption2 = clusterOption2, 
    plotColors = NULL, legendColors = NULL, plotSymbols = NULL, 
    legendSymbols = NULL, WIDTH = 1000, HEIGHT = 1414, POINTSIZE = 24, 
    MAXARRAY = 41, normalization.m = NULL, fileName) 
{
    if (is.null(experimentFactor)) 
        stop("The 'experimentFactor' parameter must be specified")
    if (is.null(plotColors)) 
        stop("the 'plotColors' parameter is required")
    if (is.null(legendColors)) 
        stop("the 'legendColors' parameter is required")
    if (is.null(plotSymbols)) 
        stop("the 'plotSymbols' parameter is required")
    if (is.null(legendSymbols)) 
        stop("the 'legendSymbols' parameter is required")
    if (normalized) {
        Type = "NORM"
        main <- paste("Cluster dendrogram of", normalization.m, 
            "normalized data")
    }
    else {
        if (normalized == FALSE) 
            Type <- "RAW"
        main <- "Cluster dendrogram of raw data"
    }
    if (length(sampleNames(x.lumi)) < 3) {
        warning("Only ", length(sampleNames(x.lumi)), " sample(s) in dataset, no clustering plot made")
    }
    else {
        switch(tolower(clusterOption1), pearson = {
            correl <- cor.dist(t(exprs(x.lumi)), abs = FALSE)
        }, spearman = {
            correl <- spearman.dist(t(exprs(x.lumi)), abs = FALSE)
        }, euclidean = {
            correl <- euc(t(exprs(x.lumi)))
        })
        clust <- hclust(correl, method = tolower(clusterOption2))
        png(file = paste(fileName, "_", Type, "_", "dataCluster_", 
            clusterOption1, "_", clusterOption2, ".png", sep = ""), 
            width = WIDTH, height = HEIGHT, pointsize = POINTSIZE)
        if (length(sampleNames(x.lumi)) < MAXARRAY) {
            cexval1 <- 0.75
            cexval2 <- 1.23
            cexval3 <- 0.55
        }
        else {
            cexval1 <- 0.55
            cexval2 <- 1.6
            cexval3 <- 0.41
        }
        par(cex = cexval1, oma = c(14, 1, 0, 0))
        par(cex.axis = cexval2, cex.lab = cexval2, cex.main = cexval2)
        plot(clust, hang = -1, main = main, xlab = paste("distance:", 
            clusterOption1), sub = paste(" cluster method:", 
            clusterOption2))
        points(1:length(clust$order), rep(0, length(clust$order)), 
            pch = 15, col = "white", cex = 1.5)
        points(1:length(clust$order), rep(0, length(clust$order)), 
            pch = plotSymbols[clust$order], col = plotColors[clust$order])
        if (length(levels(experimentFactor)) > 1) {
            legend("topright", levels(experimentFactor), pch = legendSymbols, 
                col = legendColors, ncol = 3)
        }
        par(cex = cexval3)
        dev.off()
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

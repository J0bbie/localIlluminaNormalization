\name{pcaFun}
\alias{pcaFun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a PCA of the expressions, grouped on the groups defined in the description file.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pcaFun(x.lumi, normalized = FALSE, experimentFactor = NULL, scaled_pca = TRUE, plotColors = NULL, legendColors = NULL, plotSymbols = NULL, legendSymbols = NULL, namesInPlot = FALSE, groupsInLegend = NULL, WIDTH = 1000, HEIGHT = 1414, POINTSIZE = 24, normalization.m = NULL, fileName)
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
  \item{scaled_pca}{
%%     ~~Describe \code{scaled_pca} here~~
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
  \item{namesInPlot}{
%%     ~~Describe \code{namesInPlot} here~~
}
  \item{groupsInLegend}{
%%     ~~Describe \code{groupsInLegend} here~~
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
    scaled_pca = TRUE, plotColors = NULL, legendColors = NULL, 
    plotSymbols = NULL, legendSymbols = NULL, namesInPlot = FALSE, 
    groupsInLegend = NULL, WIDTH = 1000, HEIGHT = 1414, POINTSIZE = 24, 
    normalization.m = NULL, fileName) 
{
    if (groupsInLegend & is.null(experimentFactor)) 
        stop("The 'experimentFactor' parameter must be specified when groupsInLegend is true")
    if (is.null(plotColors)) 
        stop("the 'plotColors' parameter is required")
    if (groupsInLegend & is.null(legendColors)) 
        stop("the 'legendColors' parameter is required when groupsInLegend is true")
    if (is.null(plotSymbols)) 
        stop("the 'plotSymbols' parameter is required")
    if (groupsInLegend & is.null(legendSymbols)) 
        stop("the 'legendSymbols' parameter is required when groupsInLegend is true")
    if (groupsInLegend & !is.null(experimentFactor)) {
        if (length(levels(experimentFactor)) <= 1) {
            warning("groupsInLegend set to true, but no groups indicated; groups not added to legend")
        }
    }
    if (length(sampleNames(x.lumi)) < 3) {
        warning("Only", length(sampleNames(x.lumi)), "sample(s) in dataset, no PCA plot made")
    }
    else {
        if (normalized) {
            Type = "NORM"
            tmain <- paste("PCA analysis of", normalization.m, 
                "normalized data", sep = " ")
        }
        else {
            Type <- "RAW"
            tmain <- "PCA analysis of raw data"
        }
        pca1 <- NULL
        try(pca1 <- prcomp(t(exprs(x.lumi)[apply(exprs(x.lumi), 
            1, function(r) {
                sum(is.na(r)) == 0
            }), ]), retx = T, center = T, scale = scaled_pca), 
            TRUE)
        if (is.null(pca1) & scaled_pca) {
            try(pca1 <- prcomp(t(exprs(Data)[apply(exprs(x.lumi), 
                1, function(r) {
                  sum(is.na(r)) == 0
                }), ]), retx = T, center = T, scale = FALSE), 
                TRUE)
            if (!is.null(pca1)) 
                warning("pca with scaling unsuccessful, successfully retried without scaling")
        }
        if (!is.null(pca1)) {
            perc_expl1 <- round(((pca1$sdev[1:3]^2)/sum(pca1$sdev^2)) * 
                100, 2)
            cex.circle <- 1.5
            cex.text <- 0.7
            cex.legend <- 0.75
            tcol <- "#444444"
            png(file = paste(fileName, Type, "dataPCA_analysis.png", 
                sep = "_"), width = WIDTH + 200 * (!namesInPlot), 
                height = HEIGHT + 283 * (!namesInPlot), pointsize = POINTSIZE)
            if (!namesInPlot) {
                layout(rbind(c(1, 1, 2, 2, 5), c(3, 3, 4, 4, 
                  5)))
            }
            else {
                layout(rbind(c(1, 1, 2, 2), c(1, 1, 2, 2), c(3, 
                  3, 4, 4), c(3, 3, 4, 4)))
            }
            par(oma = c(20, 0, 5, 0))
            plot(pca1$x[, 1], pca1$x[, 2], cex = cex.circle, 
                pch = plotSymbols, col = plotColors, xlab = paste("PC1 (", 
                  perc_expl1[1], "\%)", sep = ""), ylab = paste("PC2 (", 
                  perc_expl1[2], "\%)", sep = ""))
            if (namesInPlot) 
                text(pca1$x[, 1], pca1$x[, 2], sampleNames(x.lumi), 
                  pos = 4, cex = cex.text, col = tcol)
            plot(pca1$x[, 1], pca1$x[, 3], cex = cex.circle, 
                pch = plotSymbols, col = plotColors, xlab = paste("PC1 (", 
                  perc_expl1[1], "\%)", sep = ""), ylab = paste("PC3 (", 
                  perc_expl1[3], "\%)", sep = ""))
            if (namesInPlot) 
                text(pca1$x[, 1], pca1$x[, 3], sampleNames(x.lumi), 
                  pos = 4, cex = cex.text, col = tcol)
            plot(pca1$x[, 2], pca1$x[, 3], cex = cex.circle, 
                pch = plotSymbols, col = plotColors, xlab = paste("PC2 (", 
                  perc_expl1[2], "\%)", sep = ""), ylab = paste("PC3 (", 
                  perc_expl1[3], "\%)", sep = ""))
            if (namesInPlot) 
                text(pca1$x[, 2], pca1$x[, 3], sampleNames(x.lumi), 
                  pos = 4, cex = cex.text, col = tcol)
            barplot((100 * pca1$sdev^2)/sum(pca1$sdev^2), xlab = "components", 
                ylab = "\% of total variance explained")
            groupsAdded <- FALSE
            if (groupsInLegend) {
                if (length(levels(experimentFactor)) > 1) {
                  groupsAdded <- TRUE
                }
            }
            if (namesInPlot) {
                if (groupsAdded) {
                  legend("topright", levels(experimentFactor), 
                    pch = legendSymbols, col = legendColors, 
                    cex = cex.legend)
                }
            }
            else {
                par(mar = c(0, 0, 0, 0))
                plot(1, type = "n", xaxt = "n", yaxt = "n", xlab = "", 
                  ylab = "", bty = "n")
                if (groupsAdded) {
                  legend("topleft", c(levels(experimentFactor), 
                    "", sampleNames(x.lumi), ncol = 3), pch = c(legendSymbols, 
                    20, plotSymbols), col = c(legendColors, "white", 
                    plotColors), cex = (cex.legend + 0.1))
                }
                else {
                  legend("topleft", sampleNames(x.lumi), pch = plotSymbols, 
                    col = plotColors, cex = 0.7, bty = "n", ncol = 3)
                }
            }
            mtext(tmain, side = 3, outer = TRUE, font = 2, cex = 1.2)
            dev.off()
        }
        else {
            warning("PCA on the", Type, "data set unsuccessful, image skipped")
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

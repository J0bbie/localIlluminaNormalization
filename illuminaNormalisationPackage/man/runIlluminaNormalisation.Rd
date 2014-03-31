\name{runIlluminaNormalisation}
\alias{runIlluminaNormalisation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Main flow to run the entire normalization process. Gets the arguments from a CURL/command line arguments call.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
runIlluminaNormalisation(commandArgsCurl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{commandArgsCurl}{
%%     ~~Describe \code{commandArgsCurl} here~~
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
function (commandArgsCurl) 
{
    ptm <- proc.time()
    SCRIPT.DIR <- dirname(sys.frame(1)$ofile)
    source(paste(SCRIPT.DIR, "getArguments.R", sep = "/"))
    userParameters <- getArguments(commandArgsCurl)
    source(paste(userParameters$scriptDir, "functions_loadPackages.R", 
        sep = "/"))
    source(paste(userParameters$scriptDir, "functions_makeImages.R", 
        sep = "/"))
    source(paste(userParameters$scriptDir, "functions_qualityControl.R", 
        sep = "/"))
    cat("\nLoading required packages.\n")
    pkgs <- c("limma", "ALL", "bioDist", "gplots", "annotate", 
        "arrayQualityMetrics", switch(userParameters$species, 
            Human = pkgs <- c("org.Hs.eg.db"), Mouse = pkgs <- c("org.Mm.eg.db"), 
            Rat = pkgs <- c("org.Rn.eg.db")), "lumi", userParameters$lib.mapping, 
        userParameters$lib.All.mapping)
    loadPackages(pkgs)
    cat("\nRequired packages succesfully loaded.\n")
    if (userParameters$normalize) {
        expData <- paste(userParameters$inputDir, userParameters$sampleProbeProfilePath, 
            sep = "")
        cat("\nLoading sample probe profile:", expData, "\n", 
            sep = "")
        rawData <- import.rawData(expData, userParameters$detectionTh, 
            userParameters$convertNuID, userParameters$checkDupId, 
            userParameters$lib.mapping, userParameters$dec, userParameters$parseColumnName, 
            userParameters$rawDataQC)
        sampleNames(rawData) <- make.names(sampleNames(rawData))
        cat("\nSuccesfully loaded the Sample Probe Profile.\n")
        cat("\nChecking if description data is valid for the given sample probe profile.\n")
        description$arraySampleNames = make.names(description[, 
            1])
        if (sum(length(sampleNames(rawData)) - length(description[, 
            1])) > 0) {
            message <- paste("Error: Number of array names in raw data file and number of array names in description file is not of the same size!")
            cat(message)
            if (userParameters$createLog) 
                sink()
            stop(message)
        }
        file_order <- match(description[, 4], sampleNames(rawData))
        if (sum(is.na(file_order)) > 0) {
            message <- paste("\nError: Assigned array names in raw data file and file names in description file do not match!\n")
            cat(message)
            if (userParameters$createLog) 
                sink()
            stop(message)
        }
        if (length(description[, 2]) != length(unique(description[, 
            2]))) {
            message <- ("Error: Assigned sampleNames are not unique!")
            cat(message)
            if (userParameters$createLog) 
                sink()
            stop(message)
        }
        rawData <- rawData[, file_order]
        cat("\nDescription data is valid.\n")
        if (userParameters$perGroup) {
            cat("\nRe-ordering raw Sample Probe Profile per group defined in the description file.\n")
            file_order2 <- match(description2[, 4], sampleNames(rawData))
            if (sum(is.na(file_order2)) > 0) {
                message <- ("Error: File names in Sample Probe Profile and file names in description file do not match!")
                cat(message)
                if (userParameters$createLog) 
                  sink()
                stop(message)
            }
            rawData <- rawData[, file_order2]
            sampleNames(rawData) <- as.character(description2[, 
                2])
            cat("\nRe-ordering succesfull.\n")
        }
        else {
            sampleNames(rawData) <- as.character(description[, 
                2])
            cat("\nSample names have been given to the arrays.\n")
        }
        if (userParameters$bgSub) {
            cat("\nSkipping background correction\n", "\nNormalizing the raw Sample Probe Profiles:", 
                userParameters$sampleProbeProfilePath, "\n", 
                sep = "")
            normData <- lumi.normData(rawData, bg.correct = FALSE, 
                userParameters$bgcorrect.m, userParameters$variance.stabilize, 
                userParameters$variance.m, userParameters$normalize, 
                userParameters$normalization.m, userParameters$normDataQC)
            cat("\nNormalization of raw data has been successfull.\n")
        }
        else {
            controlData <- paste(userParameters$inputDir, userParameters$controlProbeProfilePath, 
                sep = "")
            cat("\nPerforming background correction (", userParameters$bgcorrect.m, 
                ") on the Sample Probe Profile using the Control Probe Profile: ", 
                controlData, "\n", sep = "")
            cat("\nLoading Control Probe Profile:", userParameters$controlProbeProfilePath, 
                "\n", sep = "")
            cat("\nCombining Control data with Sample data.\n")
            rawData.ctrl <- addControlData2lumi(controlData, 
                rawData)
            controlData <- as.data.frame(getControlData(rawData.ctrl), 
                row.names = NULL)
            cat("\nNormaling (", userParameters$normType, ") the raw Sample Probe Profile using background correction:", 
                userParameters$controlProbeProfilePath, "\n", 
                sep = "")
            switch(userParameters$normType, lumi = normData <- lumi.normData(rawData.ctrl, 
                bg.correct = TRUE, userParameters$bgcorrect.m, 
                userParameters$variance.stabilize, userParameters$variance.m, 
                userParameters$normalize, userParameters$normalization.m, 
                userParameters$normDataQC), neqc = normData <- neqc.normData(rawData.ctrl, 
                controlData))
            cat("\nNormalization of raw data with background correction has been successfull.\n")
        }
        cat("\nCreating eSets expression matrix of raw data.\n")
        eset.rawData <- exprs(rawData)
        cat("\nSuccesfully created eSets expression matrix of raw data.\n")
        cat("\nCreating eSets expression matrix of normalized data.\n")
        eset.normData <- exprs(normData)
        cat("\nSuccesfully created eSets expression matrix of normalized data.\n")
        if (userParameters$rawSummary) {
            fileName <- paste(userParameters$outputDir, userParameters$studyName, 
                "_summary_rawData.txt", sep = "")
            cat("\nCreating summary file of the means and SD of the raw data: ", 
                fileName, "\n", sep = "")
            rawSum.table = createSummary(rawData, fileName)
            cat("\nSuccesfully made summary file of raw data.\n")
        }
        if (userParameters$normSummary) {
            fileName <- paste(userParameters$outputDir, userParameters$studyName, 
                "_summary_normData.txt", sep = "")
            cat("\nCreating summary file of the means and SD of the normalized data: ", 
                fileName, "\n", sep = "")
            normSum.table = createSummary(normData, fileName)
            cat("\nSuccesfully made summary file of normalized data.\n")
        }
        if (userParameters$save.rawData) {
            fileName <- paste(userParameters$outputDir, userParameters$studyName, 
                "_rawData.Rdata", sep = "")
            cat("\nSaving lumiBatch R object of the raw data in: ", 
                fileName, "\n", sep = "")
            save(rawData, file = fileName)
            cat("\nSuccesfully saved lumiBatch R object of the raw data\n")
        }
        if (userParameters$save.normData) {
            fileName <- paste(userParameters$outputDir, userParameters$studyName, 
                "_normData.Rdata", sep = "")
            cat("\nSaving lumiBatch R object of the normalized data in: ", 
                fileName, "\n", sep = "")
            save(normData, file = fileName)
            cat("\nSuccesfully saved lumiBatch R object of the normalized data\n")
        }
    }
    if (userParameters$performStatistics) {
        if (userParameters$statSubset) {
            cat(paste("\nReading in the statFile conating the subset of samples on which to perform the statistics: ,", 
                paste(userParameters$inputDir, userParameters$statFile, 
                  sep = ""), "\n", sep = ""))
            statFile <- read.table(paste(userParameters$inputDir, 
                userParameters$statFile, sep = "/"), header = F, 
                stringsAsFactors = F, sep = "\t", quote = "")
            statFile$arraySampleNames = make.names(statFile[, 
                1])
            statFile$groups <- description$FactorValue[match(statFile[, 
                1], description$SourceName)]
            cat("\nSuccesfully read in the statFile containing the samples for the subset!\n")
            if (userParameters$perGroup) {
                statFile = statFile[order(statFile$groups, statFile$arraySampleNames), 
                  ]
            }
        }
        cat("\nCreating a plot colorset for each array group.\n")
        if (userParameters$statSubset == FALSE) {
            if (userParameters$perGroup) {
                experimentFactor <- as.factor(description2[, 
                  3])
                colList <- colorsByFactor(experimentFactor)
                plotColors <- colList$plotColors
                legendColors <- colList$legendColors
                rm(colList)
            }
            else {
                experimentFactor <- as.factor(description[, 3])
                colList <- colorsByFactor(experimentFactor)
                plotColors <- colList$plotColors
                legendColors <- colList$legendColors
                rm(colList)
            }
        }
        else {
            experimentFactor <- as.factor(statFile$groups)
            colList <- colorsByFactor(experimentFactor)
            plotColors <- colList$plotColors
            legendColors <- colList$legendColors
            rm(colList)
        }
        plotSymbols <- 18 - as.numeric(experimentFactor)
        legendSymbols <- sort(unique(plotSymbols), decreasing = TRUE)
        cat("\nPlot colorset sucesfully made.\n")
        if (userParameters$rawDataQC) {
            if (userParameters$statSubset) {
                cat("\nMaking subset of samples in raw data.!\n")
                x <- sampleNames(rawData)[order(sampleNames(rawData))]
                matchedSamples <- match(statFile[, 1], x)
                rawData <- rawData[, matchedSamples]
                cat("\nSuccesfully made subset of samples in raw data!\n")
            }
            cat("\nCreating QC plots for the raw data.\n")
            fileNamePrefix <- paste(userParameters$outputDir, 
                "/", userParameters$studyName, "_RAW", sep = "")
            if (userParameters$raw.boxplot) {
                cat("\nPlot boxplot for raw intensities\n")
                gar <- box.plot(rawData, fileNamePrefix, col = plotColors, 
                  maxArray = 50)
            }
            if (userParameters$raw.density) {
                cat("\nPlot density histogram for raw intensities\n")
                gar <- density.plot(rawData, fileNamePrefix, 
                  col = plotColors, maxArray = 16)
            }
            if (userParameters$raw.cv) {
                cat("\nPlot density for coefficient of variance for raw intensities\n")
                gar <- cv.plot(rawData, fileNamePrefix, col = plotColors, 
                  maxArray = 16)
            }
            fileNamePrefix <- paste(userParameters$outputDir, 
                "/", userParameters$studyName, sep = "")
            if (userParameters$raw.sampleRelation) {
                cat("\nHierarchical clustering of raw data\n")
                gar <- clusterFun(rawData, normalized = FALSE, 
                  experimentFactor = experimentFactor, clusterOption1 = userParameters$clusterOption1, 
                  clusterOption2 = userParameters$clusterOption2, 
                  plotColors = plotColors, legendColors = legendColors, 
                  plotSymbols = plotSymbols, legendSymbols = legendSymbols, 
                  WIDTH = userParameters$img.width, HEIGHT = userParameters$img.heigth, 
                  POINTSIZE = userParameters$img.pointSize, MAXARRAY = userParameters$img.maxArray, 
                  normalization.m = userParameters$normalization.m, 
                  fileName = fileNamePrefix)
            }
            if (userParameters$raw.pca) {
                cat("\nPCA graph for raw data\n")
                groupsInLegend = !(length(unique(levels(experimentFactor)))) >= 
                  10
                gar <- pcaFun(rawData, normalized = FALSE, experimentFactor = experimentFactor, 
                  plotColors = plotColors, legendColors = legendColors, 
                  plotSymbols = plotSymbols, legendSymbols = legendSymbols, 
                  groupsInLegend = groupsInLegend, namesInPlot = ((max(nchar(sampleNames(rawData))) <= 
                    10) && (length(sampleNames(rawData)) <= (userParameters$img.maxArray/2))), 
                  WIDTH = userParameters$img.width, HEIGHT = userParameters$img.heigth, 
                  POINTSIZE = userParameters$img.pointSize, normalization.m = userParameters$normalization.m, 
                  fileName = fileNamePrefix)
            }
            if (userParameters$raw.correl) {
                cat("\nCorrelation plot for raw data\n")
                gar <- correlFun(rawData, normalized = FALSE, 
                  experimentFactor = experimentFactor, clusterOption1 = userParameters$clusterOption1, 
                  clusterOption2 = userParameters$clusterOption2, 
                  legendColors = legendColors, WIDTH = userParameters$img.width, 
                  HEIGHT = userParameters$img.heigth, POINTSIZE = userParameters$img.pointSize, 
                  MAXARRAY = userParameters$img.maxArray, normalization.m = userParameters$normalization.m, 
                  fileName = fileNamePrefix)
            }
        }
        else {
            cat("\nSkipping QC plots for the raw data.\n")
        }
        if (userParameters$normDataQC) {
            if (userParameters$loadOldNorm) {
                cat("\nLoading old normalized data\n")
                load(paste(userParameters$inputDir, userParameters$normData, 
                  sep = "/"))
            }
            if (userParameters$loadOldNorm && userParameters$perGroup) {
                cat("\nRe-ordering old normalized data per group defined in the description file.\n")
                sampleNames(normData) <- make.names(sampleNames(normData))
                matchedSamples <- match(description2[, 4], sampleNames(normData))
                if (sum(is.na(file_order2)) > 0) {
                  message <- ("Error: File names in old normalized data and file names in description file do not match!")
                  cat(message)
                  if (userParameters$createLog) 
                    sink()
                  stop(message)
                }
                normData <- normData[, file_order2]
                sampleNames(normData) <- as.character(description2[, 
                  2])
                cat("\nRe-ordering succesfull.\n")
            }
            if (userParameters$normDataQC && userParameters$statSubset) {
                cat("\nMaking subset of samples in raw data.!\n")
                x <- sampleNames(normData)[order(sampleNames(normData))]
                matchedSamples <- match(statFile[, 1], x)
                normData <- normData[, matchedSamples]
                cat("\nSuccesfully made subset of samples in raw data!\n")
            }
            cat("\nCreating QC plots for the normalized data.\n")
            fileNamePrefix <- paste(userParameters$outputDir, 
                "/", userParameters$studyName, "_NORM", sep = "")
            if (userParameters$norm.boxplot) {
                cat("\nPlot boxplot for normalized intensities\n")
                gar <- box.plot(normData, fileNamePrefix, col = plotColors, 
                  maxArray = 50)
            }
            if (userParameters$norm.density) {
                cat("\nPlot density histogram for normalized intensities\n")
                gar <- density.plot(normData, fileNamePrefix, 
                  col = plotColors, maxArray = 16)
            }
            if (userParameters$norm.cv) {
                cat("\nPlot density for coefficient of variance for normalized intensities\n")
                gar <- cv.plot(normData, fileNamePrefix, col = plotColors, 
                  maxArray = 16)
            }
            fileNamePrefix <- paste(userParameters$statisticsDir, 
                "/", userParameters$studyName, sep = "")
            if (userParameters$norm.sampleRelation) {
                cat("\nHierarchical clustering of normalized data\n")
                gar <- clusterFun(normData, normalized = TRUE, 
                  experimentFactor = experimentFactor, clusterOption1 = userParameters$clusterOption1, 
                  clusterOption2 = userParameters$clusterOption2, 
                  plotColors = plotColors, legendColors = legendColors, 
                  plotSymbols = plotSymbols, legendSymbols = legendSymbols, 
                  WIDTH = userParameters$img.width, HEIGHT = userParameters$img.heigth, 
                  POINTSIZE = userParameters$img.pointSize, MAXARRAY = userParameters$img.maxArray, 
                  normalization.m = userParameters$normalization.m, 
                  fileName = fileNamePrefix)
            }
            if (userParameters$norm.pca) {
                cat("\nPCA graph for normalized data\n")
                groupsInLegend = !(length(unique(levels(experimentFactor)))) >= 
                  10
                gar <- pcaFun(normData, normalized = TRUE, experimentFactor = experimentFactor, 
                  plotColors = plotColors, legendColors = legendColors, 
                  plotSymbols = plotSymbols, legendSymbols = legendSymbols, 
                  groupsInLegend = groupsInLegend, namesInPlot = ((max(nchar(sampleNames(normData))) <= 
                    10) && (length(sampleNames(normData)) <= 
                    (userParameters$img.maxArray/2))), WIDTH = userParameters$img.width, 
                  HEIGHT = userParameters$img.heigth, POINTSIZE = userParameters$img.pointSize, 
                  normalization.m = userParameters$normalization.m, 
                  fileName = fileNamePrefix)
            }
            if (userParameters$norm.correl) {
                cat("\nCorrelation plot for normalized data\n")
                gar <- correlFun(normData, normalized = TRUE, 
                  experimentFactor = experimentFactor, clusterOption1 = userParameters$clusterOption1, 
                  clusterOption2 = userParameters$clusterOption2, 
                  legendColors = legendColors, WIDTH = userParameters$img.width, 
                  HEIGHT = userParameters$img.heigth, POINTSIZE = userParameters$img.pointSize, 
                  MAXARRAY = userParameters$img.maxArray, normalization.m = userParameters$normalization.m, 
                  fileName = fileNamePrefix)
            }
        }
        else {
            cat("\nSkipping QC plots for the normalized data.\n")
        }
    }
    else {
        cat("\nSkipping QC plots for raw and normalized data.\n")
    }
    if (userParameters$filtering) {
        cat("\nPerforming filtering on low expressions of probes.\n")
        filtered.normData = filterFun(rawData, normData, userParameters$filter.Th, 
            userParameters$filter.dp)
        cat("\nFiltering is done!\n")
    }
    if (userParameters$createAnno) {
        cat("\nMerging annotation file with the raw expression data\n")
        anno.rawData = createAnnoFun(eset.rawData, userParameters$lib.mapping, 
            userParameters$lib.All.mapping)
        eset.anno.rawData = cbind(anno.rawData, eset.rawData)
        fileName <- paste(userParameters$outputDir, userParameters$studyName, 
            "_rawData.txt", sep = "")
        cat("\nMerging done, now saving merged file of raw data to: ", 
            fileName, "\n")
        write.table(eset.anno.rawData, file = fileName, quote = FALSE, 
            sep = "\t", row.names = F, col.names = T)
        cat("\nSaving merged raw data completed.\n")
        cat("\nMerging annotation file with the normalized expression data\n")
        anno.normData = createAnnoFun(eset.normData, userParameters$lib.mapping, 
            userParameters$lib.All.mapping)
        eset.anno.normData = cbind(anno.normData, eset.normData)
        fileName <- paste(userParameters$outputDir, userParameters$studyName, 
            "_normData_", userParameters$normType, ".txt", sep = "")
        write.table(eset.anno.normData, file = fileName, quote = FALSE, 
            sep = "\t", row.names = F, col.names = T)
        cat("\nSaving merged normalized data to the fileserver completed.\n")
    }
    else {
        cat("\nSkipping the output of the merged annotation/samples expression files.\n")
    }
    if (userParameters$createAnno && userParameters$filtering) {
        cat("\nMerging annotation file with the filtered normalized expression data\n")
        anno.filtData = createAnnoFun(filtered.normData, userParameters$lib.mapping, 
            userParameters$lib.All.mapping)
        eset.anno.filtData = cbind(anno.filtData, filtered.normData)
        fileName <- paste(userParameters$outputDir, userParameters$studyName, 
            "_normData_Filtered_", userParameters$normType, ".txt", 
            sep = "")
        write.table(eset.anno.filtData, file = fileName, quote = FALSE, 
            sep = "\t", row.names = F, col.names = T)
        cat("\nSaving merged filtered normalized data completed.\n")
    }
    else {
        cat("\nSkipping the output of the merged filtered normalized annotation/samples expression files.\n")
    }
    if (userParameters$saveHistory) {
        fileName <- paste(userParameters$outputDir, userParameters$studyName, 
            ".Rhistory", sep = "")
        savehistory(fileName)
        cat("\nSaved Rhistory in: ", fileName, "\n", sep = "")
    }
    message = ("\nIllumina Beadchip normalisation (and/or QC) has been completed!!\n")
    x <- cat(message)
    message <- paste("Completed in:", ceiling((proc.time() - 
        ptm)[3]/60), "minutes", sep = " ")
    cat(message)
    if (userParameters$createLog) 
        sink()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

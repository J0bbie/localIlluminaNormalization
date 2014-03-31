\name{getArguments}
\alias{getArguments}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gets the commandline arguments and parses these into a list using optparser.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getArguments(commandArguments)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{commandArguments}{
%%     ~~Describe \code{commandArguments} here~~
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
function (commandArguments) 
{
    inst <- installed.packages()
    notInst <- "optparse" \%in\% inst
    if (notInst) {
        library("optparse", character.only = TRUE, quietly = TRUE)
    }
    else {
        install.packages("optparse", repos = "http://R-Forge.R-project.org")
        inst <- installed.packages()
        notInst <- "optparse" \%in\% inst
        if (notInst) {
            stop("\nCould not install optparse package!")
        }
    }
    option_list <- list(make_option(c("-i", "--inputDir"), type = "character", 
        default = dirname(sys.frame(1)$ofile), help = "Path to folder where the Control_Probe_Profile, Sample_Probe_Profile and Description file are found \ndefault = [\%default] "), 
        make_option(c("-o", "--outputDir"), type = "character", 
            default = dirname(sys.frame(1)$ofile), help = "Path to folder where the output files will be stored \ndefault = [\%default] "), 
        make_option(c("-a", "--annoDir"), type = "character", 
            default = dirname(sys.frame(1)$ofile), help = "Path to folder where the annotation files for the arrays are stored. (Containing the Probes + Genes etc. on the array) \ndefault = [\%default] "), 
        make_option("--scriptDir", type = "character", default = dirname(sys.frame(1)$ofile), 
            help = "Path to folder where the scripts are stored. \ndefault = [\%default] "), 
        make_option("--species", type = "character", default = "Human", 
            help = "Species used for samples (E.g. Human, Mouse or Rat etc) \ndefault = [\%default] "), 
        make_option("--arrayType", type = "character", default = "HumanHT-12", 
            help = "Main type of array that was used (E.g. HumanHT-12 or HumanWG-6) \ndefault = [\%default] "), 
        make_option("--annoType", type = "character", default = "HumanHT-12_V4_0_R2_15002873_B", 
            help = "What annotation file should be used to check Genes+Probes etc (E.g. humanHT-12_V4_0_R2_15002873_B) \ndefault = [\%default] "), 
        make_option(c("-N", "--studyName"), type = "character", 
            default = paste(format(Sys.Date(), "\%Y-\%m-\%d"), sub(":", 
                ".", sub(":", ".", format(Sys.time(), "\%X"))), 
                sep = "_"), help = "Used to be called ns, Name of the study (Used in naming the output files) \ndefault = [\%default] "), 
        make_option("--createAnno", type = "logical", default = TRUE, 
            help = "Create an annotation file containing lumiID, probeID etc. \ndefault = [\%default] "), 
        make_option("--saveHistory", type = "logical", default = TRUE, 
            help = "Whether to save the R history to the output directory. \ndefault = [\%default] "), 
        make_option(c("-l", "--createLog"), type = "logical", 
            default = FALSE, help = "Create a log file in the output directory. Captures ALL messages from stdout and stderr. \ndefault = [\%default] "), 
        make_option(c("-B", "--statSubset"), type = "logical", 
            default = FALSE, help = "Whether statistics should be done only on a subset. (defined in --statFile) \ndefault = [\%default] "), 
        make_option(c("-f", "--loadOldNorm"), type = "logical", 
            default = FALSE, help = "Whether to load the old normalized data (given with -m/--normData). \ndefault = [\%default] "), 
        make_option(c("-s", "--sampleProbeProfilePath"), type = "character", 
            default = "Sample_Probe_Profile.txt", help = "Used to be called expFile, contains the expression values of the samples \ndefault = [\%default] "), 
        make_option(c("-c", "--controlProbeProfilePath"), type = "character", 
            default = "Control_Probe_Profile.txt", help = "Used to be called bgFile, contains the expression values of the control probes \ndefault = [\%default] "), 
        make_option(c("-d", "--descFile"), type = "character", 
            default = "descriptionFile.txt", help = "Tab-delimited file containing: arrayNames | sampleNames | sampleGroup \ndefault = [\%default] "), 
        make_option(c("-F", "--statFile"), type = "character", 
            default = "statSubsetFile.txt", help = "File containing a single column with the sampleNames or ArrayNames on which statistics should be performed. If none given, statistics is performed on all samples in the description file. \ndefault = [\%default] "), 
        make_option(c("-m", "--normData"), type = "character", 
            default = "normData.Rdata", help = "R Lumibatch object containing the normalized expression values. \ndefault = [\%default] "), 
        make_option(c("-u", "--bgSub"), type = "logical", default = FALSE, 
            help = "Whether background correction has been done already. (E.g. in Illumina GenomeStudio) \nIf FALSE -> perform bg correction using controlProbeProfileFile. If TRUE -> Skip bg correction \ndefault = [\%default] "), 
        make_option("--detectionTh", type = "double", default = 0.01, 
            help = "The p-value threshold of determining detectability of the expression. \ndefault = [\%default] "), 
        make_option("--convertNuID", type = "logical", default = TRUE, 
            help = "Determine whether to convert the probe identifier as nuID \ndefault = [\%default] "), 
        make_option("--dec", type = "character", default = ".", 
            help = "The character used in the files to indicate decimal values (E.g. '.' or ',' etc.) \ndefault = [\%default] "), 
        make_option("--parseColumnName", type = "logical", default = FALSE, 
            help = "Determine whether to parse the column names and retrieve the sample information \n(Assume the sample information is separated by a tab) \ndefault = [\%default] "), 
        make_option("--checkDupId", type = "logical", default = TRUE, 
            help = "Determine whether to check duplicated TargetIDs or ProbeIds. The duplicated ones will be averaged. \ndefault = [\%default] "), 
        make_option("--save.rawData", type = "logical", default = TRUE, 
            help = "Whether to save lumi.batch of raw data as R object in WORK.DIR \ndefault = [\%default] "), 
        make_option("--normType", type = "character", default = "lumi", 
            help = "Type of normalization to do (lumi or neqc) \ndefault = [\%default] "), 
        make_option("--bgcorrect.m", type = "character", default = "bgAdjust", 
            help = "List of the parameters for lumiExpresso{lumiB}, method for background correction. \nPossible parameters: c('none', 'bgAdjust', 'forcePositive', 'bgAdjust.affy') \ndefault = [\%default] "), 
        make_option("--variance.stabilize", type = "logical", 
            default = TRUE, help = "Whether do to variance stabilization \ndefault = [\%default] "), 
        make_option("--variance.m", type = "character", default = "log2", 
            help = "Method of variance stabilization for lumiExpresso{lumiT} package.\nPossible parameters are: c('vst', 'log2', 'cubicRoot') \ndefault = [\%default] "), 
        make_option(c("-n", "--normalize"), type = "logical", 
            default = TRUE, help = "Whether to normalize or not \ndefault = [\%default] "), 
        make_option("--normalization.m", type = "character", 
            default = "quantile", help = "List of parameters for lumiExpresso{lumiB}, method of normalization. \n parameters are: c(quantile, rsn, ssn, loess, vsn) \ndefault = [\%default] "), 
        make_option("--save.normData", type = "logical", default = TRUE, 
            help = "Whether to save lumi.batch of normalized data as R object in WORK.DIR \ndefault = [\%default] "), 
        make_option("--filtering", type = "logical", default = TRUE, 
            help = "Whether filtering of probes with a low/no expression should be performed \ndefault = [\%default] "), 
        make_option("--filter.Th", type = "double", default = 0.01, 
            help = "Threshold for probe filtering. \ndefault = [\%default] "), 
        make_option("--filter.dp", type = "integer", default = 0, 
            help = "Threshold for the count of the same probe in multiple samples \ndefault = [\%default] "), 
        make_option(c("-p", "--performStatistics"), type = "logical", 
            default = TRUE, help = "Should clustering and PCA be done alongside the normalization of the data? (Individual options are below) \ndefault = [\%default] "), 
        make_option("--rawDataQC", type = "logical", default = TRUE, 
            help = "Determine whether to do QC assessment for the raw data; if false no summary can be computed. \ndefault = [\%default] "), 
        make_option("--normDataQC", type = "logical", default = TRUE, 
            help = "Determine whether to do QC assessment for the normed data; if false no summary can be computed. \ndefault = [\%default] "), 
        make_option("--rawSummary", type = "logical", default = TRUE, 
            help = "Whether to create a summary table in WORK.DIR of the raw data \ndefault = [\%default] "), 
        make_option("--normSummary", type = "logical", default = TRUE, 
            help = "Whether to create a summary table in WORK.DIR of the normalized data \ndefault = [\%default] "), 
        make_option("--perGroup", type = "logical", default = TRUE, 
            help = "Reorder rawData lumibatch file FIRST on Group and THEN ON sampleNames (Used for the order of visualization in plots) \ndefault = [\%default] "), 
        make_option("--raw.boxplot", type = "logical", default = TRUE, 
            help = "Should a boxplot be made for the raw data. \ndefault = [\%default] "), 
        make_option("--raw.density", type = "logical", default = TRUE, 
            help = "Should a density plot be made for the raw data. \ndefault = [\%default] "), 
        make_option("--raw.cv", type = "logical", default = TRUE, 
            help = "Should a cv plot be made for the raw data. \ndefault = [\%default] "), 
        make_option("--raw.sampleRelation", type = "logical", 
            default = TRUE, help = "Should a sample relation plot be made for the raw data. \ndefault = [\%default] "), 
        make_option("--raw.pca", type = "logical", default = TRUE, 
            help = "Should a PCA be made for the raw data. \ndefault = [\%default] "), 
        make_option("--raw.correl", type = "logical", default = TRUE, 
            help = "Should a correlation plot be made for the raw data. \ndefault = [\%default] "), 
        make_option("--norm.boxplot", type = "logical", default = TRUE, 
            help = "Should a boxplot be made for the normalized data. \ndefault = [\%default] "), 
        make_option("--norm.density", type = "logical", default = TRUE, 
            help = "Should a density plot be made for the normalized data. \ndefault = [\%default] "), 
        make_option("--norm.cv", type = "logical", default = TRUE, 
            help = "Should a cv plot be made for the normalized data. \ndefault = [\%default] "), 
        make_option("--norm.sampleRelation", type = "logical", 
            default = TRUE, help = "Should a sample relation plot be made for the normalized data. \ndefault = [\%default] "), 
        make_option("--norm.pca", type = "logical", default = TRUE, 
            help = "Should a PCA be made for the normalized data. \ndefault = [\%default] "), 
        make_option("--norm.correl", type = "logical", default = TRUE, 
            help = "Should a correlation plot be made for the normalized data. \ndefault = [\%default] "), 
        make_option("--clusterOption1", type = "character", default = "Pearson", 
            help = "Distance calculation method. \nPossible parameters are c('Pearson', 'Spearman', 'Euclidean') \ndefault = [\%default] "), 
        make_option("--clusterOption2", type = "character", default = "average", 
            help = "Method of clustering. \nPossible parameters are c('Ward', 'McQuitty', 'average', 'median', 'single', 'complete', 'centroid') \ndefault = [\%default] "), 
        make_option("--img.width", type = "numeric", default = 1920, 
            help = "The max. width of the plots. \ndefault = [\%default] "), 
        make_option("--img.heigth", type = "numeric", default = 1080, 
            help = "The max. heigth of the plots. \ndefault = [\%default] "), 
        make_option("--img.pointSize", type = "numeric", default = 24, 
            help = "The size of the points on plots. \ndefault = [\%default] "), 
        make_option("--img.maxArray", type = "numeric", default = 41, 
            help = "The maximum datapoint on each plot per page. \ndefault = [\%default] "))
    userParameters <- parse_args(OptionParser(option_list = option_list), 
        commandArguments)
    userParameters <- checkUserInput(userParameters, arrayTypeList, 
        arrayAnnoList)
    return(userParameters)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

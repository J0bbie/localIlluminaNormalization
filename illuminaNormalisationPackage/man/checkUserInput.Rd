\name{checkUserInput}
\alias{checkUserInput}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
checks the validity of the user inputs, e.g. do the required files and path exist.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
checkUserInput(userParameters, arrayTypeList, arrayAnnoList)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{userParameters}{
%%     ~~Describe \code{userParameters} here~~
}
  \item{arrayTypeList}{
%%     ~~Describe \code{arrayTypeList} here~~
}
  \item{arrayAnnoList}{
%%     ~~Describe \code{arrayAnnoList} here~~
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
function (userParameters, arrayTypeList, arrayAnnoList) 
{
    userParameters$scriptDir <- correctDirectory(userParameters$scriptDir)
    userParameters$inputDir <- correctDirectory(userParameters$inputDir)
    userParameters$outputDir <- correctDirectory(userParameters$outputDir)
    userParameters$annoDir <- correctDirectory(userParameters$annoDir)
    if (userParameters$createLog) {
        fileName <- file(paste(userParameters$outputDir, userParameters$studyName, 
            "_log.txt", sep = ""))
        sink(fileName)
        sink(fileName, type = "message")
        cat("Creating log file in: ", paste(userParameters$outputDir, 
            userParameters$studyName, "_log.txt", sep = ""), 
            "\n")
    }
    if (!userParameters$loadOldNorm) {
        checkCombi <- userParameters$species \%in\% names(arrayTypeList) && 
            userParameters$arrayType \%in\% arrayTypeList[[userParameters$species]] && 
            userParameters$annoType \%in\% arrayAnnoList[[userParameters$arrayType]]
        if (!checkCombi) {
            message <- paste("\n", "The combination of species, array type and array annotation file is not correct:", 
                "\n", "- Species: ", userParameters$species, 
                "\n", "- Array type: ", userParameters$arrayType, 
                "\n", "- Annotation file: ", userParameters$annoType, 
                sep = " ")
            cat(message)
            if (userParameters$createLog) 
                sink()
            stop(message)
        }
        else {
            print("Combination of species, arrayType and annoType is OK.")
        }
        userParameters$lib.mapping = paste("lumi", userParameters$species, 
            "IDMapping", sep = "")
        userParameters$lib.All.mapping = paste("lumi", userParameters$species, 
            "All.db", sep = "")
    }
    if (file.info(userParameters$scriptDir)$isdir == FALSE) {
        message <- paste("\nThe script directory does not exist:", 
            userParameters$scriptDir, sep = "")
        cat(message)
        if (userParameters$createLog) 
            sink()
        stop(message)
    }
    if (file.info(userParameters$inputDir)$isdir == FALSE) {
        message <- paste("\nThe input directory does not exist:", 
            userParameters$inputDir, sep = "")
        cat(message)
        if (userParameters$createLog) 
            sink()
        stop(message)
    }
    if (userParameters$statSubset) {
        if (file.exists(paste(userParameters$statisticsDir, userParameters$statFile, 
            sep = "")) == FALSE) {
            message <- paste("\nNo statFile in path:", paste(userParameters$statisticsDir, 
                userParameters$statFile, sep = ""), sep = " ")
            cat(message)
            if (userParameters$createLog) 
                sink()
            stop(message)
        }
    }
    if (file.info(userParameters$outputDir)$isdir == FALSE) {
        dir.create(userParameters$outputDir)
        if (!file.info(userParameters$outputDir)$isdir) {
            if (userParameters$createLog) 
                sink()
            message <- paste("\nThe output directory does not exist and cannot be created:", 
                userParameters$outputDir, sep = "")
            cat(message)
            if (userParameters$createLog) 
                sink()
            stop(message)
        }
    }
    if (!userParameters$loadOldNorm) {
        if (file.exists(paste(userParameters$inputDir, userParameters$sampleProbeProfilePath, 
            sep = "")) == FALSE) {
            message <- paste("\nNo Sample Probe Profile in path:", 
                paste(userParameters$inputDir, userParameters$sampleProbeProfilePath, 
                  sep = ""), sep = " ")
            cat(message)
            if (userParameters$createLog) 
                sink()
            stop(message)
        }
        if (file.exists(paste(userParameters$inputDir, userParameters$controlProbeProfilePath, 
            sep = "")) == FALSE) {
            message <- paste("\nNo Control Probe Profile in path:", 
                paste(userParameters$inputDir, userParameters$controlProbeProfilePath, 
                  sep = ""), sep = " ")
            cat(message)
            if (userParameters$createLog) 
                sink()
            stop(message)
        }
        if (file.exists(paste(userParameters$outputDir, userParameters$descFile, 
            sep = "")) == FALSE) {
            message <- paste("\nNo Description file:", paste(userParameters$outputDir, 
                userParameters$descFile, sep = ""), sep = " ")
            cat(message)
            if (userParameters$createLog) 
                sink()
            stop(message)
        }
    }
    else {
        if (file.exists(paste(userParameters$inputDir, userParameters$normData, 
            sep = "")) == FALSE) {
            message <- paste("\nNo normalized R object:", paste(userParameters$inputDir, 
                userParameters$normData, sep = ""), sep = " ")
            cat(message)
            if (userParameters$createLog) 
                sink()
            stop(message)
        }
    }
    return(userParameters)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

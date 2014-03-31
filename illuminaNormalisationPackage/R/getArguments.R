#Author:                      Job van Riet
#Date of  creation:           26-2-14
#Date of modification:	26-2-14
#Version:                     1.0
#Modifications:               Original version
#Known bugs:                  None known
#Function:                    This script tries to catch all the options for the normalization used by the ArrayAnalysis scripts for Illumina normalization.
#                             Also checks if the entered parameters are all valid and begins a log file if createLog(-l) == TRUE.

#Useful short flags:
#-s       Name of the Sample_Probe_Profile
#-c       Name of the Control_Probe_Profile
#-d       Name of the descriptionFile (Tab-delimited file containing: arrayNames | sampleNames | sampleGroup)
#-N       Name of the study (Used in the naming of files)
#-l       Create a log file in the output directory (Sinks stdout and stderr)
#-i       Input directory (Where the files are found, if not the same as )
#-o       Output directory (Where the files are stored)
#-a       Directory where annotation files are kept for lumi. (Apparently not neccessary)
#-p       Whether to perform plots + PCA on data. (TRUE/FALSE)
#-n       Whether to normalize or not (TRUE/FALSE)
#-u       Whether to skip background correction (TRUE/FALSE) (TRUE means skip) 
#-m       Name of the file containing the normalized data R object   (When already normalized and wanting to do new statistics)
#-F       File containing the sampleNames/AssayNames on which statistics should be done.
#-B       Whether to perform statistics on a subset (As defined in -F)
#-f       Whether to load the old normalized data.

#Makes an optionList of the possible parameters (Using optparse library)
#Also reads the passed command-line arguments and returns these in a list
#Needs commandArgs(trailingOnly = TRUE) as input
getArguments <- function(commandArguments){
          
          #Load optparse package, if not found, try to install
          inst<-installed.packages()
          notInst <- "optparse" %in% inst
          
          #If installed, simply load
          if(notInst){
                library("optparse", character.only = TRUE, quietly=TRUE)    
          }else{
                    #Install if not installed and load then
                    install.packages("optparse", repos="http://R-Forge.R-project.org")
                    #Get a new list of all the installed packages
                    inst<-installed.packages()
                    notInst <- "optparse" %in% inst
                    if(notInst){
                              stop("\nCould not install optparse package!")
                    }
          }         
                    
          option_list <- list(
                    make_option(c("-i", "--inputDir"), type="character", default=dirname(sys.frame(1)$ofile),
                                help="Path to folder where the Control_Probe_Profile, Sample_Probe_Profile and Description file are found \ndefault = [%default] "),
                    
                    make_option(c("-o","--outputDir"), type="character", default=dirname(sys.frame(1)$ofile),
                                help = "Path to folder where the output files will be stored \ndefault = [%default] "),
                    
                    make_option(c("-a","--annoDir"), type="character", default=dirname(sys.frame(1)$ofile),
                                help="Path to folder where the annotation files for the arrays are stored. (Containing the Probes + Genes etc. on the array) \ndefault = [%default] "),
                    
                    make_option("--scriptDir", type="character", default=dirname(sys.frame(1)$ofile),
                                help="Path to folder where the scripts are stored. \ndefault = [%default] "),
                    
                    make_option("--species",  type="character", default="Human",
                                help="Species used for samples (E.g. Human, Mouse or Rat etc) \ndefault = [%default] "),
                    
                    make_option("--arrayType",  type="character", default="HumanHT-12",
                                help="Main type of array that was used (E.g. HumanHT-12 or HumanWG-6) \ndefault = [%default] "),
                    
                    make_option("--annoType",  type="character", default="HumanHT-12_V4_0_R2_15002873_B",
                                help="What annotation file should be used to check Genes+Probes etc (E.g. humanHT-12_V4_0_R2_15002873_B) \ndefault = [%default] "),
                    
                    make_option(c("-N","--studyName"),  type="character", default=paste(format(Sys.Date(), "%Y-%m-%d"), sub(":", ".", sub(":", "." ,format(Sys.time(), "%X"))), sep="_" ),
                                help="Used to be called ns, Name of the study (Used in naming the output files) \ndefault = [%default] "),
                    
                    make_option("--createAnno",  type="logical", default=TRUE,
                                help="Create an annotation file containing lumiID, probeID etc. \ndefault = [%default] "),
                    
                    make_option("--saveHistory",  type="logical", default=TRUE,
                                help="Whether to save the R history to the output directory. \ndefault = [%default] "),
                    
                    make_option(c("-l","--createLog"),  type="logical", default=FALSE,
                                help="Create a log file in the output directory. Captures ALL messages from stdout and stderr. \ndefault = [%default] "),
          
                    make_option(c("-B","--statSubset"),  type="logical", default=FALSE,
                                help="Whether statistics should be done only on a subset. (defined in --statFile) \ndefault = [%default] "),
                    
                    make_option(c("-f","--loadOldNorm"), type="logical", default=FALSE,
                                help="Whether to load the old normalized data (given with -m/--normData). \ndefault = [%default] "),
                    
                    #####################################################################################################
                    #                                       The input files paths                                       #
                    #         (E.g. /var/www/Application/studyFolder/expressionData/raw/sample_probe_profile.txt)       #
                    #####################################################################################################
                    
                    make_option(c("-s","--sampleProbeProfilePath"), type="character", default="Sample_Probe_Profile.txt",
                                help="Used to be called expFile, contains the expression values of the samples \ndefault = [%default] "),
                    
                    make_option(c("-c","--controlProbeProfilePath"), type="character", default="Control_Probe_Profile.txt",
                                help="Used to be called bgFile, contains the expression values of the control probes \ndefault = [%default] "),
                    
                    make_option(c("-d","--descFile"), type="character", default="descriptionFile.txt",
                                help="Tab-delimited file containing: arrayNames | sampleNames | sampleGroup \ndefault = [%default] "),
                    
                    make_option(c("-F","--statFile"), type="character", default="statSubsetFile.txt",
                                help="File containing a single column with the sampleNames or ArrayNames on which statistics should be performed. If none given, statistics is performed on all samples in the description file. \ndefault = [%default] "),
                    
                    make_option(c("-m","--normData"), type="character", default="normData.Rdata",
                                help="R Lumibatch object containing the normalized expression values. \ndefault = [%default] "),

                    
                    #####################################################################################################
                    #                             Variables for importing raw data with lumi                            #
                    #####################################################################################################
                    
                    make_option(c("-u", "--bgSub"), type="logical", default=FALSE,
                                help="Whether background correction has been done already. (E.g. in Illumina GenomeStudio) \nIf FALSE -> perform bg correction using controlProbeProfileFile. If TRUE -> Skip bg correction \ndefault = [%default] "),
                    
                    make_option("--detectionTh", type="double", default=0.01,
                                help="The p-value threshold of determining detectability of the expression. \ndefault = [%default] "),
                    
                    make_option("--convertNuID", type="logical", default=TRUE,
                                help="Determine whether to convert the probe identifier as nuID \ndefault = [%default] "),
                    
                    make_option("--dec", type="character", default=".",
                                help="The character used in the files to indicate decimal values (E.g. '.' or ',' etc.) \ndefault = [%default] "),          
                    
                    make_option("--parseColumnName", type="logical", default=FALSE,
                                help="Determine whether to parse the column names and retrieve the sample information \n(Assume the sample information is separated by a tab) \ndefault = [%default] "),          
                    
                    make_option("--checkDupId", type="logical", default=TRUE,
                                help="Determine whether to check duplicated TargetIDs or ProbeIds. The duplicated ones will be averaged. \ndefault = [%default] "),          
                    
                    make_option("--save.rawData", type="logical", default=TRUE,
                                help="Whether to save lumi.batch of raw data as R object in WORK.DIR \ndefault = [%default] "), 
                    
                    #####################################################################################################
                    #                             Variables lumi normalization                                          #
                    #####################################################################################################
                    
                    make_option("--normType", type="character", default="lumi",
                                help="Type of normalization to do (lumi or neqc) \ndefault = [%default] "), 
                    
                    make_option("--bgcorrect.m", type="character", default="bgAdjust",
                                help="List of the parameters for lumiExpresso{lumiB}, method for background correction. \nPossible parameters: c('none', 'bgAdjust', 'forcePositive', 'bgAdjust.affy') \ndefault = [%default] "), 
                    
                    make_option("--variance.stabilize", type="logical", default=TRUE,
                                help="Whether do to variance stabilization \ndefault = [%default] "), 
                    
                    make_option("--variance.m", type="character", default="log2",
                                help="Method of variance stabilization for lumiExpresso{lumiT} package.\nPossible parameters are: c('vst', 'log2', 'cubicRoot') \ndefault = [%default] "), 
                    
                    make_option(c("-n", "--normalize"), type="logical", default=TRUE,
                                help="Whether to normalize or not \ndefault = [%default] "), 
                    
                    make_option("--normalization.m", type="character", default="quantile",
                                help="List of parameters for lumiExpresso{lumiB}, method of normalization. \n parameters are: c(quantile, rsn, ssn, loess, vsn) \ndefault = [%default] "),
                    
                    make_option("--save.normData", type="logical", default=TRUE,
                                help="Whether to save lumi.batch of normalized data as R object in WORK.DIR \ndefault = [%default] "), 
                    
                    #####################################################################################################
                    #                                       Filtering options                                           #
                    #####################################################################################################
                    
                    make_option("--filtering", type="logical", default=TRUE,
                                help="Whether filtering of probes with a low/no expression should be performed \ndefault = [%default] "), 
                    
                    make_option("--filter.Th", type="double", default=0.01,
                                help="Threshold for probe filtering. \ndefault = [%default] "), 
                    
                    make_option("--filter.dp", type="integer", default=0,
                                help="Threshold for the count of the same probe in multiple samples \ndefault = [%default] "), 
                    
                    #####################################################################################################
                    #                                       QC/Statistics options                                       #
                    #####################################################################################################
                    
                    make_option(c("-p", "--performStatistics"), type="logical", default=TRUE,
                                help="Should clustering and PCA be done alongside the normalization of the data? (Individual options are below) \ndefault = [%default] "), 
                    
                    make_option("--rawDataQC", type="logical", default=TRUE,
                                help="Determine whether to do QC assessment for the raw data; if false no summary can be computed. \ndefault = [%default] "), 
          
                    make_option("--normDataQC", type="logical", default=TRUE,
                                help="Determine whether to do QC assessment for the normed data; if false no summary can be computed. \ndefault = [%default] "), 
                    
                    make_option("--rawSummary", type="logical", default=TRUE,
                                help="Whether to create a summary table in WORK.DIR of the raw data \ndefault = [%default] "), 
                    
                    make_option("--normSummary", type="logical", default=TRUE,
                                help="Whether to create a summary table in WORK.DIR of the normalized data \ndefault = [%default] "), 
                    
                    make_option("--perGroup", type="logical", default=TRUE,
                                help="Reorder rawData lumibatch file FIRST on Group and THEN ON sampleNames (Used for the order of visualization in plots) \ndefault = [%default] "), 
                    
                    make_option("--raw.boxplot", type="logical", default=TRUE,
                                help="Should a boxplot be made for the raw data. \ndefault = [%default] "), 
                    
                    make_option("--raw.density", type="logical", default=TRUE,
                                help="Should a density plot be made for the raw data. \ndefault = [%default] "), 
                    
                    make_option("--raw.cv", type="logical", default=TRUE,
                                help="Should a cv plot be made for the raw data. \ndefault = [%default] "), 
                    
                    make_option("--raw.sampleRelation", type="logical", default=TRUE,
                                help="Should a sample relation plot be made for the raw data. \ndefault = [%default] "), 
                    
                    make_option("--raw.pca", type="logical", default=TRUE,
                                help="Should a PCA be made for the raw data. \ndefault = [%default] "), 
                    
                    make_option("--raw.correl", type="logical", default=TRUE,
                                help="Should a correlation plot be made for the raw data. \ndefault = [%default] "), 
                    
                    make_option("--norm.boxplot", type="logical", default=TRUE,
                                help="Should a boxplot be made for the normalized data. \ndefault = [%default] "), 
                    
                    make_option("--norm.density", type="logical", default=TRUE,
                                help="Should a density plot be made for the normalized data. \ndefault = [%default] "), 
                    
                    make_option("--norm.cv", type="logical", default=TRUE,
                                help="Should a cv plot be made for the normalized data. \ndefault = [%default] "), 
                    
                    make_option("--norm.sampleRelation", type="logical", default=TRUE,
                                help="Should a sample relation plot be made for the normalized data. \ndefault = [%default] "), 
                    
                    make_option("--norm.pca", type="logical", default=TRUE,
                                help="Should a PCA be made for the normalized data. \ndefault = [%default] "), 
                    
                    make_option("--norm.correl", type="logical", default=TRUE,
                                help="Should a correlation plot be made for the normalized data. \ndefault = [%default] "), 
                    
                    make_option("--clusterOption1", type="character", default="Pearson",
                                help="Distance calculation method. \nPossible parameters are c('Pearson', 'Spearman', 'Euclidean') \ndefault = [%default] "), 
                    
                    make_option("--clusterOption2", type="character", default="average",
                                help="Method of clustering. \nPossible parameters are c('Ward', 'McQuitty', 'average', 'median', 'single', 'complete', 'centroid') \ndefault = [%default] "),
                    
                    #####################################################################################################
                    #                             Display parameters for the images                                     #
                    #####################################################################################################
                    
                    make_option("--img.width", type="numeric", default=1920,
                                help="The max. width of the plots. \ndefault = [%default] "), 
                    
                    make_option("--img.heigth", type="numeric", default=1080,
                                help="The max. heigth of the plots. \ndefault = [%default] "), 
                    
                    make_option("--img.pointSize", type="numeric", default=24,
                                help="The size of the points on plots. \ndefault = [%default] "), 
                    
                    make_option("--img.maxArray", type="numeric", default=41,
                                help="The maximum datapoint on each plot per page. \ndefault = [%default] ")
          )
          
          
          #Get a list of the named parameters that were given when running this script
          userParameters <- parse_args(OptionParser(option_list = option_list), commandArguments)
          
          
          #Check the validity of the parameters and also create the log file if needed
          userParameters <- checkUserInput(userParameters, arrayTypeList, arrayAnnoList)
          
          #Return the parameters if all were valid
          return(userParameters)
}
#Check if the required parameters are all valid. Create a log file is createLog == TRUE
#Also check if the given species, array type and array annotation combi is valid.
checkUserInput <-function(userParameters, arrayTypeList, arrayAnnoList) {
          #Check if the directories exist, also clean their path if not properly closed of with last /
          userParameters$scriptDir      <- correctDirectory(userParameters$scriptDir)
          userParameters$inputDir       <- correctDirectory(userParameters$inputDir)
          userParameters$outputDir      <- correctDirectory(userParameters$outputDir)
          userParameters$annoDir        <- correctDirectory(userParameters$annoDir)
          
          #Create a logFile in the outputdirectory
          if(userParameters$createLog){
                    fileName <- file(paste(userParameters$outputDir, userParameters$studyName, "_log.txt", sep = ""))
                    sink(fileName)
                    sink(fileName, type="message")
                    cat("Creating log file in: ", paste(userParameters$outputDir, userParameters$studyName, "_log.txt", sep = "") ,"\n")
          }          
          
          #Dont check arguments if only the statistics is being done
          if(!userParameters$loadOldNorm){
                    #Check if combination of species, arrayType and arrayAnnotation is valid.
                    checkCombi <- userParameters$species %in% names(arrayTypeList) && userParameters$arrayType %in% arrayTypeList[[userParameters$species]] && userParameters$annoType %in% arrayAnnoList[[userParameters$arrayType]]
                    if (!checkCombi) {
                              message <- paste('\n' , "The combination of species, array type and array annotation file is not correct:", '\n' ,
                                               "- Species: ", userParameters$species, '\n' ,
                                               "- Array type: ", userParameters$arrayType, '\n' ,
                                               "- Annotation file: ", userParameters$annoType, sep=" ")
                              cat(message)
                              if(userParameters$createLog) sink()
                              stop (message)
                              
                              
                    } 
                    else {
                              print("Combination of species, arrayType and annoType is OK.") 
                    }       
                    
                    #Add correct library for mapping based on species
                    userParameters$lib.mapping = paste( "lumi", userParameters$species, "IDMapping", sep="");
                    userParameters$lib.All.mapping = paste( "lumi", userParameters$species, "All.db", sep="");
          }
          if (file.info(userParameters$scriptDir)$isdir == FALSE){
                    message <- paste("\nThe script directory does not exist:",userParameters$scriptDir, sep="")
                    cat(message)
                    if(userParameters$createLog) sink()
                    stop(message)
          }
          
          if (file.info(userParameters$inputDir)$isdir == FALSE){
                    message <- paste("\nThe input directory does not exist:",userParameters$inputDir, sep="")
                    cat(message)
                    if(userParameters$createLog) sink()
                    stop(message)
          }
          
          #If statistics should only be performed on a smaller subset of samples, check if the file containing the sampleNames of this subset exist.
          if(userParameters$statSubset){
                    if (file.exists(paste(userParameters$statisticsDir, userParameters$statFile, sep="")) == FALSE){
                              message <- paste("\nNo statFile in path:", paste(userParameters$statisticsDir, userParameters$statFile, sep=""), sep=" ")
                              cat(message)
                              if(userParameters$createLog) sink()
                              stop(message)
                    }        
          }
          
          #Make output directory if not yet exist
          if (file.info(userParameters$outputDir)$isdir == FALSE){
                    dir.create(userParameters$outputDir)
                    if(!file.info(userParameters$outputDir)$isdir){
                              if(userParameters$createLog) sink()
                              message <- paste("\nThe output directory does not exist and cannot be created:",userParameters$outputDir, sep="")
                              cat(message)
                              if(userParameters$createLog) sink()
                              stop(message)
                    }
          }
          if(!userParameters$loadOldNorm){
                    #Check if the paths to the input files are all valid
                    if (file.exists(paste(userParameters$inputDir, userParameters$sampleProbeProfilePath, sep="")) == FALSE){
                              message <- paste("\nNo Sample Probe Profile in path:", paste(userParameters$inputDir, userParameters$sampleProbeProfilePath, sep=""), sep=" ")
                              cat(message)
                              if(userParameters$createLog) sink()
                              stop(message)
                    }
                    if (file.exists(paste(userParameters$inputDir, userParameters$controlProbeProfilePath, sep="")) == FALSE){
                              message <- paste("\nNo Control Probe Profile in path:", paste(userParameters$inputDir, userParameters$controlProbeProfilePath, sep=""), sep=" ")
                              cat(message)
                              if(userParameters$createLog) sink()
                              stop(message)
                    }
                    if (file.exists(paste(userParameters$outputDir, userParameters$descFile, sep="")) == FALSE){
                              message <- paste("\nNo Description file:", paste(userParameters$outputDir, userParameters$descFile, sep="") , sep=" ")
                              cat(message)
                              if(userParameters$createLog) sink()
                              stop(message)
                    }
          }else{
                    if (file.exists(paste(userParameters$inputDir, userParameters$normData, sep="")) == FALSE){
                              message <- paste("\nNo normalized R object:", paste(userParameters$inputDir, userParameters$normData, sep="") , sep=" ")
                              cat(message)
                              if(userParameters$createLog) sink()
                              stop(message)
                    }         
          }
          return(userParameters)
}

#Amend paths of directories if not started or closed off correctly with /
correctDirectory <- function(dirName) {
          lastChar <- substr(dirName,nchar(dirName)-1,nchar(dirName))
          if(lastChar != "/"){
                    dirName <- paste(dirName,"/",sep="")
          }
          return(dirName)
}

#####################################################################################################
#                             Lists of possible arrays and annotations                              #
#####################################################################################################

arrayTypeList = list(
          Human = c( "HumanHT-12", "HumanRef-8", "HumanWG-6"),
          Mouse = c( "MouseRef-8", "MouseWG-6"),
          Rat   = c( "RatRef-8")
)

arrayAnnoList = list(
          `HumanHT-12` = c("HumanHT-12_V4_0_R2_15002873_B_WGDASL", 
                           "HumanHT-12_V4_0_R2_15002873_B", 
                           "HumanHT-12_V4_0_R1_15002873_B", 
                           "HumanHT-12_V3_0_R2_11283641_A",
                           "HumanHT-12_V3_0_R3_11283641_A"),
          `HumanRef-8` = c("HumanRef-8_V3_0_R3_11282963_A", 
                           "HumanRef-8_V3_0_R2_11282963_A", 
                           "HUMANREF-8_V3_0_R1_11282963_A_WGDASL", 
                           "HumanRef-8_V2_0_R4_11223162_A"),
          `HumanWG-6`  = c("HumanWG-6_V2_0_R4_11223189_A", 
                           "HumanWG-6_V3_0_R2_11282955_A",
                           "HumanWG-6_V3_0_R3_11282955_A"),
          `MouseRef-8` = c("MouseRef-8_V1_1_R4_11234312_A", 
                           "MouseRef-8_V2_0_R2_11278551_A",
                           "MouseRef-8_V2_0_R3_11278551_A"),
          `MouseWG-6` = c( "MouseWG-6_V1_1_R4_11234304_A", 
                           "MouseWG-6_V2_0_R2_11278593_A",
                           "MouseWG-6_V2_0_R3_11278593_A"),
          `RatRef-12` = c( "RatRef-12_V1_0_R5_11222119_A")
)
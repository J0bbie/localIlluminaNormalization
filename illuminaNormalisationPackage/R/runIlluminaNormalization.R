#Author:                      Job van Riet + ArrayAnalysis.org
#Date of  creation:           26-2-14
#Date of modification:        26-2-14
#Version:                     1.0
#Modifications:               Original version
#Known bugs:                  None known
#Function:                    This script functions as the main script, calling the functions of the other scripts.
#                             This pipeline is used to normalize Illumina microarray data. (From ArrayAnalysis.org)
#                             This script is meant to be run by openCPU using the TNO dashboard.

#Main function to run the normalisation, needs the command line arguments (Curl call)
runIlluminaNormalisation <- function (commandArgsCurl){
          
          #####################################################################################################
          #                                       Main Process flow                                          #
          #####################################################################################################
          
          #####################################################################################################
          #                             Load additional scripts, packages and parameters                      #
          #####################################################################################################
          
          #Keep track of the running time of this script.
          ptm <- proc.time()
          
          #Path to folder where the R scripts are found for the normalization of Illumina expression data. (Ideally, where this script is located also)
          #Can also be given as parameters, but this is to load the getArguments script
          SCRIPT.DIR <- dirname(sys.frame(1)$ofile)
          
          #Functions to get the passed parameters and set the default values of all parameters used in this pipeline
          source(paste(SCRIPT.DIR,"getArguments.R",sep="/"))
          
          #Get the command-line parameters that were given to this script (Parameters defined in getArguments.R)
          #Also check the validity of these parameters and directories
          userParameters <- getArguments(commandArgsCurl)
          
          #Sample for hardcoding arguments:
          #userParameters <- getArguments("-h")
          #or
          #userParameters <- getArguments(c("-o", "C:/Users/rietjv/AppData/Local/My Local Documents/output2", "-i" ,"//tsn.tno.nl/Data/Users/rietjv/Home/Documents/Database/extra/Report","-s","Sample_Probe_Profile_102259-2.txt","-c", "Control_Probe_Profile_102259-2.txt","-d","headers.txt"))
          
          #Function to install missing libraries
          source(paste(userParameters$scriptDir,"functions_loadPackages.R",sep="/"))
          
          #Functions for the creation of the plots
          source(paste(userParameters$scriptDir,"functions_makeImages.R",sep="/"))
          
          #Functions for QCing the raw and normalized data
          source(paste(userParameters$scriptDir,"functions_qualityControl.R",sep="/"))
          
          cat("\nLoading required packages.\n")
          
          #Create a list of the mandatory packages needed for this pipeline.
          pkgs <- c( "limma", "ALL","bioDist", "gplots",
                     "annotate", "arrayQualityMetrics",
                     switch(userParameters$species,
                            Human = pkgs <- c("org.Hs.eg.db"),
                            Mouse = pkgs <- c("org.Mm.eg.db"),
                            Rat =   pkgs <- c("org.Rn.eg.db") 
                     ), "lumi",
                     userParameters$lib.mapping,
                     userParameters$lib.All.mapping
          )
          
          #Install any missing R libraries if needed
          loadPackages(pkgs)
          
          cat("\nRequired packages succesfully loaded.\n")
          
          ##################################################################################
          ##        Load description file (arrayNames | sampleNames | sampleGroup)        ##
          ##################################################################################
          
          #If normalization is true:
          if(userParameters$normalize){
                    ###############################################################################
                    ## Load RAW data, perform pre-processing and normalization                   ##
                    ###############################################################################
                    
                    expData <- paste(userParameters$inputDir, userParameters$sampleProbeProfilePath, sep="")
                    
                    cat("\nLoading sample probe profile:", expData,"\n", sep="")
                    
                    rawData <- import.rawData(expData, userParameters$detectionTh ,userParameters$convertNuID, 
                                              userParameters$checkDupId, userParameters$lib.mapping, userParameters$dec, 
                                              userParameters$parseColumnName, userParameters$rawDataQC);
                    
                    #create generic sampleNames with function make.names
                    sampleNames(rawData)<- make.names(sampleNames(rawData))
                    
                    cat("\nSuccesfully loaded the Sample Probe Profile.\n")
          
                    
                    ##################################################################################
                    ##                            Check description file                            ##
                    ##################################################################################
                    
                    cat("\nChecking if description data is valid for the given sample probe profile.\n")
                    
                    #Create new column with format sampleNames as read-in with make.names
                    description$arraySampleNames = make.names(description[,1])
                    
                    #Check if all the arrays in the Sample_Probe_File have been named in the descriptionFile
                    if( sum( length(sampleNames(rawData)) - length(description[,1]))  > 0){
                              message <- paste("Error: Number of array names in raw data file and number of array names in description file is not of the same size!")
                              cat(message)
                              if(userParameters$createLog) sink()
                              stop(message)
                    }
                    
                    #Match sampleNames from datafile with first column from description file
                    file_order <- match(description[,4],sampleNames(rawData))
                    
                    #Check on NA values in file_order; if na in file_order stop
                    if(sum(is.na(file_order)) > 0){
                              message <- paste("\nError: Assigned array names in raw data file and file names in description file do not match!\n")
                              cat(message)
                              if(userParameters$createLog) sink()
                              stop(message)
                    }
                    
                    #Check if every sampleName is unique in description file
                    if(length(description[,2]) != length(unique(description[,2])) ){
                              message <- ("Error: Assigned sampleNames are not unique!")
                              cat(message)
                              if(userParameters$createLog) sink()
                              stop(message)
                    }
                    
                    #Change order of rawData in order of file_order
                    rawData <- rawData[,file_order]
                    
                    cat("\nDescription data is valid.\n")
                    
                    ##################################################################################
                    ##        Reorder rawData lumibatch file on Group and sampleNames               ##
                    ##################################################################################
                    
                    #Reorder the samples per defined group, this makes sure the samples in a group are shown together in the plots.
                    if(userParameters$perGroup){
                              cat("\nRe-ordering raw Sample Probe Profile per group defined in the description file.\n")
                              
                              #Match sampleNames from datafile with first column from description file
                              file_order2 <- match(description2[,4],sampleNames(rawData))
                              
                              #If not all the array have a sample name
                              if(sum(is.na(file_order2)) > 0) {
                                        message <- ("Error: File names in Sample Probe Profile and file names in description file do not match!")
                                        cat(message)
                                        if(userParameters$createLog) sink()
                                        stop(message)
                              }
                              #Reorder the raw expression data
                              rawData <- rawData[,file_order2]
                              
                              #Change sampleNames into reordered description file
                              sampleNames(rawData)<- as.character(description2[,2]) 
                              
                              cat("\nRe-ordering succesfull.\n")
                    }else {
                              #Change sampleNames into loaded description file
                              sampleNames(rawData)<- as.character(description[,2])
                              cat("\nSample names have been given to the arrays.\n")
                    }
                    
                    ##################################################################################
                    ##                            Background correction                             ##
                    ##################################################################################
                    
                    #If background corrections has already been performed
                    if(userParameters$bgSub) { 
                              cat("\nSkipping background correction\n", "\nNormalizing the raw Sample Probe Profiles:", userParameters$sampleProbeProfilePath, "\n", sep="")
                              
                              #Normalize lumiBatch object of raw data
                              normData <- lumi.normData(rawData, 
                                                        bg.correct=FALSE, userParameters$bgcorrect.m,
                                                        userParameters$variance.stabilize, userParameters$variance.m,
                                                        userParameters$normalize, userParameters$normalization.m, 
                                                        userParameters$normDataQC);
                              
                              cat("\nNormalization of raw data has been successfull.\n")   
                              
                    }else{    
                              #Background correction has not been done already
                              controlData <- paste(userParameters$inputDir, userParameters$controlProbeProfilePath, sep="")
                              cat("\nPerforming background correction (", userParameters$bgcorrect.m, ") on the Sample Probe Profile using the Control Probe Profile: ",controlData, "\n", sep="")
                              cat("\nLoading Control Probe Profile:", userParameters$controlProbeProfilePath, "\n", sep="")
                              
                              #Checks headers of controlData file with rawData object
                              #Add control data to the rawData lumiBatch object file
                              cat("\nCombining Control data with Sample data.\n")
                              rawData.ctrl <- addControlData2lumi(controlData, rawData)
                              
                              #Get control data in a data.frame
                              controlData <- as.data.frame(getControlData(rawData.ctrl), row.names = NULL )
                              
                              #Normalize lumi batch raw data object using 'lumi' or 'neqc' function
                              cat("\nNormaling (", userParameters$normType ,") the raw Sample Probe Profile using background correction:", userParameters$controlProbeProfilePath, "\n", sep="")
                              
                              switch (userParameters$normType,
                                      lumi = normData <- lumi.normData(rawData.ctrl, 
                                                                       bg.correct=TRUE, userParameters$bgcorrect.m,
                                                                       userParameters$variance.stabilize, userParameters$variance.m,
                                                                       userParameters$normalize, userParameters$normalization.m, 
                                                                       userParameters$normDataQC),
                                      neqc = normData <- neqc.normData(rawData.ctrl, controlData)
                              )
                              
                              cat("\nNormalization of raw data with background correction has been successfull.\n")  
                    }
                    
                    #Create rawData exprs eset table
                    cat("\nCreating eSets expression matrix of raw data.\n") 
                    eset.rawData <- exprs(rawData)
                    cat("\nSuccesfully created eSets expression matrix of raw data.\n") 
                    
                    #Create normData exprs eset table
                    cat("\nCreating eSets expression matrix of normalized data.\n") 
                    eset.normData <- exprs(normData)
                    cat("\nSuccesfully created eSets expression matrix of normalized data.\n") 
          
                    ##################################################################################
                    ##                            Summary tables                                    ##
                    ##################################################################################
                    
                    #Make a summary of the raw data
                    if (userParameters$rawSummary){
                              fileName <- paste(userParameters$outputDir, userParameters$studyName,"_summary_rawData.txt", sep="") 
                              cat("\nCreating summary file of the means and SD of the raw data: ", fileName ,"\n", sep="")
                              rawSum.table = createSummary(rawData, fileName);
                              cat("\nSuccesfully made summary file of raw data.\n")          
                    }
                    
                    #Make a summary of the normalized data
                    if(userParameters$normSummary){
                              fileName <- paste(userParameters$outputDir, userParameters$studyName,"_summary_normData.txt", sep="") 
                              cat("\nCreating summary file of the means and SD of the normalized data: ", fileName ,"\n", sep="")
                              normSum.table = createSummary(normData, fileName );
                              cat("\nSuccesfully made summary file of normalized data.\n") 
                    }                   
          
                    ##################################################################################
                    ##                  Save lumiBatch files as a R object                          ##
                    ##################################################################################
                    
                    #Save lumiBatch R object of rawData
                    if(userParameters$save.rawData) {
                              fileName <- paste(userParameters$outputDir, userParameters$studyName,"_rawData.Rdata", sep="")
                              cat("\nSaving lumiBatch R object of the raw data in: ", fileName ,"\n", sep="")
                              save(rawData, file=fileName )
                              cat("\nSuccesfully saved lumiBatch R object of the raw data\n")
                    }
                    
                    #Save lumiBatch R object of normData
                    if(userParameters$save.normData) {
                              fileName <- paste(userParameters$outputDir, userParameters$studyName,"_normData.Rdata", sep="")
                              cat("\nSaving lumiBatch R object of the normalized data in: ", fileName ,"\n", sep="")
                              save(normData, file=fileName )
                              cat("\nSuccesfully saved lumiBatch R object of the normalized data\n")
                    }
          } #End background correction/normalization
          
          ##################################################################################
          ##                  Make QC plots of raw and/or normalized data                 ##
          ##################################################################################
          
          if(userParameters$performStatistics){
                    
                    #################################################################################
                    #                 Open subset file if user selected a subset                     #
                    #################################################################################
                    
                    if(userParameters$statSubset){
                              cat(paste("\nReading in the statFile conating the subset of samples on which to perform the statistics: ,",paste(userParameters$inputDir, userParameters$statFile, sep=""), "\n", sep=""))
                              statFile <- read.table(paste(userParameters$inputDir, userParameters$statFile, sep="/"),
                                                     header=F,  
                                                     stringsAsFactors = F,
                                                     sep='\t',
                                                     quote="")
                              
                              #Get the indexes of the matched samples to the samples in the normData object
                              #All other samples will not be used in the statistics
                              statFile$arraySampleNames = make.names(statFile[,1])
                              
                              #Get the groups from the description file and add to the statFile dataFrame
                              statFile$groups <- description$FactorValue[match(statFile[,1], description$SourceName)]
                              
                              cat("\nSuccesfully read in the statFile containing the samples for the subset!\n")
                              
                              if(userParameters$perGroup){
                                        #Order the groups in ascending order
                                        statFile =  statFile[order(statFile$groups, statFile$arraySampleNames), ]
                              }
                    }
                                        
                    ###############################################################################
                    # Create array groups, array names and plot variables                         #
                    ###############################################################################
                    
                    cat("\nCreating a plot colorset for each array group.\n")
                    
                    if(userParameters$statSubset == FALSE){
                              #Create colorset for the array groups
                              if(userParameters$perGroup){
                                        #Use reordered  description file ordered per group
                                        experimentFactor <- as.factor(description2[,3])
                                        colList          <- colorsByFactor(experimentFactor)
                                        plotColors       <- colList$plotColors
                                        legendColors     <- colList$legendColors
                                        rm(colList)
                              } else {
                                        #Use originaly loaded description file
                                        experimentFactor <- as.factor(description[,3])
                                        colList          <- colorsByFactor(experimentFactor)
                                        plotColors       <- colList$plotColors
                                        legendColors     <- colList$legendColors
                                        rm(colList)
                              }
                    }else{
                              #Use the samples and groups from the subset
                              experimentFactor <- as.factor(statFile$groups)
                              colList          <- colorsByFactor(experimentFactor)
                              plotColors       <- colList$plotColors
                              legendColors     <- colList$legendColors
                              rm(colList)
                    }
                    
                    #Create symbolset for the array groups
                    plotSymbols <- 18-as.numeric(experimentFactor)
                    legendSymbols <- sort(unique(plotSymbols), decreasing=TRUE)
                    
                    cat("\nPlot colorset sucesfully made.\n")
          
                    #################################################################################
                    #                             QC plots of the raw data                          #
                    #################################################################################
                    
                    if(userParameters$rawDataQC){
                              
                              #Make subset of raw data if needed
                              if(userParameters$statSubset){
                                        cat("\nMaking subset of samples in raw data.!\n")
                                        
                                        x <- sampleNames(rawData)[order(sampleNames(rawData))]
                                        matchedSamples <- match(statFile[,1], x)
                                        
                                        #Make subset of samples
                                        rawData <- rawData[, matchedSamples]
                                        
                                        cat("\nSuccesfully made subset of samples in raw data!\n")   
                              }  
                              
                              cat("\nCreating QC plots for the raw data.\n")
                              fileNamePrefix <- paste(userParameters$outputDir, "/",  userParameters$studyName , "_RAW" ,sep="")
                              if(userParameters$raw.boxplot) {
                                        cat("\nPlot boxplot for raw intensities\n")
                                        gar <-box.plot(rawData, fileNamePrefix , col=plotColors, maxArray=50)
                              }
                              
                              if(userParameters$raw.density) {
                                        cat("\nPlot density histogram for raw intensities\n")
                                        gar <- density.plot(rawData, fileNamePrefix , col=plotColors, maxArray=16)
                              }
                              
                              if(userParameters$raw.cv) {
                                        cat("\nPlot density for coefficient of variance for raw intensities\n")
                                        gar <- cv.plot(rawData, fileNamePrefix, col=plotColors, maxArray=16)                    }
                              
                              fileNamePrefix <- paste(userParameters$outputDir, "/", userParameters$studyName,sep="")
                              if(userParameters$raw.sampleRelation) {
                                        cat("\nHierarchical clustering of raw data\n") 
                                        gar <- clusterFun(rawData, normalized=FALSE,  experimentFactor=experimentFactor,
                                                          clusterOption1=userParameters$clusterOption1, clusterOption2=userParameters$clusterOption2,
                                                          plotColors=plotColors, legendColors=legendColors,
                                                          plotSymbols=plotSymbols, legendSymbols=legendSymbols,
                                                          WIDTH=userParameters$img.width,HEIGHT=userParameters$img.heigth,
                                                          POINTSIZE=userParameters$img.pointSize,MAXARRAY=userParameters$img.maxArray,
                                                          normalization.m=userParameters$normalization.m, fileName=fileNamePrefix) 
                              }
                              
                              if(userParameters$raw.pca) {  
                                        cat("\nPCA graph for raw data\n")
                                        groupsInLegend =  !( length(unique(levels(experimentFactor))) ) >=10
                                        
                                        gar <- pcaFun(rawData, normalized=FALSE ,experimentFactor=experimentFactor, 
                                                      plotColors=plotColors, legendColors=legendColors, plotSymbols=plotSymbols,
                                                      legendSymbols=legendSymbols, groupsInLegend=groupsInLegend,
                                                      namesInPlot=((max(nchar(sampleNames(rawData)))<=10)&& (length(sampleNames(rawData))<=(userParameters$img.maxArray/2))),
                                                      WIDTH=userParameters$img.width,HEIGHT=userParameters$img.heigth,
                                                      POINTSIZE=userParameters$img.pointSize, normalization.m=userParameters$normalization.m, 
                                                      fileName=fileNamePrefix)
                              }
                              
                              if(userParameters$raw.correl){  
                                        cat("\nCorrelation plot for raw data\n")
                                        gar <- correlFun(rawData, normalized=FALSE, experimentFactor=experimentFactor, 
                                                         clusterOption1=userParameters$clusterOption1, clusterOption2=userParameters$clusterOption2,
                                                         legendColors=legendColors,
                                                         WIDTH=userParameters$img.width,HEIGHT=userParameters$img.heigth,
                                                         POINTSIZE=userParameters$img.pointSize,MAXARRAY=userParameters$img.maxArray,
                                                         normalization.m=userParameters$normalization.m, fileName=fileNamePrefix)   
                              }
                    }else{
                              cat("\nSkipping QC plots for the raw data.\n")   
                    }
                    
                    if(userParameters$normDataQC){
                              
                              #################################################################################
                              #                   Perform statistics on old normalized data                   #
                              #################################################################################
                              #Load the "old" normalized data in if there should be new statistics performed
                              if(userParameters$loadOldNorm){
                                        cat("\nLoading old normalized data\n")
                                        load(paste(userParameters$inputDir, userParameters$normData, sep="/"))
                              }
                              
                              #Reorder the samples per defined group, this makes sure the samples in a group are shown together in the plots.
                              if(userParameters$loadOldNorm && userParameters$perGroup){
                                        cat("\nRe-ordering old normalized data per group defined in the description file.\n")
                                        
                                        #create generic sampleNames with function make.names
                                        sampleNames(normData) <- make.names(sampleNames(normData))
                                        
                                        #Match sampleNames from datafile with first column from description file
                                        matchedSamples<- match(description2[,4],sampleNames(normData))
                                        
                                        #If not all the array have a sample name
                                        if(sum(is.na(file_order2)) > 0) {
                                                  message <- ("Error: File names in old normalized data and file names in description file do not match!")
                                                  cat(message)
                                                  if(userParameters$createLog) sink()
                                                  stop(message)
                                        }
                                        #Reorder the normed expression data
                                        normData <- normData[,file_order2]
                                        
                                        #Change sampleNames into reordered description file
                                        sampleNames(normData) <- as.character(description2[,2]) 
                                        
                                        cat("\nRe-ordering succesfull.\n")
                              }
                              
                              #Make subset of normed data if needed
                              if(userParameters$normDataQC && userParameters$statSubset){
                                        cat("\nMaking subset of samples in raw data.!\n")
                                        
                                        x <- sampleNames(normData)[order(sampleNames(normData))]
                                        matchedSamples <- match(statFile[,1], x)
                                        
                                        #Make subset of samples
                                        normData <- normData[, matchedSamples]
                                        
                                        cat("\nSuccesfully made subset of samples in raw data!\n")   
                              }
                              
                              #################################################################################
                              #                             QC plots of the normalized data                   #
                              #################################################################################
                              
                              cat("\nCreating QC plots for the normalized data.\n")
                              fileNamePrefix <- paste(userParameters$outputDir, "/",  userParameters$studyName , "_NORM" ,sep="")
                              if(userParameters$norm.boxplot) {
                                        cat("\nPlot boxplot for normalized intensities\n")
                                        gar <- box.plot(normData, fileNamePrefix , col=plotColors, maxArray=50)
                              }
                              
                              if(userParameters$norm.density) {
                                        cat("\nPlot density histogram for normalized intensities\n")
                                        gar <- density.plot(normData, fileNamePrefix , col=plotColors, maxArray=16)
                              }
                              
                              if(userParameters$norm.cv) {
                                        cat("\nPlot density for coefficient of variance for normalized intensities\n")
                                        gar <- cv.plot(normData, fileNamePrefix , col=plotColors, maxArray=16)
                              }
                              
                              fileNamePrefix <- paste(userParameters$statisticsDir, "/", userParameters$studyName,sep="")
                              if(userParameters$norm.sampleRelation) {
                                        cat("\nHierarchical clustering of normalized data\n") 
                                        gar <- clusterFun(normData, normalized=TRUE,  experimentFactor=experimentFactor,
                                                          clusterOption1=userParameters$clusterOption1, clusterOption2=userParameters$clusterOption2,
                                                          plotColors=plotColors, legendColors=legendColors,
                                                          plotSymbols=plotSymbols, legendSymbols=legendSymbols,
                                                          WIDTH=userParameters$img.width,HEIGHT=userParameters$img.heigth,
                                                          POINTSIZE=userParameters$img.pointSize,MAXARRAY=userParameters$img.maxArray,
                                                          normalization.m=userParameters$normalization.m, fileName=fileNamePrefix)
                                        
                              }
                              
                              if(userParameters$norm.pca) {  
                                        cat("\nPCA graph for normalized data\n")
                                        groupsInLegend =  !( length(unique(levels(experimentFactor))) ) >=10
                                        
                                        gar <- pcaFun(normData, normalized=TRUE ,experimentFactor=experimentFactor, 
                                                      plotColors=plotColors, legendColors=legendColors, plotSymbols=plotSymbols,
                                                      legendSymbols=legendSymbols, groupsInLegend=groupsInLegend,
                                                      namesInPlot=((max(nchar(sampleNames(normData)))<=10)&& (length(sampleNames(normData))<=(userParameters$img.maxArray/2))),
                                                      WIDTH=userParameters$img.width,HEIGHT=userParameters$img.heigth,
                                                      POINTSIZE=userParameters$img.pointSize, normalization.m=userParameters$normalization.m, 
                                                      fileName=fileNamePrefix)
                                        
                              }
                              
                              if(userParameters$norm.correl){  
                                        cat("\nCorrelation plot for normalized data\n")
                                        gar <- correlFun(normData, normalized=TRUE, experimentFactor=experimentFactor, 
                                                         clusterOption1=userParameters$clusterOption1, clusterOption2=userParameters$clusterOption2,
                                                         legendColors=legendColors,
                                                         WIDTH=userParameters$img.width,HEIGHT=userParameters$img.heigth,
                                                         POINTSIZE=userParameters$img.pointSize,MAXARRAY=userParameters$img.maxArray,
                                                         normalization.m=userParameters$normalization.m, fileName=fileNamePrefix)
                              }
                    }else{
                              cat("\nSkipping QC plots for the normalized data.\n")
                    }
          }else{
                    cat("\nSkipping QC plots for raw and normalized data.\n")
          }
          
          
          #################################################################################
          #                             Perform filtering                                 #
          #################################################################################
          
          #Will create a matrix a filter on expression lower than the user defined expression threshold and bead count
          #This table will be later be used to create a filtered annotation matrix
          if(userParameters$filtering){
                    cat("\nPerforming filtering on low expressions of probes.\n")
                    filtered.normData = filterFun(rawData, normData, userParameters$filter.Th, userParameters$filter.dp)
                    cat("\nFiltering is done!\n")
          }
          
          
          #################################################################################
          #                             Create annotation files                           #
          #################################################################################
          
          #Make a annotationfile for the raw data. (Containing gene-expressions per gene for each sample)
          if(userParameters$createAnno){
                    
                    cat("\nMerging annotation file with the raw expression data\n")
                    
                    #Merge samples to genes
                    anno.rawData = createAnnoFun(eset.rawData, userParameters$lib.mapping, userParameters$lib.All.mapping);
                    eset.anno.rawData  = cbind(anno.rawData, eset.rawData);
                    
                    fileName <- paste(userParameters$outputDir,   userParameters$studyName,"_rawData.txt", sep="")
                    
                    cat("\nMerging done, now saving merged file of raw data to: ", fileName ,"\n")
                    
                    write.table(eset.anno.rawData, file = fileName, quote= FALSE, sep='\t', row.names= F, col.names= T )
                    
                    cat("\nSaving merged raw data completed.\n")
                    
                    cat("\nMerging annotation file with the normalized expression data\n")
                    
                    #Merge norm data with anno
                    anno.normData = createAnnoFun(eset.normData, userParameters$lib.mapping, userParameters$lib.All.mapping);
                    eset.anno.normData = cbind(anno.normData, eset.normData);
                    
                    fileName <- paste(userParameters$outputDir,  userParameters$studyName, "_normData_", userParameters$normType, ".txt", sep="")
                    
                    write.table(eset.anno.normData, file = fileName, quote= FALSE, sep='\t', row.names= F, col.names= T )
                    
                    cat("\nSaving merged normalized data to the fileserver completed.\n")
                    
          }else{
                    cat("\nSkipping the output of the merged annotation/samples expression files.\n")
          }
          
          if(userParameters$createAnno && userParameters$filtering) {
                    ##filtered data with anno
                    cat("\nMerging annotation file with the filtered normalized expression data\n")
                    
                    anno.filtData = createAnnoFun(filtered.normData, userParameters$lib.mapping, userParameters$lib.All.mapping);
                    eset.anno.filtData = cbind(anno.filtData, filtered.normData);
                    
                    fileName <- paste(userParameters$outputDir, userParameters$studyName, "_normData_Filtered_", userParameters$normType, ".txt", sep="")
                    
                    write.table(eset.anno.filtData, file = fileName, quote= FALSE, sep='\t', row.names= F, col.names= T )
                    
                    cat("\nSaving merged filtered normalized data completed.\n")
          }else{
                    cat("\nSkipping the output of the merged filtered normalized annotation/samples expression files.\n")
          }
          
          #Save the R history
          if(userParameters$saveHistory){
                    fileName <- paste(userParameters$outputDir, userParameters$studyName, ".Rhistory", sep="")
                    savehistory(fileName)
                    cat("\nSaved Rhistory in: ", fileName,"\n", sep="")
          }
          
          #Yay, everything worked! (Or at least no errors were found ;) ) 
          message = ("\nIllumina Beadchip normalisation (and/or QC) has been completed!!\n")
          x<- cat(message)
          
          message <- paste("Completed in:", ceiling( (proc.time() -ptm)[3]/60) ,"minutes", sep=" ")
          cat(message)
          
          #Stop recording log
          if(userParameters$createLog) sink()
          
}#End function runIlluminaNormalisation
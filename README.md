localIlluminaNormalization
==========================

###Pipeline to perform local illuminina normalization, QC and statistics.

Scripts for locally running background correction, normalization, quality assessment and statistics of Illumina Beadchip data using the methods found in the limma and lumi packages. If possible, missing packages are installed.

Scripts based on arrayAnalysis developed by: Department of Bioinformatics - BiGCaT Bioinformatics and Systems Biology Research Group Maastricht University - The Netherlands 

The main files which are needed are:
- Sample Probe Profile (`-s/--sampleProbeProfilePath`)
  - Contains the expression values of the samples on the array
- Control Probe Profile (`-c/--controlProbeProfilePath`)
  - Contains the expression values of the control probes, used when background correcting
- descriptionFile	(`-d/--descFile`)
  - (Tab-delimited file containing: arrayNames (Name of the array as defined in the Sample Probe profile) | sampleNames (User specified name for the sample) | sampleGroup (If samples should be clustered, provide identical identifiers here) )
- statSubsetFile  *(Optional)* (`-f/--statFile`)
  - This file contains the sampleNames that should be used in the subset where statistics are run on if -B TRUE. The sampleNames should correspond with the sampleNames in the descriptionFile.
- R Lumibatch object  *(Optional)* (`-m/--normData`)
  - An R LumiBatch object from a previous normalization containing the normalized expression values if -f is true.	

##Usage
The pipeline takes in command-line arguments from Rscript or RStudio. (Hardcoded commandArgs call with a vector)

The scripts can be run using:
- Rscript:
  - `Rscript runIlluminaNormalization -i \<inputDir\> -o \<outputDir\> -s \<sampleProfile.txt\> -c  \<controlProfile.txt\> -d  \<descriptionFile.txt\> etc. `
- RStudio:
  - Arguments can be given by adding this line in `runIlluminaNormalization.R` and supplying a vector of arguments: `userParameters <- getArguments(c("-o", "~/outputDir", "-i" ,"~/inputDir","-s","Sample_Probe_Profile.txt","-c", "Control_Probe_Profile.txt","-d","descFile.txt"))`

##Parameters
All ~55 possible arguments can also be seen by running runIlluminaNormalization.R -h.

Useful short flags:

- `-s`  Name of the Sample_Probe_Profile
- `-c`  Name of the Control_Probe_Profile
- `-d`  Name of the descriptionFile (Tab-delimited file containing: arrayNames | sampleNames | sampleGroup)
- `-N`  Name of the study (Used in the naming of files)
- `-l`  Create a log file in the output directory (Sinks stdout and stderr) (Useful when running as a background process or in batch modus)
- `-i`  Input directory (Where the files are found, if not the same as )
- `-o`  Output directory (Where the files are stored)
- `-p`  Whether to perform plots + PCA on data. (TRUE/FALSE)
- `-n`  Whether to normalize or not (TRUE/FALSE)
- `-u`  Whether to skip background correction (TRUE/FALSE) (TRUE means skip background correction)
- `-m`  Name of the file containing the normalized data R object   (When already normalized and wanting to do new statistics)
- `-F`  File containing the sampleNames/AssayNames on which statistics should be done.
- `-B`  Whether to perform statistics on a subset (As defined in -F)
- `-f`  Whether to load the old normalized data, thereby skipping normalization and only perform statistics.

Full list of parameters:

Options:

	-i INPUTDIR, --inputDir=INPUTDIR
		Path to folder where the Control_Probe_Profile, Sample_Probe_Profile and Description file are found 
    default = [\<path where script is run\>] 

	-o OUTPUTDIR, --outputDir=OUTPUTDIR
		Path to folder where the output files will be stored 
    default = [\<path where script is run\>] 

	-a ANNODIR, --annoDir=ANNODIR
		Path to folder where the annotation files for the arrays are stored. (Containing the Probes + Genes etc. on the array) 
	  default = [\<path where script is run\>] 

	--scriptDir=SCRIPTDIR
		Path to folder where the scripts are stored. 
	  default = [\<path where script is run\>] 

	--species=SPECIES
		Species used for samples (E.g. Human, Mouse or Rat etc) 
	  default = [Human] 

	--arrayType=ARRAYTYPE
		Main type of array that was used (E.g. HumanHT-12 or HumanWG-6) 
  	default = [HumanHT-12] 

	--annoType=ANNOTYPE
		What annotation file should be used to check Genes+Probes etc (E.g. humanHT-12_V4_0_R2_15002873_B) 
	  default = [HumanHT-12_V4_0_R2_15002873_B] 

	-N STUDYNAME, --studyName=STUDYNAME
		Used to be called ns, Name of the study (Used in naming the output files) 
	  default = [2014-03-17_11.52.13] 

	--createAnno=CREATEANNO
		Create an annotation file containing lumiID, probeID etc. 
	  default = [TRUE] 

	--saveHistory=SAVEHISTORY
		Whether to save the R history to the output directory. 
	  default = [TRUE] 

	-l CREATELOG, --createLog=CREATELOG
		Create a log file in the output directory. Captures ALL messages from stdout and stderr. 
	  default = [FALSE] 

	-B STATSUBSET, --statSubset=STATSUBSET
		Whether statistics should be done only on a subset. (defined in --statFile) 
	  default = [FALSE] 

	-f LOADOLDNORM, --loadOldNorm=LOADOLDNORM
		Whether to load the old normalized data (given with -m/--normData). 
		default = [FALSE] 

	-s SAMPLEPROBEPROFILEPATH, --sampleProbeProfilePath=SAMPLEPROBEPROFILEPATH
		Used to be called expFile, contains the expression values of the samples 
	  default = [Sample_Probe_Profile.txt] 

	-c CONTROLPROBEPROFILEPATH, --controlProbeProfilePath=CONTROLPROBEPROFILEPATH
		Used to be called bgFile, contains the expression values of the control probes 
	  default = [Control_Probe_Profile.txt] 

	-d DESCFILE, --descFile=DESCFILE
		Tab-delimited file containing: arrayNames | sampleNames | sampleGroup 
	  default = [descriptionFile.txt] 

	-F STATFILE, --statFile=STATFILE
		File containing a single column with the sampleNames or ArrayNames on which statistics should be performed. If none given, statistics is performed on all samples in the description file. 
  	default = [statSubsetFile.txt] 

	-m NORMDATA, --normData=NORMDATA
		R Lumibatch object containing the normalized expression values. 
  	default = [normData.Rdata] 

	-u BGSUB, --bgSub=BGSUB
		Whether background correction has been done already. (E.g. in Illumina GenomeStudio) 
		If FALSE -> perform bg correction using controlProbeProfileFile. If TRUE -> Skip bg correction 
		default = [FALSE] 

	--detectionTh=DETECTIONTH
		The p-value threshold of determining detectability of the expression. 
		default = [0.01] 

	--convertNuID=CONVERTNUID
		Determine whether to convert the probe identifier as nuID 
		default = [TRUE] 

	--dec=DEC
		The character used in the files to indicate decimal values (E.g. '.' or ',' etc.) 
		default = [.] 

	--parseColumnName=PARSECOLUMNNAME
		Determine whether to parse the column names and retrieve the sample information 
		(Assume the sample information is separated by a tab) 
		default = [FALSE] 

	--checkDupId=CHECKDUPID
		Determine whether to check duplicated TargetIDs or ProbeIds. The duplicated ones will be averaged. 
		default = [TRUE] 

	--save.rawData=SAVE.RAWDATA
		Whether to save lumi.batch of raw data as R object in WORK.DIR 
		default = [TRUE] 

	--normType=NORMTYPE
		Type of normalization to do (lumi or neqc) 
		default = [lumi] 

	--bg.correct=BG.CORRECT
		Whether to do a background correction using the controlProbeProfile (If not done already) 
  	default = [TRUE] 

	--bgcorrect.m=BGCORRECT.M
		List of the parameters for lumiExpresso{lumiB}, method for background correction. 
		Possible parameters: c('none', 'bgAdjust', 'forcePositive', 'bgAdjust.affy') 
  	default = [bgAdjust] 

	--variance.stabilize=VARIANCE.STABILIZE
		Whether do to variance stabilization 
	default = [TRUE] 

	--variance.m=VARIANCE.M
		Method of variance stabilization for lumiExpresso{lumiT} package.
		Possible parameters are: c('vst', 'log2', 'cubicRoot') 
  	default = [log2] 

	-n NORMALIZE, --normalize=NORMALIZE
		Whether to normalize or not 
	  default = [TRUE] 

	--normalization.m=NORMALIZATION.M
		List of parameters for lumiExpresso{lumiB}, method of normalization. 
		parameters are: c(quantile, rsn, ssn, loess, vsn) 
  	default = [quantile] 

	--save.normData=SAVE.NORMDATA
		Whether to save lumi.batch of normalized data as R object in WORK.DIR 
  	default = [TRUE] 

	--filtering=FILTERING
		Whether filtering of probes with a low/no expression should be performed 
  	default = [TRUE] 

	--filter.Th=FILTER.TH
		Threshold for probe filtering. 
	  default = [0.01] 

	--filter.dp=FILTER.DP
		Threshold for the count of the same probe in multiple samples 
	  default = [0] 

	-p PERFORMSTATISTICS, --performStatistics=PERFORMSTATISTICS
		Should clustering and PCA be done alongside the normalization of the data? (Individual options are below) 
	  default = [TRUE] 

	--rawDataQC=RAWDATAQC
		Determine whether to do QC assessment for the raw data; if false no summary can be computed. 
	  default = [TRUE] 

	--normDataQC=NORMDATAQC
		Determine whether to do QC assessment for the normed data; if false no summary can be computed. 
  	default = [TRUE] 

	--rawSummary=RAWSUMMARY
		Whether to create a summary table in WORK.DIR of the raw data 
	  default = [TRUE] 

	--normSummary=NORMSUMMARY
		Whether to create a summary table in WORK.DIR of the normalized data 
	  default = [TRUE] 

	--perGroup=PERGROUP
		Reorder rawData lumibatch file FIRST on Group and THEN ON sampleNames (Used for the order of visualization in plots) 
	  default = [TRUE] 

	--raw.boxplot=RAW.BOXPLOT
		Should a boxplot be made for the raw data. 
	  default = [TRUE] 

	--raw.density=RAW.DENSITY
		Should a density plot be made for the raw data. 
	  default = [TRUE] 

	--raw.cv=RAW.CV
		Should a cv plot be made for the raw data. 
	  default = [TRUE] 

	--raw.sampleRelation=RAW.SAMPLERELATION
		Should a sample relation plot be made for the raw data. 
    default = [TRUE] 

	--raw.pca=RAW.PCA
		Should a PCA be made for the raw data. 
	  default = [TRUE] 

	--raw.correl=RAW.CORREL
		Should a correlation plot be made for the raw data. 
	  default = [TRUE] 

	--norm.boxplot=NORM.BOXPLOT
		Should a boxplot be made for the normalized data. 
	  default = [TRUE] 

	--norm.density=NORM.DENSITY
		Should a density plot be made for the normalized data. 
	  default = [TRUE] 

	--norm.cv=NORM.CV
		Should a cv plot be made for the normalized data. 
	  default = [TRUE] 

	--norm.sampleRelation=NORM.SAMPLERELATION
		Should a sample relation plot be made for the normalized data. 
	  default = [TRUE] 

	--norm.pca=NORM.PCA
		Should a PCA be made for the normalized data. 
	  default = [TRUE] 

	--norm.correl=NORM.CORREL
		Should a correlation plot be made for the normalized data. 
	  default = [TRUE] 

	--clusterOption1=CLUSTEROPTION1
		Distance calculation method. 
    Possible parameters are c('Pearson', 'Spearman', 'Euclidean') 
	  default = [Pearson] 

	--clusterOption2=CLUSTEROPTION2
		Method of clustering. 
    Possible parameters are c('Ward', 'McQuitty', 'average', 'median', 'single', 'complete', 'centroid') 
	  default = [average] 

	--img.width=IMG.WIDTH
		The max. width of the plots. 
	  default = [1920] 

	--img.heigth=IMG.HEIGTH
		The max. heigth of the plots. 
	  default = [1080] 

	--img.pointSize=IMG.POINTSIZE
		The size of the points on plots. 
	  default = [24] 

	--img.maxArray=IMG.MAXARRAY
		The maximum datapoint on each plot per page. 
	  default = [41] 

	-h, --help
		Show this help message and exit

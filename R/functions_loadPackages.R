#Author:                      Varshna Goelela
#Date of  creation:           26-2-14
#Date of modification:        26-2-14 (Job van Riet)
#Version:                     1.1
#Modifications:               Original version
#Known bugs:                  None known
#Function:                    This script contains the functions load mandatory packages and install them if they are not found on the system.

###############################################################################
# calls function checkInstallPackages                                         #
# loads all the mandatory library packages                                    #
# Needs a character list which contains names of R packages                   #
###############################################################################

loadPackages <- function(man) {
          #Function to check if all mandatory packages are installed
          checkInstalledPackages(man)
          cat("Loading library packages:\n")
          #For each factor (package) in fMan load package in quiet mode
          for (i in 1:length(man) ) {
                    library(man[i] , character.only = TRUE)
                    print( paste (" '",man[i] ,"'" , " has been loaded", sep= "") )
          }
}

###############################################################################
# check installation of mandatory packages                                    #
# perform installation of mandatory packages                                  #
###############################################################################

checkInstalledPackages <- function(man) {
          
          #Get a list of all the installed packages
          inst<-installed.packages()
          
          #Check if mandatory packages are listed in installed list
          notInst<-man[!(man %in% inst[,1])]
          
          #If the # of installed packages != the same as the # of defined mandatory packages
          if ( length( notInst)!=0) {
                    print("Some of the mandatory packages are missing, we are now going to install them", sep="")
                    
                    #Install packages that are not installed
                    source("http://www.bioconductor.org/biocLite.R")
                    biocLite( notInst)
                    
                    #Get a new list of all the installed packages
                    inst<-installed.packages()
                    
                    #Check for second time if mandatory packages are listed in new installed list
                    notInst2 <- man[!(man %in% inst[,1])]
                    
                    if ( length( notInst2)==0) {
                              print( paste("Ok, required packages are installed and loaded: ", man, sep=" ") )
                    } 
                    else{
                              res <- ( paste("!!!Check if followwing package:", notInst2, 
                                             "is available in the repository.",
                                             "And try a manual install of the package", sep= " ") )
                    }
          } 
          else{
                    res <-( paste("Ok, required package installed:" , man, sep=" ") )
          }
          return(res)
}
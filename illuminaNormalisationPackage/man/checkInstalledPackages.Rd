\name{checkInstalledPackages}
\alias{checkInstalledPackages}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
checks if the packages that are given as parameters are all installed, if not try to install these packages.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
checkInstalledPackages(man)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{man}{
%%     ~~Describe \code{man} here~~
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
function (man) 
{
    inst <- installed.packages()
    notInst <- man[!(man \%in\% inst[, 1])]
    if (length(notInst) != 0) {
        print("Some of the mandatory packages are missing, we are now going to install them", 
            sep = "")
        source("http://www.bioconductor.org/biocLite.R")
        biocLite(notInst)
        inst <- installed.packages()
        notInst2 <- man[!(man \%in\% inst[, 1])]
        if (length(notInst2) == 0) {
            print(paste("Ok, required packages are installed and loaded: ", 
                man, sep = " "))
        }
        else {
            res <- (paste("!!!Check if followwing package:", 
                notInst2, "is available in the repository.", 
                "And try a manual install of the package", sep = " "))
        }
    }
    else {
        res <- (paste("Ok, required package installed:", man, 
            sep = " "))
    }
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{createAnnoFun}
\alias{createAnnoFun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates an annotation file of the raw/normed data containing the expressions and added entrez gene names etc.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
createAnnoFun(df, lib.mapping, lib.All.mapping)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
%%     ~~Describe \code{df} here~~
}
  \item{lib.mapping}{
%%     ~~Describe \code{lib.mapping} here~~
}
  \item{lib.All.mapping}{
%%     ~~Describe \code{lib.All.mapping} here~~
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
function (df, lib.mapping, lib.All.mapping) 
{
    nuIDs = rownames(df)
    annoFile = as.data.frame(nuIDs)
    if (require(lib.All.mapping, character.only = TRUE) & require(annotate)) {
        annoFile$ILMN_GENE <- nuID2IlluminaID(nuIDs, lib.mapping, 
            idType = c("Symbol"))
        annoFile$geneSymbol <- getSYMBOL(nuIDs, lib.All.mapping)
        annoFile$PROBE_ID <- nuID2IlluminaID(nuIDs, lib.mapping, 
            idType = c("Probe"))
        annoFile$ENTREZ_GENE_ID <- unlist(lookUp(nuIDs, lib.All.mapping, 
            "ENTREZID"))
        annoFile$GENE_NAME <- sapply(lookUp(nuIDs, lib.All.mapping, 
            "GENENAME"), function(x) x[1])
        annoFile$ACCESSION <- unlist(lookUp(nuIDs, lib.All.mapping, 
            "ACCNUM"))
        annoFile$PROBE_SEQUENCE <- unlist(id2seq(nuIDs))
    }
    return(annoFile)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

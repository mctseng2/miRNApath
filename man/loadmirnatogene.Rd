\name{loadmirnatogene}
\alias{loadmirnatogene}
\title{Load miRNA to gene associations for miRNApath}
\description{
This method loads associations between miRNAs to the genes they
affect.
}
\usage{
loadmirnatogene(mirnafile, mirnaobj, mirnacol="miRNA Name",
genecol="Entrez Gene ID", columns=NA)
}
\arguments{
 \item{mirnafile}{
     The tab-delimited miRNA results file to be loaded. The file
     is expected to be in tall-skinny format.
 }
 \item{mirnaobj}{
     An object of type mirnapath containing data resulting from
     the \code{loadmirnapath} method.
 }
 \item{mirnacol}{
     The name of the column header which contains the miRNA names
     being assayed. That is, the name of the column header in the
     file being read.
 }
 \item{genecol}{
     The name of the column header which contains the gene names
     being assayed.
 }
 \item{columns}{
     The names of any additional columns in the file being read
     which should equate with the mirnapath object.
 }
}
\details{
  The data is expected to have miRNA names which exactly match
  those in the mirnaTable item of the mirnapath object. Also, the
  gene names are expected to match exactly with those gene names
  loaded by \code{loadmirnapathways}.
}
\value{
The method returns an object of type mirnapath, a list with
components:

  \item{ mirnaTable }{
      data.frame containing the miRNA results data
  }
  \item{ columns }{
      list containing the names of required column headers
      associated to the actual column header supplied in the
      dataset contained in mirnaTable. Required headers:
      mirnacol, assayidcol. Optional headers: groupcol,
      pvaluecol, foldchangecol, expressioncol,
      filterflagcol
  }
  \item{ groupcount }{
      the number of groups contained in mirnaTable using the
      groupcol, if supplied
  }
  \item{ state }{
      the current state of the object, using the following
      values in order of progress through the typical workflow:
      unfiltered, filtered, enriched.
  }
}
\references{
    John Cogswell (2008) Identification of miRNA changes
    in Alzheimer's disease brain and CSF yields putative 
    biomarkers and insights into disease pathways, Journal of
    Alzheimer's Disease 14, 27-41.
}
\author{ James M. Ward \email{jmw86069@gmail.com} }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{
    \code{\link{loadmirnapath}}
    \code{\link{filtermirnapath}}
    \code{\link{loadmirnatogene}}
    \code{\link{loadmirnapathways}}
}
\examples{

## Load miRNA expression data from AD miRNA paper
## This data contains miRNA expression data, 
data(mirnaobj);

## Display the state, which should generally be "unfiltered"
## at this point
mirnaobj@state;

## Display summary information about the object
mirnaobj;

## Annotate hits by filtering by P-value 0.05
mirnaobj <- filtermirnapath( mirnaobj, pvalue = 0.05,
    expression = NA, foldchange = NA );

## Write a file as example of required input
write.table(mirnaobj@mirnaGene, file = "mirnaGene.txt", 
    quote = FALSE, row.names = FALSE, col.names = TRUE, na = "",
    sep = "\t");

## Load the miRNA to gene associations
mirnaobj <- loadmirnatogene( mirnafile = "mirnaGene.txt",
    mirnaobj = mirnaobj, mirnacol = "miRNA Name",
    genecol = "Entrez Gene ID", 
    columns = c(assayidcol = "ASSAYID") );

## Display summary, noting the number of genes reported
mirnaobj;
}

\keyword{ IO }
\keyword{ manip }
\keyword{ attribute }

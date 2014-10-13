\name{loadmirnapathways}
\alias{loadmirnapathways}
\title{Load gene to pathway associations for miRNApath}
\description{
 This method loads associations between genes and the pathways to which they belong.
}
\usage{
loadmirnapathways(mirnaobj, pathwayfile, genecol="Entrez Gene ID",
pathwaycol="PATHWAY", columns=c(), pathwayidcol=NA)
}
\arguments{
 \item{mirnaobj}{
     An object of type mirnapath containing data resulting from
     the \code{loadmirnapath} method.
 }
 \item{pathwayfile}{
     The file containing the gene to pathway associations.
 }
 \item{genecol}{
     The name of the column header which contains the gene names
     associated with pathway data.
 }
 \item{pathwaycol}{
     The name of the column header which contains the pathway
     names.
 }
 \item{columns}{
     The names of any additional columns in the file being read
     which should equate with the mirnapath object.
 }
 \item{pathwayidcol}{
     The (optional) column header for IDs associated with the
     pathway names.
 }
}
\details{
  The data loaded is expected to have gene names which exactly
  match those gene names loaded by \code{loadmirnatogene}.
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
\seealso{
    \code{\link{loadmirnapath}},
    \code{\link{filtermirnapath}},
    \code{\link{loadmirnatogene}},
    \code{\link{loadmirnapathways}},
    \code{\link{runEnrichment}}
}
\examples{

## Load miRNA expression data from AD miRNA paper
## This data contains miRNA expression data, 
data(mirnaobj);

## Write a file as example of required input
write.table(mirnaobj@mirnaPathways, file = "mirnaPathways.txt", 
    quote = FALSE, row.names = FALSE, col.names = TRUE, na = "",
    sep = "\t");

## Load the gene to pathway associations
mirnaobj <- loadmirnapathways( mirnaobj = mirnaobj, 
    pathwayfile = "mirnaPathways.txt", 
    pathwaycol = "Pathway Name", genecol = "Entrez Gene ID");

## Display summary, noting the number of genes reported
mirnaobj;
}

\keyword{ IO }
\keyword{ manip }
\keyword{ attribute }

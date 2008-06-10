\name{loadmirnapath}
\alias{loadmirnapath}
\title{Load miRNApath Data}
\description{
    This method loads data from a tab-delimited flatfile into an
    object of type mirnapath to be used for further miRNA
    analysis.
}
\usage{
loadmirnapath(mirnafile="mirna_input.txt", mirnacol="miRNA Name",
    assayidcol="ASSAYID", groupcol="GROUP", 
    filterflagcol="FILTERFLAG", expressioncol=NA, 
    foldchangecol=NA, pvaluecol=NA)
}
\arguments{
 \item{mirnafile}{
     The tab-delimited miRNA results file to be loaded. The file is expected to be in tall-skinny format.
 }
 \item{mirnacol}{
     The name of the column header which contains the miRNA names
     being assayed.
 }
 \item{assayidcol}{
     The name of the column containing values which distinguish
     different assays for the same miRNA.
 }
 \item{groupcol}{
     The (optional) name of the column which contains sample
     group information. Enrichment is run separately for each
     sample group, defining a unique universe for the basis of
     the enrichment.
 }
 \item{filterflagcol}{
     The column header which does or will contain a flag
     distinguishing hits from non-hits. This column is typically
     not supplied and is created during the 
     \code{filtermirnapath} step.
 }
 \item{expressioncol}{
     The (optional) column header for values containing the
     expression abundances of the miRNAs assayed.
 }
 \item{foldchangecol}{
     The (optional) column header for values containing the fold changes
     of the miRNAs assayed.
 }
 \item{pvaluecol}{
     The (pvaluecol) column header for values containing the
     P-values of the miRNAs assayed.
 }
}
\details{
This method is the primary means for loading data into the
miRNApath package.

Data is not assumed to have any particular numerical values,
however the basic column types are typically used: expression
abundance, fold change, and P-value. Should one or more columns
be specified and available, it will be available for filtering
later on with \code{\link{filtermirnapath}}.

The group column assumes there is one column containing all
sample group information.

The assayid column is used to distinguish multiple assays for the
same miRNA, such as different vendors, or even different
preparations of the same miRNA assay.
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
    \code{\link{loadmirnapathways}}
}
\examples{

## Start with miRNA data from this package
data(mirnaobj);

## Write a file as example of required input
write.table(mirnaobj@mirnaTable, file = "mirnaobj.txt", 
    quote = FALSE, row.names = FALSE, col.names = TRUE, na = "",
    sep = "\t");

## Now essentially load it back, but assign column headers
mirnaobj <- loadmirnapath( mirnafile = "mirnaobj.txt",
    pvaluecol = "P-value", groupcol = "GROUP", 
    mirnacol = "miRNA Name", assayidcol = "ASSAYID" );

## Display summary information for the resulting object
mirnaobj;
}

\keyword{ IO }
\keyword{ manip }

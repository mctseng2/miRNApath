\name{filtermirnapath}
\alias{filtermirnapath}
\title{ Filter miRNApath data }
\description{
    This method filters the miRNApath data to denote hits versus
    non-hits, the required distinction for running the enrichment
    algorithm. Data is expected to have been loaded by the
    \code{loadmirnapath} method.
}
\usage{
filtermirnapath(mirnaobj, pvalue=NA, expression=NA, foldchange=NA)
}
\arguments{
 \item{mirnaobj}{
    An object of type mirnapath containing data resulting from
    the \code{loadmirnapath} method.
 }
 \item{pvalue}{
    If a p-value column has been defined in the mirnapath object,
    this value is used to define a subset of entries within the
    dataset which will be denoted as hits.
 }
 \item{expression}{
    If an expression column has been defined in the mirnapath
    object, this value will be used to define entries with
    expression above this expression level as hits.
 }
 \item{foldchange}{
    If a fold change column has been defined in the mirnapath
    object, this value is used to require hits to have a fold
    change greater than or equal to this value. The fold change
    is evaluated in both the positive and the negative, such
    that a foldchange=2 will allow foldchange=2 and
    foldchange=-2.
 }
}
\details{
    This method takes a mirnapath object and assigns a flag for
    hits and non-hits, depending upon what filter criteria was
    provided. If multiple criteria are provided, they will all be
    collectively applied such that all criteria must be
    fulfilled. To that end, multiple calls to this function
    on the same mirnapath object should successively shrink the
    list of hits dependent upon the given criteria.
}
\value{
  Object of type mirnapath. The state of the object will reflect
  that the data has been filtered.
}
\references{
    John Cogswell (2008) Identification of miRNA changes
    in Alzheimer's disease brain and CSF yields putative 
    biomarkers and insights into disease pathways, Journal of
    Alzheimer's Disease 14, 27-41.
}
\author{ James M. Ward \email{jmw86069@gmail.com} }
\note{
    The method attempts to convert fold change columns
    appropriately so that filtering by 2 will properly mark
    entries greater than 2, less than -2, or less than 0.5,
    as the case may be.
}
\seealso{
    \code{\link{loadmirnapath}},
    \code{\link{filtermirnapath}},
    \code{\link{loadmirnatogene}},
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

## Display summary, noting the state is "filtered"
mirnaobj;
}

\keyword{ manip }

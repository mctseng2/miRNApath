\name{miRNApath-package}
\alias{miRNApath-package}
\docType{package}
\title{
miRNApath: Pathway Enrichment for miRNA Expression Data

}
\description{
This package provides methods for assessing the statistical 
over-representation of miRNA effects on gene sets, using
supplied miRNA-to-gene associations. Because these associations
are notably many-to-many (one miRNA to many genes; one gene
affected by many miRNAs) the assessment is complex and warrants
perhaps different approaches than are classically performed on
differential gene expression datasets.
}
\details{
\tabular{ll}{
Package: \tab miRNApath\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2008-04-02\cr
License: \tab LGL-2.1, see COPYING.LIB\cr
}
}
\author{
James M. Ward

Maintainer: James M. Ward <jmw86069@gmail.com>

}
\references{
    John Cogswell (2008) Identification of miRNA changes
    in Alzheimer's disease brain and CSF yields putative 
    biomarkers and insights into disease pathways, Journal of
    Alzheimer's Disease 14, 27-41.
}
\keyword{ package }
\keyword{ htest }
\seealso{
    \code{\link{loadmirnapath}},
    \code{\link{filtermirnapath}},
    \code{\link{loadmirnatogene}},
    \code{\link{loadmirnapathways}},
    \code{\link{runEnrichment}}
}
\examples{
\dontrun{
## Start with miRNA data from this package
data(mirnaobj);

## Write a file as example of required input
write.table(mirnaobj@mirnaTable, file = "mirnaTable.txt", 
    quote = FALSE, row.names = FALSE, col.names = TRUE, na = "",
    sep = "\t");

## Now essentially load it back, but assign column headers
mirnaobj <- loadmirnapath( mirnafile = "mirnaTable.txt",
    pvaluecol = "P-value", groupcol = "GROUP", 
    mirnacol = "miRNA Name", assayidcol = "ASSAYID" );

## Start with miRNA data from this package
data(mirnaobj);

## Write a file as example of required input
write.table(mirnaobj@mirnaGene, file = "mirnaGene.txt", 
    quote = FALSE, row.names = FALSE, col.names = TRUE, na = "",
    sep = "\t");

## Load the miRNA to gene associations
mirnaobj <- loadmirnatogene( mirnafile = "mirnaGene.txt",
    mirnaobj = mirnaobj, mirnacol = "miRNA Name",
    genecol = "Entrez Gene ID", 
    columns = c(assayidcol = "ASSAYID") );

## Write a file as example of required input
write.table(mirnaobj@mirnaPathways, file = "mirnaPathways.txt", 
    quote = FALSE, row.names = FALSE, col.names = TRUE, na = "",
    sep = "\t");

## Load the gene to pathway associations
mirnaobj <- loadmirnapathways( mirnaobj = mirnaobj, 
    pathwayfile = "mirnaPathways.txt", 
    pathwaycol = "Pathway Name", genecol = "Entrez Gene ID");

## Annotate hits by filtering by P-value 0.05
mirnaobj <- filtermirnapath( mirnaobj, pvalue = 0.05,
    expression = NA, foldchange = NA );

## Now run enrichment test
mirnaobj <- runEnrichment( mirnaobj=mirnaobj, Composite=TRUE,
   groups=NULL, permutations=0 );

## Print out a summary table of significant results
finaltable <- mirnaTable( mirnaobj, groups=NULL, format="Tall", 
    Significance=0.1, pvalueTypes=c("pvalues") );
finaltable[1:4,];

## Example which calls heatmap function on the resulting data
widetable <- mirnaTable( mirnaobj, groups=NULL, format="Wide", 
    Significance=0.1, na.char=NA, pvalueTypes=c("pvalues") );
## Assign 1 to NA values, assuming they're all equally
## non-significant
widetable[is.na(widetable)] <- 1;

## Display a heatmap of the result across sample groups
pathwaycol <- mirnaobj@columns["pathwaycol"];
pathwayidcol <- mirnaobj@columns["pathwayidcol"];
rownames(widetable) <- apply(widetable[,c(pathwaycol,
   pathwayidcol)], 1, function(i)paste(i, collapse="-"));
wt <- as.matrix(widetable[3:dim(widetable)[2]], mode="numeric")
heatmap(wt, scale="col");

## Show results where pathways are shared in four or more
## sample groups
pathwaySubset <- apply(wt, 1, function(i)
{
   length(i[i < 1]) >= 4;
} )
heatmap(wt[pathwaySubset,], scale="row");
}
}

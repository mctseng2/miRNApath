\name{mirnaTable}
\alias{mirnaTable}
\title{ Create miRNA Enrichment Summary Table as data.frame }
\description{
  This function takes an miRNApath object which has been
  evaluated by runEnrichment(), and provides a data.frame
  summary.
}
\usage{
mirnaTable(mirnaobj, groups=NULL, format="Tall", 
   Significance=0.2, na.char=NA, pvalueTypes=c("pvalues",
   "permpvalues"), maxStringLength=NA)
}
\arguments{
 \item{mirnaobj}{
     An object of type mirnapath containing data resulting from
     the \code{loadmirnapath} method.
 }
 \item{groups}{
     List of groups to include in the data.frame, or NULL to
     include all groups in the miRNApath object.
 }
 \item{format}{
     This parameter tells the method to return "Tall",
     "SuperTall", or "Wide" data. See details below for a
     description of each format.
 }
 \item{Significance}{
     A numerical value specifying the P-value cutoff to use to
     subset the data returned in the data.frame. To avoid
     subsetting the data, provide a value of 1.
 }
 \item{na.char}{
     Value to use for NA instead of leaving NA as-is, potentially
     useful for text output.
 }
 \item{pvalueTypes}{
     Defines which P-value columns should be returned, more
     useful for the Wide format which could otherwise have two
     sets of P-value columns if permutation adjustment were
     used.
 }
 \item{maxStringLength}{
     Defines the maximum length per character string, after
     being determined by \code{nchar}. Strings and
     column headers are both truncated to this length.
 }
}
\details{
  This function simply combines the various results from the
  runEnrichment method into one data.frame suitable for plotting
  or printing in a table. Due to potentially large data volume,
  the subset feature even when used liberally can substantially
  reduce the returned dataset size.

  The \code{maxStringLength} value is particularly useful, often
  critical, for displaying a summary table in text format, since
  pathway names sample group names can be quite long.
  Although there is no default, a recommended value of 50 seems
  to fit the appropriate balance of being short enough to fit
  within a table, and yet be long enough to describe the
  pathway. The Wide format will contain sample group names as
  column headers, and a value of 50 should not in theory affect
  the name, except where it wouldn't be readable in a table
  anyway.
}
\value{
  data.frame

  For Tall data, the columns contain P-values and other values
  useful for discriminating potential hits, the rows contain each
  miRNA-group combination tested which meets the P-value cutoff.
  The miRNAs and genes contributing to the enrichment results
  are concatenated to be summarized in one row and can be rather
  large.

  For SuperTall data, the Tall table as described above is
  returned, except that the concatenated miRNA-gene values are
  separated to one row each. Every individual miRNA and gene
  value is represented on its own row, which can facilitate some
  summary views or data filtering techniques (e.g. Excel or
  Spotfire.)

  For Wide data, the columns contain the group names, the rows
  contain the pathway name, and the cells contain the P-value.
  Note that the column names will have the P-value column header
  prepended to the column name, e.g. "pvalue.GroupName".

  An important note when supplying string na.char values, be
  sure to convert the data to a numeric matrix before calling
  functions such as heatmap, taking care to remove
  string values or convert strings to 1.0 beforehand.
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
    \code{\link{runEnrichment}},
}
\examples{
## Start with miRNA data from this package
data(mirnaobj);

## Now run enrichment test
mirnaobj <- runEnrichment( mirnaobj=mirnaobj, Composite=TRUE,
   groups=NULL, permutations=0 );

## Print out a summary table of significant results
finaltable <- mirnaTable( mirnaobj, groups=NULL, format="Tall", 
    Significance=0.1, pvalueTypes=c("pvalues") );
finaltable[1:20,];

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

\keyword{ manip }

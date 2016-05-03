\name{runEnrichment}
\alias{runEnrichment}
\title{ Perform gene set enrichment analysis on a miRNApath 
   object }
\description{
   This method performs a hypergeometric enrichment analysis on
   a miRNApath object, which should already contain miRNA data,
   miRNA-gene associations, gene-pathway associations, and some
   criteria for filtering miRNA hits from the full tested set.
}
\usage{
runEnrichment(mirnaobj, Composite=TRUE, groups=NULL,
permutations=0)
}
\arguments{
 \item{mirnaobj}{
     An object of type mirnapath containing data resulting from
     the \code{loadmirnapath} method.
 }
 \item{Composite}{
     Defines whether the enrichment treats miRNA-gene as the
     enriched entity, or uses only the gene.
 }
 \item{groups}{
     List of groups to include in the analysis, although each
     group is analyzed independent from the other groups.
 }
 \item{permutations}{
     The number of permutations to use in calculating an adusted
     P-value. Value of 0 will perform no permutations.
 }
}
\details{
  The composite flag indicates whether to treat the fully 
  expanded miRNA-gene combinations as separate enrichment events
  (TRUE), or whether to treat all effects on one gene as one
  collective event. The latter case reverts to the classic
  un-ordered hypergeometric enrichment technique.

  However the expansion of combinations is the current method
  chosen to represent the multiple predicted effects of miRNAs
  to one gene, and the predicted effect of one miRNA to multiple
  genes. The algorithm will identify statistically significantly
  enriched results when the combination of these effects is
  greater than would be anticipated by random chance.

  The adjusted P-value is calculated using the rank of
  unadjusted P-values divided by the number of permutations
  minus one (such that the best rank from 1,000 permutations
  yields an adjusted P-value of 0.001.) The default value 0 was
  put in place to save time, since most adjustments resulted in
  stronger "hits" and weaker "non-hits" in terms of pathways
  enriched. Thus the results are not substantially changed, and
  permutation adjustment is saved for the final result set.
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
      mirnacol, assayidcol, groupcol, filterflagcol.
  }
  \item{ groupcount }{
      The number of groups contained in mirnaTable using the
      groupcol, if supplied
  }
  \item{ state }{
      The current state of the object, in this case "enriched".
  }
  \item{ mirnaGene }{
      data.frame containing associations between miRNAs and
      genes.
  }
  \item{ mirnaPathways }{
      data.frame containing gene-pathway associations.
  }
  \item{ pathwaycount }{
      Numerical value indicating how many pathways are
      available in the data, provided for convenience.
  }
  \item{ filters }{
      List of filters applied to the data, which may include:
      "P-value", "Fold change", and/or "Expression".
  }
  \item{ enrichment }{
      Enrichment summary data in the form of a list of elements
      for each sample group (the sample group is the name of
      each element.) Each list element is itself a list with
      enrichment result data for each sample group, as
      independently calculated: "pvalues" - list of P-values
      named by pathway ID. "Measured pathway mirnaGenes" -
      total number of miRNA-gene-pathway combinations measured,
      which gives some idea of the overall coverage of pathways.
      The general point is that miRNAs have the potential to
      cover many genes and pathways. "Total mirnaGenes" - number
      of miRNA-gene combinations represented in the data.
      "Enriched pathway mirnaGenes" - number of miRNA-gene values
      enriched in the pathway tested. "Enriched by miRNA" - list
      of miRNAs involved in the pathway tested, with the list of
      genes in parentheses per miRNA. "Enriched by Gene" - same
      as previous except switching gene and miRNA. "Total
      enriched mirnaGenes" - the total number of miRNA-gene
      values involved in any pathway enrichment (significant or
      not.) The total values are useful when comparing across
      sample groups, looking particularly for groups with few
      changes or those with a uniquely high number of changes.
      Lastly, with permutations > 0 "Permutation P-value" will
      contain the rank-adjusted P-value as described in the
      details section.
  }
  \item{ pathwayList }{
      Named list of pathways contained in the
      mirnaobj\$mirnaPathways object, named by the pathway ID
      values found in the pathwayidcol column. This list
      facilitates converting the data in the enrichment element
      to pathway names, since those values are named by the
      pathway ID to conserve memory.
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
\dontrun{
## Start with miRNA data from this package
data(mirnaobj);

## Now run enrichment test
mirnaobj <- runEnrichment( mirnaobj=mirnaobj, Composite=TRUE,
   groups=NULL, permutations=0 );
}
}

\keyword{ htest }

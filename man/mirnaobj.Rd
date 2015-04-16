\name{mirnaobj}
\alias{mirnaobj}
\docType{data}
\title{ mirnaobj, an S4 object of class "mirnapath" }
\description{
  An example miRNApath data object containing miRNA data, gene 
  and pathway associations. The object represents the end result
  of the miRNApath workflow, and serves a convenient source for
  example data.
}
\usage{data(mirnaobj)}
\format{
  The format is an S4 class \code{"mirnapath"} with
   \code{slotNames} as follows:

  \describe{\item{ \code{mirnaTable} }{
      data.frame containing the miRNA results data, expected to
      contain columns with miRNA name, gene name, and ideally
      some column(s) for filtering hits versus background, e.g.
      fold change, expression abundance, P-value. Once the data
      is filtered (see state below) there will be a column with
      a flag indicating which entries are hits and which are
      considered background. This column is found in 
      \code{mirnaobj@columns["filterflagcolumn"]} and is typically
      \code{"FILTERFLAG"}.
  }
  \item{ \code{columns} }{
      Named list of column headers used throughout the analysis.
      The purpose of the names is partly to retain the original
      headers in the mirnaTable data.frame, and partly to
      coordinate the names with the miRNA-gene and gene-pathway
      tables used later in the analysis. The recognized headers:
      \code{mirnacol, assayidcol, genecol, pvaluecol, foldchangecol,
      pathwaycol, pathwayidcol, groupcol, mirnagene}. See the
      documentation for the mirnapath object type for more
      details about usage.
  }
  \item{ \code{groupcount} }{
      Numerical value indicating how many sample groups are
      available in the data, provided for convenience.
  }
  \item{ \code{state} }{
      Character value indicating the current analysis state,
      with values: \code{"unfiltered"} if results are loaded but not
      yet filtered; \code{"filtered"} if results are loaded and hits
      are defined with the filterflagcol column; \code{"enriched"} if
      the data is loaded, filtered, and analyzed for enrichment.
      One can load mirna-gene and gene-pathway data at any point
      which necessitates using the \code{mirnaobj@mirnaGene} or
      \code{mirnaobj@mirnaPathways} object elements to determine if
      that data has been loaded.
  }
  \item{ \code{mirnaGene} }{
      data.frame containing associations between miRNA and
      genes. The data should contain one miRNA-to-gene
      relationship per row, and should contain only those two
      columns. Additional columns are maintained but ignored.
      Note that one can use any values in the genecol column,
      provided they match exactly with values found in the
      \code{mirnaobj@mirnaPathways} element (see below.) Therefore,
      if desired one can use transcript or gene associations,
      or other integration methods as desired.
  }
  \item{ \code{mirnaPathways} }{
      data.frame containing gene-pathway associations. The data
      should contain only one gene-to-pathway association per
      row of data. The data can have pathway ID values, which
      may facilitate comparisons to pathway databases (and may
      allow substantial data volume reduction if necessary.) If
      there is no pathwayidcol column, then one will be created
      using a numerical assignments of the pathway names. Note
      that this conversion is not sensitive to pathway sources,
      so care should be taken to include pathway source in the
      pathway name if two sources share the same pathway name.
      The same is true for pathway ID values, should they be
      purely numerical and have shared values across pathway
      sources.
  }
  \item{ \code{pathwaycount} }{
      Numerical value indicating how many pathways are
      available in the data, provided for convenience.
  }
  \item{ \code{filters} }{
      List of filters applied to the data, which may include:
      \code{"P-value", "Fold change"}, and/or \code{"Expression"}.
  }
  \item{ \code{enrichment} }{
      Enrichment summary data in the form of a list of elements
      for each sample group (the sample group is the name of
      each element.) Each list element is itself a list with
      enrichment result data for each sample group, as
      independently calculated: \code{"pvalues"} - list of P-values
      named by pathway ID; \code{"Measured pathway mirnaGenes"} -
      total number of miRNA-gene-pathway combinations measured,
      which gives some idea of the overall coverage of pathways.
      The general point is that miRNAs have the potential to
      cover many genes and pathways; \code{"Total mirnaGenes"} - number
      of miRNA-gene combinations represented in the data;
      \code{"Enriched pathway mirnaGenes"} - number of miRNA-gene values
      enriched in the pathway tested; \code{"Enriched by miRNA"} - list
      of miRNAs involved in the pathway tested, with the list of
      genes in parentheses per miRNA; \code{"Enriched by Gene"} - same
      as previous except switching gene and miRNA; \code{"Total
      enriched mirnaGenes"} - the total number of miRNA-gene
      values involved in any pathway enrichment (significant or
      not.) The total values are useful when comparing across
      sample groups, looking particularly for groups with few
      changes or those with a uniquely high number of changes.
  }
  \item{ \code{pathwayList} }{
      Named list of pathways contained in the
      \code{mirnaobj@mirnaPathways} object, named by the pathway ID
      values found in the pathwayidcol column. This list
      facilitates converting the data in the enrichment element
      to pathway names, since those values are named by the
      pathway ID to conserve memory.
  }}

}
\source{
  Journal of Alzheimers Disease 14, 27-41.
}
\references{
    John Cogswell (2008) Identification of miRNA changes
    in Alzheimer's disease brain and CSF yields putative 
    biomarkers and insights into disease pathways, Journal of
    Alzheimer's Disease 14, 27-41.
}
\examples{
## Load the data
data(mirnaobj)

## Print the default summary
mirnaobj;
}
\keyword{datasets}

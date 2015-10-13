setMethod("show", signature(object="mirnapath"),
definition <- function(object)
{
   cat("mirnapath object:\n");
   y = summary(object);
   print(y);
   cat(paste( c("\nColumns specified:", sapply(names(object@columns), function(i){ paste( c(i, " = \"", object@columns[i], "\""), collapse="" )}), "\n"), collapse="\n   "));
   if ("filters" %in% names(object))
   {
      cat(paste( c("Filters Applied:", sapply(names(object@filters), function(i){ paste( c(i, " = \"", object@filters[i], "\""), collapse="" )}), "\n"), collapse="\n   "));
   } else {
      cat("Filters Applied:\n   none\n\n");
   }
   # number of miRNA's
   mirnacount = length(unique(object@mirnaTable[, object@columns["mirnacol"] ]));
   cat("Number of miRNAs:", mirnacount, "\n");

   # number of individual miRNA-gene predictions
   if ("mirnaGene" %in% names(object))
   {
      mirnagenecount = length(levels(as.factor(object@mirnaGene[, object@columns["mirnagene"] ])));
      cat("Number of individual miRNA-gene predictions:", mirnagenecount, "\n");

      # number of genes represented
      genecount = length(levels(as.factor(object@mirnaGene[, object@columns["genecol"] ])));
      cat("Total number of genes represented:", genecount, "\n");
   }

   cat("Number of sample groups:", object@groupcount, "\n");
   cat("Number of pathways:", ifelse(is.null(object@pathwaycount), 0, object@pathwaycount), "\n");
   cat("State:", object@state, "\n");
}
);

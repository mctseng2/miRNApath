`mirnaTable` <-
function
( mirnaobj, groups = NULL,
  format = "Tall",
  Significance = 0.2,
  na.char = NA,
  pvalueTypes = c("pvalues", "permpvalues"),
  maxStringLength = NA
 )
{
   # some column headers
   pathwaycol = mirnaobj@columns["pathwaycol"];
   pathwayidcol = mirnaobj@columns["pathwayidcol"];

   # Make sure pvalueTypes are present in the data
   pvalueTypes <- pvalueTypes[pvalueTypes %in% names(mirnaobj@enrichment[[1]]) ];
   if (length(pvalueTypes) == 0)
   {
      stop("The pvalueTypes argument does not match the enrichment\nP-value columns.");
   }

   # Determine which groups to return
   if (is.null(groups))
   {
      groups = names(mirnaobj@enrichment);
   } else {
      groups = groups[groups %in% names(mirnaobj@enrichment)];
   }

   # Return a data.frame with the relevant information
   fulltable <- do.call(rbind, lapply(as.character(groups), function(group)
   {
      t1 <- data.frame( check.names = FALSE, mirnaobj@enrichment[[group]]
         [c(pvalueTypes, "Measured pathway mirnaGenes", 
         "Enriched pathway mirnaGenes", "Genes Enriched", 
         "miRNAs Enriched", "Total mirnaGenes", "Total filtered mirnaGenes")],
         "Group" = rep(group, length(mirnaobj@enrichment[[group]][[1]]) ) );
      # removed "Enriched by miRNA", "Enriched by Gene" for space considerations

      # Filter by P-value
      pvaltest <- t1[,pvalueTypes] < Significance;
      # If more than one P-value column, let any of them meet the criteria
      if (length(pvalueTypes) > 1)
      {
         pvaltest <- apply(pvaltest, 1, function(i){ any(i) } );
      }
      t1 <- t1[pvaltest,];
      # Add some counts for miRNA-gene hits, instead of very long text strings
#      enrichMiRNACount = sapply(t1[,"Enriched by miRNA"],
#         function(i){nchar(gsub("[^;]", "", i))+1});
#      enrichGeneCount = sapply(t1[,"Enriched by Gene"],
#         function(i){nchar(gsub("[^;]", "", i))+1});
#      enrichEntityCount = sapply(t1[,"Enriched by Gene"],
#         function(i){nchar(gsub("[^;,]", "", i))+1});
#      t1[,"# miRNA enriched"] = enrichMiRNACount;
#      t1[,"# Genes enriched"] = enrichGeneCount;
#      t1[,"# miRNA-Genes enriched"] = enrichEntityCount;
      t1 <- data.frame( check.names=FALSE, t1, pathwayidcol = rownames(t1), "Pathway Name"=mirnaobj@pathwayList[rownames(t1)] );
      colnames(t1)[colnames(t1) == "pathwayidcol"] = pathwayidcol;
      t1[,pathwayidcol] = as.character(t1[,pathwayidcol]);
      t1;
   } ) ) ;
   if (!is.na(maxStringLength))
   {
      fulltable[,!colnames(fulltable) %in% pathwayidcol] = 
      apply(fulltable[,!colnames(fulltable) %in% pathwayidcol], 2, function(i)
      {
         if (all(is.character(i)))
         {
            j = sapply(i, function(k)
            {
               k = substr(k, 1, min(nchar(k), maxStringLength));
            } )
         }
      } )
   }
   if (format == "Tall")
   {
      fulltable;
   } else if (format == "SuperTall")
   {
      # Create a format with mirnaGene values fully expanded
      # useful for bar chart histogram views and for drill-down
      SuperTall <- do.call(rbind, apply(fulltable, 1, function(i)
      {
         pathwayid <- i[pathwayidcol][[1]];
         group <- i["Group"][[1]];
         GMtable <- mirnaobj@enrichment[[group]]["Enriched miRNA-Genes"][[1]][pathwayid][[1]];
         SuperTall1 <- data.frame(check.names=FALSE, do.call(rbind,
               lapply(1:dim(GMtable)[1], function(k2) i) ),
               GMtable);
         SuperTall1;
      } ) )
      SuperTall;
   } else if (format == "Wide")
   {
      widetable <- reshape (fulltable, direction="wide", idvar=pathwayidcol,
            timevar="Group", drop=colnames(fulltable)[!colnames(fulltable) %in% 
            c(pathwayidcol, "Group", pvalueTypes)] );
      # Replace NA with another character if supplied
      if (!is.na(na.char))
      {
         widetable[is.na(widetable)] <- na.char;
      }
      # Add in pathway names using the pathway IDs
      widetable[,pathwaycol] <- mirnaobj@pathwayList[widetable[,1]];
      # Then put the pathway IDs and names first in the table
      widetable <- widetable[,c(pathwayidcol, pathwaycol, colnames(widetable)
         [!colnames(widetable) %in% c(pathwayidcol, pathwaycol)] )];
      widetable;
   }
}


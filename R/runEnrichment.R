`runEnrichment` <-
function
( mirnaobj,
  Composite = TRUE,
  groups = NULL,
  permutations = 0)
{
	# The Composite argument tells whether the assayid and gene id should be concatenated
	# in order to preserve the original library size faithfully

   # The groups argument specifies which subsets of groups to analyze, sending NA will
   # run all groups in the mirnaobj@mirnaTable

   mirnacol <- mirnaobj@columns["mirnacol"];
   genecol <- mirnaobj@columns["genecol"];
   mirnagene <- mirnaobj@columns["mirnagene"];
   groupcol <- mirnaobj@columns["groupcol"];
   pathwaycol <- mirnaobj@columns["pathwaycol"];
   pathwayidcol <- mirnaobj@columns["pathwayidcol"];
   filterflagcol <- mirnaobj@columns["filterflagcol"];
   if ( is.na(mirnacol) ) {
      stop("mirnaobj@columns has no value for mirnacol");
   }
   if ( is.na(genecol) ) {
      stop("mirnaobj@columns has no value for genecol");
   }
   if ( is.na(mirnagene) ) {
      stop("mirnaobj@columns has no value for mirnagene");
   }
   if ( is.na(groupcol) ) {
      stop("mirnaobj@columns has no value for groupcol");
   }
   if ( is.na(pathwaycol) ) {
      stop("mirnaobj@columns has no value for pathwaycol");
   }
   if ( is.na(pathwayidcol) ) {
      stop("mirnaobj@columns has no value for pathwayidcol");
   }
   if ( is.na(filterflagcol) ) {
      stop("mirnaobj@columns has no value for filterflagcol");
   }

   if ( is.null(groups) )
   {
      groups <- levels(as.factor(mirnaobj@mirnaTable[, mirnaobj@columns["groupcol"] ]));
   }

   groupResults <- lapply( groups, function(group)
   {
      # iterate through each group as supplied
      groupData <- mirnaobj@mirnaTable[ mirnaobj@mirnaTable[,groupcol] %in% group , ];
      
      groupmirnas <- unique(as.character(groupData[,mirnacol]));
      groupmirnas <- groupmirnas[ !groupmirnas %in% "" ];
      groupmirnaGene <- as.matrix(mirnaobj@mirnaGene[mirnaobj@mirnaGene[,mirnacol] %in% groupmirnas,]);
      groupgenes <- unique(as.character(groupmirnaGene[,genecol]));
      groupPathways <- as.matrix(mirnaobj@mirnaPathways[ mirnaobj@mirnaPathways[,genecol] %in% groupgenes ,]);
      # Repair the pathwayidcol numerical values
      groupPathways[,pathwayidcol] = gsub("^[ ]+", "", groupPathways[,pathwayidcol]);

      permresults <- lapply( 0:permutations, function(perm)
      {
         if (perm == 0)
         {
            # This iteration is non-random, not part of permutation correction
            # filteredmirnas contains the list of actual miRNA hits,
            # upon which everything else is based.
            filteredmirnas <- unique( as.character( groupData[ groupData[, filterflagcol] == 1, mirnacol] ) );
            filteredmirnas <- filteredmirnas[ !filteredmirnas %in% "" ];
         } else {
            # Randomize the miRNA hits, using the same number of miRNA hits
            # Randomizing the miRNA hits and replacing them in filteredmirnas
            # should provide a reasonable permutation model.
            filteredmirnas <- unique( as.character( groupData[ groupData[, filterflagcol] == 1, mirnacol] ) );
            filteredmirnas <- filteredmirnas[ !filteredmirnas %in% "" ];
            filteredmirnas <- sample(groupmirnas, length(filteredmirnas));
         }
         filteredmirnaGene <- mirnaobj@mirnaGene[mirnaobj@mirnaGene[,mirnacol] %in% filteredmirnas,];

         # Test for, and remove, any groups with only one element
         # (otherwise "second argument must be a list" error)
         # singletGroups = tapply(groupPathways[,genecol], groupPathways[, pathwayidcol], length) == 1;
         
         # Generate a table storing mirna, gene, and mirnagene information for each Pathway ID
         pathsizes <- do.call(rbind, tapply( as.character(groupPathways[, genecol]) , groupPathways[, pathwayidcol], function(pathGenes)
         {
            # Make sure the pathway has one gene only once each
            pathGenes = unique(pathGenes);
            
            # If the following has 'unique' it collapses miRNA having multiple binding sites on the same gene to 1 entry
            univmirnaGene <- groupmirnaGene[groupmirnaGene[, genecol] %in% pathGenes,];
            if (length(univmirnaGene) > dim(groupmirnaGene)[2])
            {
               univmirnaGene <- matrix(ncol=dim(groupmirnaGene)[2], univmirnaGene);
            } else {
               univmirnaGene <- matrix(ncol=dim(groupmirnaGene)[2], byrow=TRUE,
                  univmirnaGene);
            }
            colnames(univmirnaGene) <- colnames(groupmirnaGene);
            univmirnagenecount <- dim(univmirnaGene)[1];

            enrichmirnaGene <- univmirnaGene[univmirnaGene[, mirnacol] %in% filteredmirnas,];
            if (length(enrichmirnaGene) > dim(groupmirnaGene)[2])
            {
               enrichmirnaGene <- matrix(ncol=dim(groupmirnaGene)[2], enrichmirnaGene);
            } else {
               enrichmirnaGene <- matrix(ncol=dim(groupmirnaGene)[2], byrow=TRUE,
                  enrichmirnaGene);
            }
            enrichedmirnagenecount <- dim(enrichmirnaGene)[1];
            numenrichmirna = 0;
            numenrichgenes = 0;
            if (enrichedmirnagenecount > 0)
            {
               colnames(enrichmirnaGene) <- colnames(univmirnaGene);
               rownames(enrichmirnaGene) <- 1:dim(enrichmirnaGene)[1];
               numenrichmirna = length(unique(enrichmirnaGene[, mirnacol]));
               numenrichgenes = length(unique(enrichmirnaGene[, genecol]));
            }

            r1 <- list("Pathway miRNA-Gene Count"=univmirnagenecount, "Enriched miRNA-Gene Count"=enrichedmirnagenecount,
               "Enriched miRNA-Genes"=enrichmirnaGene, "Genes Enriched"=numenrichgenes,
               "miRNAs Enriched"=numenrichmirna);
            return(r1);
         } ) )
         # End pathsizes tapply

         # For simplicity, remove pathways with no significant mirnaGenes
         pathsizes = pathsizes[sapply(pathsizes[,"Enriched miRNA-Gene Count"], function(i){i[[1]] > 0;}), ];

         # Target P-value command
         nFilteredmirnaGene = dim(filteredmirnaGene)[1]; # total number of target mirna-genes represented by the filtered miRNAs
         nmirnaGene = dim(groupmirnaGene)[1]; # total number of mirna-genes tested (filtered plus unfiltered miRNAs)
         useCts = do.call(c, pathsizes[, "Enriched miRNA-Gene Count"]); # total enriched mirnaGenes for each pathway
         gCounts = do.call(c, pathsizes[, "Pathway miRNA-Gene Count"]); # total mirna-genes in each pathway's universe
         pvs <- phyper(useCts - 1, nFilteredmirnaGene, nmirnaGene - nFilteredmirnaGene, gCounts, lower.tail = FALSE)
         names(pvs) = rownames(pathsizes);

         ord <- order(pvs)

         if (perm == 0)
         {
            r1 = list(pvalues = pvs[ord], "Measured pathway mirnaGenes" = gCounts[ord], 
               "Enriched pathway mirnaGenes" = useCts[ord], "Total mirnaGenes" = nmirnaGene, 
               "Total filtered mirnaGenes" = nFilteredmirnaGene,
               "Enriched miRNA-Genes" = pathsizes[ord, "Enriched miRNA-Genes"],
               "Genes Enriched"=sapply(pathsizes[ord, "Genes Enriched"],c),
               "miRNAs Enriched"=sapply(pathsizes[ord, "miRNAs Enriched"],c));
            r1;
         } else {
            r1 = list(pvalues = pvs);
            r1;
         }
      } )
      # End permresults lapply

      # Determine permutation P-value (if exists)
      if (permutations > 0)
      {
         permtable = do.call(rbind, lapply(1:length(permresults), function(j)
         {
            permresults[[j]]$pvalues;
         } ) );
         permps = apply(permtable, 2, function(i)
         {
            (rank(i)/(length(i)-1))[1];
         } )
         permresults[[1]]$permpvalues = permps;
      }
      permresults[[1]];
   } ); # End lapply(groups...)
   names(groupResults) = groups;

   mirnaobj@enrichment = groupResults;
   mirnaobj@state = "enriched";

   # Double-check that the pathwayList is present, just in case
   if (!"pathwayList" %in% slotNames(mirnaobj) | length(mirnaobj@pathwayList) == 0 )
   {
      pathwayTable = unique(mirnaobj@mirnaPathways[,c(pathwayidcol, pathwaycol)]);
      pathwayList = pathwayTable[,pathwaycol];
      names(pathwayList) = pathwayTable[,pathwayidcol];
      mirnaobj@pathwayList = pathwayList;
   }

   return(mirnaobj);
}

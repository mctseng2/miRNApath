`loadmirnapathways` <-
function
( mirnaobj,
  pathwayfile,
  genecol = "Entrez Gene ID",
  pathwaycol = "PATHWAY",
  columns = c(),
  pathwayidcol = NA)
{
   mirnaPathways <- read.table(pathwayfile, header=TRUE, sep="\t", fill=FALSE, comment.char = "", check.names = FALSE, quote="\""); #\"
   checkColumns( data = mirnaPathways, mandatory = c(pathwaycol, genecol) );

   mirnaGene <- mirnaobj@mirnaGene;

   if (is.na(pathwayidcol))
   {
      pathwayids = as.numeric(as.factor(mirnaPathways[,pathwaycol]));
      pathwayidcol = "PATHWAY_ID";
      mirnaPathways[,pathwayidcol] = pathwayids;
   }

   columns = c( columns, genecol = genecol );
   columns = columns[!is.na(columns)];
   for ( column in names(columns) )
   {
      if ( length(mirnaobj@columns[names(mirnaobj@columns) %in% column]) == 0 )
      {
         stop( paste( c("The following column specified in the loadmirnatogene method is not present in mirnaobj:", column, "\n"), collapse="\n") );
      }
      if ( length(colnames(mirnaGene)[colnames(mirnaGene) %in% mirnaobj@columns[column] ]) == 0 )
      {
         stop( paste( c("The following column specified in the loadmirnatogene method is not present in mirnaGene:", column, "\n"), collapse="\n") );
      }
      tempcol = mirnaobj@columns[column];
      colnames(mirnaPathways)[colnames(mirnaPathways) %in% columns[column]] = tempcol;
   }
   mirnaobj@columns["pathwaycol"] = pathwaycol;
   mirnaobj@columns["pathwayidcol"] = pathwayidcol;
   mirnaobj@mirnaPathways = mirnaPathways;
   pathwayTable = unique(mirnaPathways[,c(pathwayidcol, pathwaycol)]);
   pathwayList = pathwayTable[,pathwaycol];
   names(pathwayList) = pathwayTable[,pathwayidcol];
   mirnaobj@pathwayList = as.character(pathwayList);
   pathwaycount = length(levels(as.factor(mirnaPathways[,pathwaycol])));
   mirnaobj@pathwaycount = pathwaycount;
   return(mirnaobj);
}


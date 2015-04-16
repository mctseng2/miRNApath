`loadmirnatogene` <-
function
( mirnafile,
  mirnaobj,
  mirnacol = "miRNA Name",
  genecol = "Entrez Gene ID",
  columns = NA )
{
   mirnaGene <- read.table(mirnafile, header=TRUE, sep="\t", fill=FALSE, comment.char = "", check.names = FALSE, quote="\"");
   checkColumns( data = mirnaGene, mandatory = c(mirnacol, genecol) );

   columns = c( columns, mirnacol = mirnacol );
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
      colnames(mirnaGene)[colnames(mirnaGene) %in% columns[column]] = tempcol;
   }
   mirnacol = mirnaobj@columns["mirnacol"];
   #CompositeLabel = apply(mirnaobj@mirnaGene[, c(mirnacol, genecol)], 1, function(i)
   CompositeLabel = apply(mirnaGene[, c(mirnacol, genecol)], 1, function(i)
   {
      paste( c(i[1], " (", i[2], ")"), collapse="" );
   });
   mirnaGene[,"miRNA-Gene"] = CompositeLabel;
   mirnaobj@columns["mirnagene"] = "miRNA-Gene";
   mirnaobj@columns["genecol"] = genecol;
   mirnaobj@mirnaGene = mirnaGene;
   return(mirnaobj);
}


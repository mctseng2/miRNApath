setClass("mirnapath", 
   representation(
      mirnaTable="data.frame", 
      columns="character", 
      groupcount="numeric", 
      state="character", 
      mirnaGene="data.frame", 
      mirnaPathways="data.frame", 
      pathwaycount="numeric", 
      filters="numeric", 
      enrichment="list", 
      pathwayList="character"),
   prototype(
      mirnaTable=data.frame(),
      columns=character(0), 
      groupcount=numeric(0), 
      state=character(0), 
      mirnaGene=data.frame(), 
      mirnaPathways=data.frame(), 
      pathwaycount=numeric(0), 
      filters=numeric(0), 
      enrichment=list(), 
      pathwayList=character(0)),
   sealed=TRUE,
   package="miRNApath" )

`loadmirnapath` <-
function
( mirnafile = "mirna_input.txt",
  mirnacol = "miRNA Name",
  assayidcol = "ASSAYID",
  groupcol = "GROUP",
  filterflagcol = "FILTERFLAG",
  expressioncol = NA,
  foldchangecol = NA,
  pvaluecol = NA)
{
   mirnaTable <- read.table(mirnafile, header=TRUE, sep="\t", fill=FALSE, comment.char = "", check.names = FALSE, quote="\"");
   # Check columns, that they exist as specified
   checkColumns(data = mirnaTable, mandatory = c( mirnacol, groupcol ), numeric = NA );
   columns = c( "mirnacol" = mirnacol, "assayidcol" = assayidcol, "groupcol" = groupcol, "filterflagcol" = filterflagcol );
   if (!is.na(groupcol))
   {
      # Count the number of groups in the dataset
      groupcount = length(levels(as.factor(mirnaTable[,groupcol])));
   } else {
      groupcount = 0;
   }

   # check whether the data is already prepared for enrichment with hits and full library
   if ( any( colnames(mirnaTable) %in% filterflagcol ) )
   {
      state = "filtered";
   } else {
      state = "unfiltered";
   }

   #mirnaobj <- list(mirnaTable = mirnaTable, columns = columns, groupcount = groupcount, state = state );
   mirnaobj <- new("mirnapath", mirnaTable = mirnaTable, columns = columns, groupcount = groupcount, state = state );
   #class(mirnaobj) = "mirnapath";
   return(mirnaobj);
}


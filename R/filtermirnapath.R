`filtermirnapath` <-
function
( mirnaobj,
  pvalue = NA,
  expression = NA,
  foldchange = NA )
{
   filterflagvals = rep(2, dim(mirnaobj@mirnaTable)[1] );
   if ( length( mirnaobj@columns[names(mirnaobj@columns) %in% "filterflagcol"] ) > 0 )
   {
      filterflagcol = mirnaobj@columns["filterflagcol"];
   } else {
      filterflagcol = "FILTERFLAG";
      mirnaobj@columns = c(mirnaobj@columns, "filterflagcol" = filterflagcol);
   }
   if ( class(mirnaobj) != "mirnapath" )
   {
      stop("The mirnaobj object is not the proper class of mirnapath.");
   }
   filters = c();
   if ( !is.na(pvalue) )
   {
      if ( length( mirnaobj@columns[names(mirnaobj@columns) %in% "pvaluecol"] ) > 0)
      {
         pvaluecol = mirnaobj@columns[names(mirnaobj@columns) %in% "pvaluecol"];
         filterflagvals[mirnaobj@mirnaTable[,pvaluecol] <= pvalue & filterflagvals != 0 & !is.na(mirnaobj@mirnaTable[,pvaluecol]) ] = 1;
         filters = c(filters, "P-value" = pvalue);
      } else {
         stop("A P-value filtering value was supplied, but no pvaluecol exists in mirnaobj.");
      }
   }
   if ( !is.na(expression) )
   {
      if ( length( mirnaobj@columns[names(mirnaobj@columns) %in% "expressioncol"] ) > 0)
      {
         expressioncol = mirnaobj@columns[names(mirnaobj@columns) %in% "expressioncol"];
         filterflagvals[ (mirnaobj@mirnaTable[,expressioncol] >= expression | !is.na(mirnaobj@mirnaTable[,expressioncol]) ) & filterflagvals != 0 ] = 1;
         filters = c(filters, "Expression" = expression);
      } else {
         stop("An expression filtering value was supplied, but no expressioncol exists in mirnaobj.");
      }
   }
   if ( !is.na(foldchange) )
   {
      if ( length( mirnaobj@columns[names(mirnaobj@columns) %in% "foldchangecol"] ) > 0)
      {
         foldchangecol = mirnaobj@columns[names(mirnaobj@columns) %in% "foldchangecol"];
         # Convert fold change style into "directional fold change" (not fraction)
         foldchanges = convertFoldChange( mirnaobj@mirnaTable[,foldchangecol], style = "directional");
         filterflagvals[foldchanges >= expression & filterflagvals != 0 & !is.na(foldchanges) ] = 1;
         filters = c(filters, "Fold change" = foldchange);
      } else {
         stop("An expression filtering value was supplied, but no expressioncol exists in mirnaobj.");
      }
   }
   if (length(filters) == 0)
   {
      print("No filters were applied.");
      mirnaobj@status = "unfiltered";
   } else {
      mirnaobj@mirnaTable[,filterflagcol] = filterflagvals;
      mirnaobj@filters = filters;
      mirnaobj@state = "filtered";
   }
   return(mirnaobj);
}


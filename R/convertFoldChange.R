`convertFoldChange` <-
function
( foldchanges, style )
{
   if ( any(foldchanges < 0) )
   {
      if (any(foldchanges < 1 & foldchanges > 0))
      {
         obsstyle = "log directional";
         if (style == "directional")
         {
            foldchanges = 2^foldchanges;
            foldchanges[foldchanges < 1] = -1 / (foldchanges[foldchanges < 1]);
         } else if (style == "fractional")
         {
            foldchanges = 2^foldchanges;
         }
      } else {
         obsstyle = "directional";
         if (style == "log directional")
         {
            foldchanges[foldchanges < 0] = -1 / foldchanges[foldchanges < 0];
            foldchanges = log2(foldchanges);
         } else if (style == "fractional")
         {
            foldchanges[foldchanges < 0] = -1 / foldchanges[foldchanges < 0];
         }
      }
   } else {
      obsstyle = "fractional";
      if (style == "directional")
      {
         foldchanges[foldchanges < 1] = -1 / (foldchanges[foldchanges < 1]);
      } else if (style == "log directional")
      {
         foldchanges = log2(foldchanges);
      }
   }
   return(foldchanges);
}


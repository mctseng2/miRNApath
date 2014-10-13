`checkColumns` <-
function
( data,
  mandatory = NA,
  numeric = NA )
{
   if ( all(!is.na(mandatory)) & length(mandatory) > 0 )
   {
      mand = sapply( mandatory, function(i)
      {
         if ( any(colnames(data) %in% i) )
         {
            T;
         } else {
            F;
         }
      } )
   } else { mand = NA }
   if ( all(!is.na(numeric)) & length(numeric) > 0 )
   {
      nume = sapply( numeric, function(i)
      {
         if ( any(colnames(data) %in% i) & all(is.numeric(data[,i])) )
         {
            T;
         } else {
            F;
         }
      } )
   } else { nume = NA }
   msg = NA;
   if ( length(mand[!mand]) > 0 & all(!is.na(mand)) )
   {
      msg = paste( c("The following mandatory columns are not found in the data:", paste(c(names(mand)[!mand]), collapse=", ")), collapse=" ");
   }
   if ( length(nume[!nume]) > 0 & all(!is.na(nume)) )
   {
      msg2 = paste( c("The following columns are required to be numeric:", paste(c(names(nume)[!nume]), collapse=", ")), collapse=" ") 
      if ( is.na(msg) )
      {
         msg = msg2;
      } else {
         msg = paste( c(msg, msg2), collapse="\n" );
      }
   }
   if ( !is.na(msg) )
   {
      stop(msg);
   }
}


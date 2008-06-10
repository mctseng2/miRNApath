\name{convertFoldChange}
\alias{convertFoldChange}
\title{ Converts Fold Change to Consistent Format }
\description{
  The script tries to determine if the fold changes
  are in an expected format, and if not, does its best
  to convert correctly. Internal function.
}
\usage{
convertFoldChange(foldchanges, style)
}
\arguments{
  \item{foldchanges}{ numerical values with fold changes }
  \item{style}{ type of fold change to return. }
}
\details{
  Allowed formats are:
  "log directional" (contains continuous numbers positive and
  negative, including numbers between -1 and 1),
  "directional", includes positive and negative values, but no
  numbers between -1 and 1, and
  "fractional", includes values greater than zero, especially
  those less than one.
}
\value{
  numerical values in the defined format
}
\author{ James M. Ward \email{jmw86069@gmail.com} }
\examples{
## No examples currently.
}

\keyword{ manip }
\keyword{ internal }

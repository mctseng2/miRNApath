\name{checkColumns}
\alias{checkColumns}
\title{ Check data.frame for Required Columns }
\description{
  Internal function designed to assist checking for required
  columns in a data.frame.
}
\usage{
checkColumns(data, mandatory=NA, numeric=NA)
}
\arguments{
  \item{data}{ data.frame }
  \item{mandatory}{ list of required columns }
  \item{numeric}{ list of columns expected to be numeric }
}
\value{
  boolean, or stop() with columns not meeting the criteria
}
\author{ James M. Ward \email{jmw86069@gmail.com} }
\examples{
## No examples currently
}

\keyword{ manip }
\keyword{ internal }

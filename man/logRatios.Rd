\name{logRatios}
\alias{logRatios}
\alias{logRatios,cghSet-method}
\alias{cghSet-class}

\title{extractor for cghSet assay data}
\description{extractor for cghSet assay data}
\usage{
logRatios(cghSet)
}
\arguments{
  \item{cghSet}{ instance of \linkS4class{cghSet} class}
}
\details{
gets the AssayData element
}
\value{
matrix
}
%\references{  }
\author{Vince Carey <stvjc@channing.harvard.edu>}
%\note{  }


%\seealso{  }

\examples{
data(neveCGHmatch)
logRatios(neveCGHmatch)[1:4,1:4]
}
\keyword{ models }

\name{neveCGHmatch}
\alias{neveCGHmatch}
\alias{neveRMAmatch}
\alias{neveExCGH}
\docType{data}
\title{ Neve Cancer Cell 2006 expression plus CGH data }
\description{
ExpressionSet and cghSet
}
\usage{data(neveCGHmatch)
data(neveRMAmatch)
data(neveExCGH)
}
\format{
The individual datasets (aCGH and expression assays)
take the form of
a cghSet for neveCGHmatch and an ExpressionSet
for neveRMAmatch.  There are only 50 samples because
only 50 could be aligned on the given sample name
tokens in the caArrayDB data as of June 9 2007.
Those sample name tokens are very mangled in the CEL files.

The combined data structure \code{neveExCGH}
has a special container class \code{cghExSet}.

}
\source{
links are provided in the pdf of the Cancer Cell paper;
see the PMID of \code{experimentData(neveCGHmatch)}
}
\references{
PMID 17157791
}
\examples{
data(neveCGHmatch)
neveCGHmatch
logRatios(neveCGHmatch)[1:4,1:4]
data(neveRMAmatch)
neveRMAmatch
}
\keyword{datasets}

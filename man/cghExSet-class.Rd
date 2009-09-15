\name{cghExSet-class}
\docType{class}
\alias{cghExSet-class}
\alias{cloneMeta,cghExSet-method}
\alias{cloneNames,cghExSet-method}
\alias{exprs,cghExSet-method}
\alias{initialize,cghExSet-method}
\alias{logRatios,cghExSet-method}
\alias{show,cghExSet-method}
\alias{cloneMeta}
\alias{cloneNames}
\alias{make_cghExSet}

\title{Class "cghExSet" for combining CGH and expression data }
\description{	 combination of an ExpressionSet and CGH assay results }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("cghExSet",
  phenoData, experimentData, annotation, exprs, logRatios, cloneMeta)}. 
}
\usage{
make_cghExSet(exprs, logRatios, cloneMeta, pd, mi, anno) # pd is AnnotatedDataFrame, mi is MIAME
}
\arguments{
  \item{exprs}{ matrix of expression assay results }
  \item{logRatios}{ matrix of aCGH assay results }
  \item{cloneMeta}{AnnotatedDataFrame for aCGH clone descriptions}
  \item{pd}{ AnnotatedDataFrame for sample level data }
  \item{mi}{ \code{\link[Biobase:class.MIAME]{MIAME}} instance for experiment
    documentation } 
  \item{anno}{character string with annotation platform descriptor for
    expression data} 
}
\section{Slots}{
	 \describe{
    \item{\code{cghAssays}:}{Object of class \code{"AssayData"} rectangular
representation of logRatio data from CGH }
    \item{\code{cloneMeta}:}{Object of class \code{"AnnotatedDataFrame"} information
on chromosome and offset of clones}
    \item{\code{assayData}:}{Object of class \code{"AssayData"} expression assay results }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} sample level data }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} reporter
level metadata for expression assay results }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"} container for
experiment documentation }
    \item{\code{annotation}:}{Object of class \code{"character"} identifiers for
expression and CGH platforms, as a named vector with elements named 'exprs'
and 'logRatios' }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}  }
  }
}
\section{Extends}{
Class \code{\link[Biobase:class.eSet]{eSet-class}}, directly.
Class
\code{\link[Biobase:class.VersionedBiobase]{VersionedBiobase-class}}, by
class "eSet", distance 2. 
Class \code{\link[Biobase:class.Versioned]{Versioned-class}}, by class
"eSet", distance 3. 
}
\section{Methods}{
  \describe{
    \item{cloneMeta}{\code{signature(cghSet = "cghExSet")}: extract annotated data
frame on clone locations for CGH component }
    \item{cloneNames}{\code{signature(cghSet = "cghExSet")}: extract character vector
of clone IDs for CGH component }
    \item{exprs}{\code{signature(object = "cghExSet")}: extract expression assay results }
    \item{initialize}{\code{signature(.Object = "cghExSet")}: infrastructure }
    \item{logRatios}{\code{signature(cghSet = "cghExSet")}: extract CGH assay results }
    \item{show}{\code{signature(object = "cghExSet")}: display object in concise form }
	 }
}
\references{R. M. Neve Cancer Cell Dec 2006 }
\author{ V Carey <stvjc@channing.harvard.edu> }
%\note{ ~~further notes~~ }
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%
%\seealso{
%	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%	or \code{\linkS4class{CLASSNAME}} for links to other classes
%}
\examples{
showClass("cghExSet")
data(neveExCGH)
logRatios(neveExCGH)[1:4,]
exprs(neveExCGH)[1:4,]
}
\keyword{classes}

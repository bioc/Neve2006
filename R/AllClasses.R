setClass("cghSet", contains="eSet")

setClass("cghExSet", representation(
    cghAssays="AssayData", cloneMeta="AnnotatedDataFrame"), contains="eSet",
    prototype = prototype(cghAssays=assayDataNew(), cloneMeta=new("AnnotatedDataFrame")))

setGeneric("logRatios", function(cghSet) standardGeneric("logRatios"))
setMethod("exprs", "cghExSet", function(object) get("exprs",object@assayData))


setMethod("logRatios", "cghSet", function(cghSet)
   assayDataElement(cghSet, "exprs"))

setMethod("initialize", "cghExSet",
          function(.Object,
                   phenoData = new("AnnotatedDataFrame"),
                   experimentData = new("MIAME"),
                   annotation = character(),
                   exprs = new("matrix"),
                   logRatios = new("matrix"),
		   cloneMeta = new("AnnotatedDataFrame")) {
            .Object = callNextMethod(.Object,
                           assayData = assayDataNew(
                             exprs=exprs),
                           phenoData = phenoData,
                           experimentData = experimentData,
                           annotation = annotation)
            .Object@cghAssays = assayDataNew(logRatios =logRatios)
            .Object@cloneMeta = cloneMeta
            .Object
          })

setMethod("logRatios", "cghExSet", function(cghSet)
   assayDataElement(cghSet@cghAssays, "logRatios"))
setGeneric("cloneNames", function(cghSet) standardGeneric("cloneNames"))
setMethod("cloneNames", "cghExSet", function(cghSet) rownames(logRatios(cghSet)))
setGeneric("cloneMeta", function(cghSet) standardGeneric("cloneMeta"))
setMethod("cloneMeta", "cghExSet", function(cghSet) cghSet@cloneMeta)

make_cghExSet = function(exprs, logRatios, cloneMeta, pd, mi, anno) {
    if (!is(exprs, "matrix"))
        stop("exprs must be of class matrix")
    if (!is(logRatios, "matrix"))
        stop("racs must be of class matrix")
    if (!is(pd, "phenoData") & !is(pd, "AnnotatedDataFrame"))
        stop("pd must be of class phenoData or AnnotatedDataFrame")
    new("cghExSet", exprs=exprs, logRatios=logRatios, cloneMeta=cloneMeta, 
        phenoData = pd, experimentData = mi,
        annotation = anno)
}

setMethod("show", "cghExSet", function(object) {
    cat("cghExSet instance (CGH logRatios + expression)\n")
    cat("CGH assayData:\n")
  cat("  Storage mode:", storageMode(object@cghAssays), "\n")
  nms <- selectSome(ccc <- cloneNames(object))
  cat("  featureNames:", paste(nms, collapse=", "))
  if ((len <- length(ccc)) > length(nms))
    cat(" (", len, " total)", sep="")
  cat("\n  Dimensions:\n")
  print(Biobase:::assayDataDims(object@cghAssays))
  cat("\nexpression assayData\n")
  cat("  Storage mode:", storageMode(object), "\n")
  nms <- selectSome(featureNames(object))
  cat("  featureNames:", paste(nms, collapse=", "))
  if ((len <- length(featureNames(object))) > length(nms))
    cat(" (", len, " total)", sep="")
  cat("\n  Dimensions:\n")
  print(dims(object))
  cat("\nphenoData\n")
  show(phenoData(object))
  cat("\n")
  show(experimentData(object))
  cat("\nAnnotation\n")
  show(annotation(object))
    })

setMethod("[", "cghExSet", function(x, i, j, ..., drop=FALSE ) {
    if (missing(drop)) drop <- FALSE
    if (missing(i) && missing(j)) return(callNextMethod())
    if (!missing(j)) {
      nass = assayDataNew("lockedEnvironment", logRatios= x@cghAssays$logRatios[,j,drop=drop])
      x@cghAssays = nass
      x@assayData = assayDataNew("lockedEnvironment", exprs=exprs(x)[,j])
      x@phenoData = phenoData(x)[j,]
      }
    if (!missing(i)) {
      warning("when using first subscript, subsetting expression features only")
      x@assayData = assayDataNew("lockedEnvironment", exprs=exprs(x)[i,])
    }
   x
})
    

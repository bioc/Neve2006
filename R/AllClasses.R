setClass("cghSet", contains="eSet")

setGeneric("logRatios", function(cghSet) standardGeneric("logRatios"))

setMethod("logRatios", "cghSet", function(cghSet)
   assayDataElement(cghSet, "exprs"))

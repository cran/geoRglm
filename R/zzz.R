".onAttach" <- function(lib, pkg)
{
  messages.screen <- ifelse(is.null(getOption("geoR.messages")), TRUE, getOption("geoR.messages"))
  packageStartupMessage("---------------------------------------------------------\n")
  pkg.info <- drop(read.dcf(file=system.file("DESCRIPTION", package="geoRglm"), fields=c("Title","Version","Date")))
  packageStartupMessage(pkg.info["Title"])
  packageStartupMessage("\n")
  packageStartupMessage(paste("geoRglm version ", pkg.info["Version"], " (", pkg.info["Date"], ") is now loaded\n", sep=""))
  packageStartupMessage("-----------------------------------------------------------\n")
  packageStartupMessage("\n")
  return(invisible(0))
}


".bilinearformXAY"<-geoR:::.bilinearformXAY 
".check.locations"<-geoR:::.check.locations
".cond.sim"<-geoR:::.cond.sim
".cor.number"<-geoR:::.cor.number
".diagquadraticformXAX"<-geoR:::.diagquadraticformXAX
## ".geoR_inout"<-geoR:::.geoR_inout
".ldots.set"<-geoR:::.ldots.set
".prepare.graph.krige.bayes"<-geoR:::.prepare.graph.krige.bayes
".solve.geoR"<-geoR:::.solve.geoR
"hist.krige.bayes"<-geoR:::hist.krige.bayes
## "image.kriging"<-geoR:::image.kriging

".First.lib" <- function(lib, pkg)
{
  messages.screen <- ifelse(is.null(getOption("geoR.messages")), TRUE, getOption("geoR.messages"))
  if(messages.screen){
    cat("-----------------------------------------------------------\n")
    if(!require(geoR)){
      cat("\n")
      cat("Package geoR is required by geoRglm\n")
      cat("It should also be installed and loaded\n")
    }
  }
  library.dynam("geoRglm", package=pkg, lib.loc=lib)
  if(messages.screen){
    pkg.info <- packageDescription("geoRglm", lib.loc = lib, fields=c("Title","Version","Date"))
    cat(pkg.info$Title)
    cat("\n")
    cat(paste("geoRglm version ", pkg.info$Version, " (", pkg.info$Date, ") is now loaded\n", sep=""))
    cat("-----------------------------------------------------------\n")
    cat("\n")
  }
  return(invisible(0))
}

".First.lib" <-
  function(lib, pkg)
{
  cat("-----------------------------------------------------------\n")
  if(!require(geoR)){
    cat("\n")
    cat("Package geoR is required by geoRglm\n")
    cat("It should also be installed and loaded\n")
  }
  library.dynam("geoRglm", package=pkg, lib.loc=lib)
  cat(package.description("geoRglm", lib = lib, field="Title"))
  cat("\n")
  ver <- package.description("geoRglm", lib = lib, field="Version")
  dat <- package.description("geoRglm", lib = lib, field="Date")
  cat(paste("geoRglm version ", ver, " (", dat, ") is now loaded\n", sep=""))
  cat("-----------------------------------------------------------\n")
  cat("\n")
  return(invisible(0))
}

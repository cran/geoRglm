".First.lib" <-
  function(lib, pkg)
{
  cat("-----------------------------------------------------------\n")
  if(is.R() && !require(geoR)){
    cat("\n")
    cat("Package geoR is required by geoRglm\n")
    cat("It should also be installed and loaded\n")
  }
  if(is.R()){
    library.dynam("geoRglm", package=pkg, lib.loc=lib)
    cat(package.description("geoRglm", lib = lib, field="Title"))
    cat("\n")
    ver <- package.description("geoRglm", lib = lib, field="Version")
    cat(paste("geoRglm version", ver,  "is now loaded\n"))
    #dumm <- package.description("geoR", lib = lib, field="Version")
    #cat(paste("geoR test", dumm,  "\n"))
  }
  else{
    cat("geoSglm: a package for generalised linear spatial models \n")
    cat("geoSglm is loaded\n")
  }
  cat("-----------------------------------------------------------\n")
  cat("\n")
  return(invisible(0))
}

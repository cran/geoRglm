"geoRglmdefunct" <- function()
{
  cat("\n")
  cat("The following functions are no longer used in geoRglm:")
  cat("---------------------------------------------------")
  cat(" pois.log.krige: use pois.krige instead")
  cat(" y50: use p50 instead")
  cat("\n")
}

cite.geoRglm <- function()
{
  cat("\n")
  cat("To cite geoR in publications, use\n\n")
  msg <- "CHRISTENSEN, O.F. & RIBEIRO Jr., P.J. (2002) geoRglm: A package for generalised linear spatial models. R-NEWS, Vol 2, No 2, 26-28. ISSN 1609-3631."
  writeLines(strwrap(msg, prefix = "  "))
  cat("\n")
  msg <- paste("Please cite geoRglm when using it for data analysis!")
  writeLines(strwrap(msg))
  cat("\nA BibTeX entry for LaTeX users is\n\n")
  cat("  @Article{,\n")
  cat("     title	   = {Christensen, O.F. and Ribeiro Jr., P.J.},\n")
  cat("     author        = {{geoRglm}: A package for generalised linear spatial models},\n")
  cat("     journal       = {R-NEWS},\n")
  cat("     year	   = {2002},\n")
  cat("     volume	   = {2},\n")
  cat("     number	   = {2},\n")
  cat("     pages	   = {26--28},\n")
  cat("     issn          = {1609-3631},\n")
  cat("     url           = {http://cran.R-project.org/doc/Rnews}\n")
  cat("   }\n\n")
}

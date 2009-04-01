.First.lib <- function (lib, pkg){
#.initAdegenetClasses()
#.initAdegenetUtils()
  library.dynam("adegenet", pkg, lib)
  cat("   ==========================\n")
  cat("    adegenet 1.2-3 is loaded  \n")
  cat("   ==========================\n\n")

  cat(" - to start, type '?adegenet'\n")
  cat(" - to browse adegenet website, type 'adegenetWeb()'\n")
  cat(" - to post questions/comments: adegenet-forum@lists.r-forge.r-project.org\n\n")
}

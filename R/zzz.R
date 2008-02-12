.First.lib <- function (lib, pkg){
#.initAdegenetClasses()
#.initAdegenetUtils()
  library.dynam("adegenet", pkg, lib)
  cat("   ##########################\n")
  cat("   ### adegenet is loaded ### \n")
  cat("   ##########################\n\n")
  
  cat(" - to start, type '?adegenet'\n")
  cat(" - to browse adegenet website, type 'adegenetWeb()'\n\n")
}

# Script checks target R-packages from 'temp_CRAN_packages.txt' and 'temp_BiocManager_packages.txt' among intalled R-packages

check_packages <- function(package_list, source='cran') {
  list.of.packages <- package_list
  if (source == 'cran') {
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  }
  if (source == 'BiocManager') {
    new.packages <- list.of.packages[!(list.of.packages %in% row.names(installed.packages()))]
  }
  return(new.packages)
}


#get args from CMD
sysArgs <- commandArgs(trailingOnly = TRUE) 
CRAN <- as.integer(sysArgs[1])  #1 or 0
BiocManager <- as.integer(sysArgs[2])  #1 or 0
 
new.packages = c()
if (CRAN) {

    CRAN_packages <- read.table(file = "temp_CRAN_packages.txt", header = FALSE, sep = ',')
    new.packages <- c(new.packages, check_packages(CRAN_packages))

}
if (BiocManager) {

    BiocManager_packages <- read.table(file = "temp_BiocManager_packages.txt", header = FALSE, sep = ',')
    new.packages <- c(new.packages, check_packages(BiocManager_packages, source='BiocManager'))

}

if (length(new.packages)) {

    print('Packages not found:')    
    print(unlist(new.packages, use.names = FALSE))
    print('They will be installed. This may take a long time (if you use google collab, then up to 20 minutes, usually 2-7 minutes)')
}  

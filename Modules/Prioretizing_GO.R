# Script for prioretizing GO-terms
# see: GOxploreR doi:10.1038/s41598-020-73326-3
# work with R.4-3-3
# inputs:
#	input_name: path to csv file with filename. Ex: "C:/User/data.csv"
#	organism: name of organism. Ex: "Human"	
#	domain: name of domain in GO-graph. Available inputs: "BP" - Biological Process
#														  "CC" - Cellular Component
#														  "MF" - Molecular Functions
# output: csv-file contains prioritized GO-terms list 	

# Check packages and install
check_packages <- function(package_list, source='cran') {
  list.of.packages <- package_list
  if (source == 'cran') {
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    if(length(new.packages)) install.packages(new.packages, repos=c("http://rstudio.org/_packages",
                                                                    "http://cran.rstudio.com"),
                                              dependencies=TRUE)
  }
  if (source == 'BiocManager') {
    new.packages <- list.of.packages[!(list.of.packages %in% row.names(installed.packages()))]
    if(length(new.packages)) {
      chooseBioCmirrorror(graphics=FALSE, id=1)
      BiocManager::install(new.packages)
    }
  }
}

check_packages(c("data.table", "BiocManager", "utils", 	"ggplot2"))
check_packages(c("GO.db", "annotate", "biomaRt"), source='BiocManager')
check_packages(c("GOxploreR"))

# import packages
suppressMessages({
  library(data.table)
  library(GO.db)
  library(GOxploreR)
})

#get args from CMD
sysArgs <- commandArgs(trailingOnly = TRUE)
input_name <- sysArgs[1] #filename with absolute path
organism <- sysArgs[2]   #organism, Example: "Human"
domain <- sysArgs[3]     #domain of GO-graph: "BP", "CL" or "MF"

output_name <- paste(substr(input_name, 1, nchar(input_name)-4), '_output.csv', sep='')
#get terms from csv
terms <- read.csv(input_name, header = TRUE)[[1]]

prioretizingGO <- function(terms, organism, domain) {
  ordered_bioproc_go <- distRankingGO(terms)
  prior_bp_terms <- prioritizedGOTerms(lst = ordered_bioproc_go$`GO-terms_ranking`,
                                       organism = organism, 
                                       sp = TRUE, 
                                       domain = domain)
  return(prior_bp_terms$HF) # список приорет терминов
}

prior_terms <- as.data.table(prioretizingGO(terms, organism=organism, domain=domain))
colnames(prior_terms) <- c('Term')
write.csv(prior_terms, output_name, row.names = FALSE)
print('file sucessfully saved in:')
print(output_name)
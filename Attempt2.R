#
#
# more learning
#setup working directory pathway
working.dir <- getwd()

# create all folders for this project so that it is organized & reproducible 
# store future file names in an object... These are base folder names
output.folder.names <- c("Clean Data", "Figures",
                         "Raw Data","Analysis","Results")

# make all folders for this project
# and make the folders if they don't exit yet.
for(i in 1:length(output.folder.names))
  if(file.exists(output.folder.names[i]) == FALSE)
    dir.create(output.folder.names[i])
#these are the pathways which are necessary to send graphs to the folders.
path.clean <- paste(working.dir, "/", output.folder.names[1], "/", sep = "")
path.figures <- paste(working.dir, "/", output.folder.names[2], "/", sep = "")
path.raw.data <- paste(working.dir, "/", output.folder.names[3], "/", sep = "")
path.analysis <- paste(working.dir, "/", output.folder.names[4], "/", sep = "")
path.results <- paste(working.dir, "/", output.folder.names[5], "/", sep = "")

# to install the core packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

# Run the first time to install:
source("https://bioconductor.org/biocLite.R")
biocLite("ComplexHeatmap")
# need the library complexheatmap
library(ComplexHeatmap)
# yay instalations are complete!


# read in practice data set 
# this is my second attempt at understanding and applying the basic 
# prinicals of creating heatmaps in R





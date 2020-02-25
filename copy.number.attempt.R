#
#
BiocInstaller::biocLite("ComplexHeatmap")



library("BiocManager")
Bioconductor version 3.9 (BiocManager 1.30.4), ?BiocManager::install

library(ComplexHeatmap)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
source("https://bioconductor.org/biocLite.R")
BiocManager::install(c(
  "GenomicFeatures", "AnnotationDbi"), update = TRUE, ask = FALSE)

# frmo biocLite
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install(version = "3.90")
library(BiocManager)
# import data set from lesson 8 of practice data example 


# import copy number data
copy_number_data <- read.csv("~/Downloads/Data and scripts for all videos/Lesson-08/copy_number_data.csv", header=FALSE)
# name it somethinga  bit easier to work with!
cn_data <- copy_number_data
# now let's view the structure of our data
View(cn_data)
head(cn_data)
# now let's check the dimensions of the data set
dim(cn_data) 
# 54, and 103
nrow(cn_data) 
# the number of rows in the data set 54
# in this case nrow is the locations (aka bins) in the genome
ncol(cn_data) 
# the number of columns in the data set 103
# in this case that is the number of cells


# now let's check the class fo the data
# this is important becasue in order to create a heat map we need to have a
# data matrix
# use the function class which will tell us the type of object that is specified
?class
class(cn_data) 
# we have a data frame 

# Make the data into a matrix so that it can be compatible with the heatmap
# function
# however, there are 3 columns that we do not want to have in the matrix
# we are going to leave out the columns (chrom,start,end) since 
# they aren't relevant in the heatmap 
# because these are the first 3 columns in the data set we can select
# for all of the rows, and only columns 4

m_cn_data <- as.matrix(cn_data[  ,c(4:100)]) # [all rows, columns 4-100]
# leave out the first 3 columns (chrom,start,end) since they don't belong in the heatmap itself

# check to see the object class
class(m_cn_data)

head(my_matrix) # woohoo it's a matrix! 
# let's check and make sure the correct # of columns are present
head(m_cn_data)
# awesome, the first 3 columns are not present in this matrix
# however, for annotating the heatmap later, we can save the chrosome column
# save it as a data frame here, but I am not sure if it would still work 
# as just an object??

class(m_cn_data)

# Save chromosome column (CHR)
# save as a fata frame
chromosome_anno <- data.frame(chrom = cn_data$V1)
chromosome_anno
# yay now we have a data frame with only one column


Heatmap(m_cn_data)

heatmap(m_cn_data)
# okay so seeing as the package won't work let's use the shitty default setting
heatmap(cn_data)
# nope okay now nothing is working
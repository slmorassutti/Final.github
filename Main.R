#
#
# Main directory for Final Project
# Created by Sam Morassutti
# Feb 2020
# R.version.string
#"R version 3.6.2 (2019-12-12)"

# libraries needed (a lot)
source("https://bioconductor.org/biocLite.R")
biocLite("ComplexHeatmap")
library(ComplexHeatmap)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install(version = "3.10")
library(BiocManager)

library(dplyr)
# libraries needed
library(devtools)
# and the package complex heatmap 
install_github("jokergoo/ComplexHeatmap")
library(circlize)
library(cluster)


# index 
#
#
#
#
#


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



# This file is the first step in my final project
# It is intended for me to familiorize myself with how to construct a heatmap
# Created by Sam Morassutti
#
#
#



#----------------------------------------------
#----- 1 Generate a random matrix---------------
#----------------------------------------------
# this random matrix will have three groups by columns
# and three groups by rows

# use the fucntion set.seed to generate random numbers
?set.seed

set.seed(123)
nr1 = 4; nr2 = 8; nr3 = 6; nr = nr1 + nr2 + nr3
nc1 = 6; nc2 = 8; nc3 = 10; nc = nc1 + nc2 + nc3 

# now let's combine all of the rows and columns in a matrix 
# let's call the matrix object mat!
mat = cbind(rbind(matrix(rnorm(nr1 * nc1, mean = 1, sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2 * nc1, mean = 0, sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3 * nc1, mean = 0, sd = 0.5), nr = nr3)),
            rbind(matrix(rnorm(nr1 * nc2, mean = 0, sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2 * nc2, mean = 1, sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3 * nc2, mean = 0, sd = 0.5), nr = nr3)),
            rbind(matrix(rnorm(nr1 * nc3, mean = 0.5, sd = 0.5), nr = nr1),
                  matrix(rnorm(nr2 * nc3, mean = 0.5, sd = 0.5), nr = nr2),
                  matrix(rnorm(nr3 * nc3, mean = 1, sd = 0.5), nr = nr3)))
# Now let's shuffle the rows and columns randomly 
mat = mat[sample(nr, nr), sample(nc, nc)] 
rownames(mat) = paste0("row", seq_len(nr))
colnames(mat) = paste0("column", seq_len(nc))
str(mat)
# Great! Now we have a random data matrix called mat
# that has 24 columns and 18 rows
mat
# woah those are some interesting numbers
# let's check the range by using the min and max function
min(mat) # -1.154584
max(mat) # 2.62052
# Now let's visualize this random data and create a heatmap with default settings
# yayayayayayay
# to get the lendeng and legend name
# use function name = "", then put the name of legend wooooo
ComplexHeatmap::Heatmap(mat, name = "mat")
# save as pdf and send to figures folder!
pdf(file = paste(path.figures,"Default Random Heatmap.pdf", sep="/"))
ComplexHeatmap::Heatmap(mat, name = "mat")
dev.off()



#---------- colours----------
# colours are a veyr important compoent of representation of data in the matrix
# colour mapping functions
#allow you to create a fucntion in which a vector of values will result
# in the correspond colours

# If using heatmap(), always use the circlcolorRam2() in order to make a color
# mapping function
# there are two arguements in the colorRamp2
# a vector of break values and their corresponding colors 


# create a color scheme that corresponds to the set values with a rang of -2 to 2
# in this map, the colors aren't affected by outliers
# the color mapping function is robost to outliers

col_fun = colorRamp2(c(-2, 0, 2), c("mediumvioletred", "white", "mediumspringgreen"))
col.fun <- col_fun(seq(-3, 3))
ComplexHeatmap::Heatmap(mat, name = "mat", col = col.fun)
# the negative values are mediumvioletred and the positive are mediumspringgreen

# save as pdf and send to figures folder!
pdf(file = paste(path.figures,"Random Heatmap2.pdf", sep="/"))
ComplexHeatmap::Heatmap(mat, name = "mat", col = col.fun)
dev.off()

# can also create a rainbow color scheme 
ComplexHeatmap::Heatmap(mat, name = "mat", col = rev(rainbow(10)))
# save as pdf and send to figures folder!
pdf(file = paste(path.figures,"Random Heatmap3.pdf", sep="/"))
ComplexHeatmap::Heatmap(mat, name = "mat", col = rev(rainbow(10)))
dev.off()

# if matrix is continuous, provide a vector of colours 
# the colours will be interpolated linearly
# however, mapping isn't robust against outliers
# becasue it starts with the min and max values
# in the matrix

# therefore, if you set the colour to the max/min
# of the matrix, it will be identical to the
# default plot


# can also have colours for a discrete numeric matrix
# and also for character matrix
# default clustering will be applied in the dsicrete plot
# and can be in character plots too if set 



# -----NA---------------
# NA is able to be in matrix
# but you can set the colour so it is identifiable
# the default colour is grey (how do I put this in legend)
# to change the colour use argument
# na_col


# let's try by making the na values black
ComplexHeatmap::Heatmap(mat, name = "mat", na_col = "black")
# it appears that there are no na values.....
# the plot is back to default colour scheme

#---------- borderds and spacing-------------
# create border by using the border true fucntion
ComplexHeatmap::Heatmap(mat, name = "mat", border = TRUE)
pdf(file = paste(path.figures,"Random Heatmap4 border.pdf", sep="/"))
ComplexHeatmap::Heatmap(mat, name = "mat", border = TRUE)
dev.off()


# add spacing between each value

#-------title-----------------
# describe what the plot is 
# can be set by column and or row
# add titles in the complexheatmap package
# so the funciton column_titles should work for the top.... 
# and the function 
# graphic paramaters can be set by the function gpar
?gpar()
# it is used by get.gpar()

# use the command row_title_gp, or column_title_gp
ComplexHeatmap::Heatmap(mat, name = "mat", col = col.fun, row_title = "Title of rows", 
        column_title = "Title of columns")
# save as pdf and send to figures folder
pdf(file = paste(path.figures,"Random Heatmap titles.1.pdf", sep="/"))
ComplexHeatmap::Heatmap(mat, name = "mat", col = col.fun, row_title = "Title of rows", 
                        column_title = "Title of columns")
dev.off()


# now try change the font and title the entire heatmap title
ComplexHeatmap::Heatmap(mat, name = "mat", column_title = "title", 
                        column_title_gp = gpar(fontsize = 15, fontface = "bold"))
# jk gpar is mean
# mmmm okay I don't know why this isn't working 

# title angles
# in this package you can also change the angle of the titles
# very useful! 



#----clustering-------------
# maybe the most important aspect/function of heatmap visulaization
# there are many clustering methods
# two default: euclidean, and pearson both with pre-determined distance methods 
# hclust, dengrodgram (object or class) that can be/ or already are clustered


# can also turn clustering off, do for rows and then columns
# depending on if clustering is needed or not
ComplexHeatmap::Heatmap(mat, name = "mat", cluster_rows = FALSE)
# save and send to figures!
pdf(file = paste(path.figures,"Random Heatmap cluster.1.pdf", sep="/"))
ComplexHeatmap::Heatmap(mat, name = "mat", cluster_rows = FALSE)
dev.off()

ComplexHeatmap::Heatmap(mat, name = "mat", cluster_columns = FALSE)
# so odd looking!
# save and send to figures
pdf(file = paste(path.figures,"Random Heatmap cluster.2.pdf", sep="/"))
ComplexHeatmap::Heatmap(mat, name = "mat", cluster_columns = FALSE)
dev.off()


# to show clustering on specififcc sides
# for example, show clustering on the right side of rows, and 
# the bottom of the columns
ComplexHeatmap:: Heatmap(mat, name = "mat", row_dend_side = "right", 
                         column_dend_side = "bottom")
# save figure and send to folder!
pdf(file = paste(path.figures,"Random Heatmap cluster.3.pdf", sep="/"))
ComplexHeatmap:: Heatmap(mat, name = "mat", row_dend_side = "right", 
                         column_dend_side = "bottom")
dev.off()



#--------------------------------------------------------------
#---------------Apply to real data set------------------------
#-------------------------------------------------------------


#------------------------------------------------------------
#--------import, save, and process---------------------------------
#----------------------------------------------------------------
# expression data set from OncoPrint
# read in data for lung adenocarcinoma
# this data set is in complexheatmap so we can read it in directory from here

mat = read.table(system.file("extdata", package = "ComplexHeatmap", 
                             "tcga_lung_adenocarcinoma_provisional_ras_raf_mek_jnk_signalling.txt"), 
                 header = TRUE, stringsAsFactors = FALSE, sep = "\t")
# save as a csv and send to raw data folder
write.csv(mat,paste(path.clean,"lung_adenocarcinoma_ras_raf_mek_jnk_signalling.csv",sep="/"))
# can also import the now saved csv to see the raw data 
View(lung_adenocarcinoma_ras_raf_mek_jnk_signalling)
# we can see just from the csv that there are 26 genes
# can notice the different types of mutations MUT, AMP, HOMDEL, MUT;AMP, NA

#------look at structure-----
head(mat)
dim(mat) # for dimensions
nrow(mat) # for number of rows 
ncol(mat) # for number of col

# check if it is matrix
class(mat)


# now, just to see what happens, let's create a plot with the imported data
# with default complexheatmap settings... only specify legend name "mat"
ComplexHeatmap:: Heatmap(mat, name = "mat")
# woah it is different this time... that is kinda awesome
# the default plot show the case ID on the left side of the plot with colours 
# coressponding to the samples in the huge legend that doesn;t fit on the page 
head(mat)
# well, that is a bit ugly... what is going on with the labels 
# let's save it for now though becasue it's so odd
pdf(file = paste(path.figures,"lung.oddness.pdf", sep="/"))
ComplexHeatmap:: Heatmap(mat, name = "mat")
dev.off()


# look at structure
str(mat)
# also a lot of NA


# THIS SECTION??????
# was described as process and examining the row and column
# names however, I am still unsure of what EXACTLY is happening and why here.....
mat[is.na(mat)] = ""
# look at row names and process the data set 
# Not completely sure why/why they are doing this section....
# this is trasnformaing the data frame into a matrix in order for it be 
# be compatible with the plotting functions

rownames(mat) = mat[, 1]
mat = mat[, -1]
mat=  mat[, -ncol(mat)]
mat = t(as.matrix(mat))
mat[1:3, 1:3]
# this removed a lot of values... but I am still not quite sure if it 
# removed all of the NA values 
# let's check the structure again
head(mat)
# now let's plot it again and see if it works 
# hopefully the error messages/warning will not show up now that it is a matix
ComplexHeatmap:: Heatmap(mat, name = "mat")
# interesting
# so now the genes are order on the right of the plot
# I think the samples on on the bottom 
# this is quite different than the first plot before the data was processed 
#----------------------------------------------------------------------
#--------mutation as colours----------------------------------------
#----------------------------------------------------------------

# let's move on to: identifying mutations with colours
# now that we know there are different mutation types/alterations
# the complex heatmap example only focused on the MUT, AMP, and HOMDEL


# let's give each one a colour so that it can be identified 
col = c("HOMDEL" = "blue", "AMP" = "red", "MUT" = "#008000")
alter_fun = list(
  background = function(x, y, w, h) {
    grid::grid.rect()(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
  HOMDEL = function(x, y, w, h) {
    grid::grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["HOMDEL"], col = NA))
  },
  # bug red
  AMP = function(x, y, w, h) {
    grid::grid.rect()(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["AMP"], col = NA))
  },
  # small green
  MUT = function(x, y, w, h) {
    grid::grid.rect()(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["MUT"], col = NA))
  }
)
# check to see if alterations are set to correct color 
col
# okay it seems to be all good yay
# can it be translated into a plot though????
# create a title for the plot so that it isn't so crazy

col = c("MUT" = "#008000", "AMP" = "red", "HOMDEL" = "blue")

oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "OncoPrint for TCGA Lung Adenocarcinoma, genes in Ras Raf MEK JNK signalling",
          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "HOMDEL", "MUT"), 
                                      labels = c("Amplification", "Deep deletion", "Mutation")))
# want to remove NA values from the matrix 
oncoPrint(mat,
          alter_fun = alter_fun, col = col, 
          remove_empty_columns = TRUE, remove_empty_rows = TRUE,
          column_title = column_title, heatmap_legend_param = heatmap_legend_param)


---------------
#-----------plotting mutations as colours???-------------
#-------------------------------------------------------

# from here on out I am pretty stuck.....
# all of the analysis that doesn't work, or is confusing is compliled 
# with annotations in the potential anyalsis folder 




























#
# Troubled Analysis
# This file was created to contain all of the parts of my Final Project code that
# either DO NOT WORK, or are confusing 
# Created by Sam Morassutti
# Feb 2020
# R.version.string
#"R version 3.6.2 (2019-12-12)"

#Index
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


#----------------------------------------------------
#---------- colours--------------------------------
#------------------------------------------------------
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


#--------------------------------------------------------------
# -----Show NA???-------------------------------------------------------
#------------------------------------------------------------
# NA is able to be in matrix
# but you can set the colour so it is identifiable
# the default colour is grey
# to change the colour use argument
# na_col


# let's try by making the na values black
ComplexHeatmap::Heatmap(mat, name = "mat", na_col = "black")
# it appears that there are no na values.....
# or this just didn't work....
# the plot is back to default colour scheme



#----------------------------------------------------------
#-------title styles???------------------------------------------------
#---------------------------------------------------------
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

# PROBLEM WITH TITLE STYLEEEEEEE
# now try change the font and title the entire heatmap title
ComplexHeatmap::Heatmap(mat, name = "mat", column_title = "title", 
                        column_title_gp = gpar(fontsize = 15, fontface = "bold"))
# jk gpar is mean
# mmmm okay I don't know why this isn't working 


# COME BACK to title angles when you can actually format the font size ect
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



---------------------------------------
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


# now, just to see what happens, let's create a plot with the imported data
# with default complexheatmap settings... only specify legend name "mat"
ComplexHeatmap:: Heatmap(mat, name = "mat")
# well, that is a bit ugly... what is going on with the labels 
# let's save it for now though becasue it's so odd
pdf(file = paste(path.figures,"lun.oddness.pdf", sep="/"))
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
rownames(mat) = mat[, 1]
mat = mat[, -1]
mat=  mat[, -ncol(mat)]
mat = t(as.matrix(mat))
mat[1:3, 1:3]
# this removed a lot of values... but I am still not quite sure if it 
# removed all of the NA values 
# let's check the structure again
str(mat)


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

#--------------------------------------------------------
#-----------plotting mutations as colours???-------------
#-------------------------------------------------------
# create title for plot 
column_title = ("Mutated genes in the Ras Raf MEK JNK signalling the Lung Adenocarcinoma")
heatmap_legend_param = list(title = "Alternations", at = c("HOMDEL", "AMP", "MUT"), 
                            labels = c("Deep deletion", "Amplification", "Mutation"))
ComplexHeatmap::oncoPrint(mat, alter_fun = alter_fun, col = col, column_title = column_title, 
                          heatmap_legend_param = heatmap_legend_param)

# this doesn;t seem to want to work...
# okay let's try without the package name in front 

oncoPrint(mat, alter_fun = alter_fun, col = col, column_title = column_title, 
          heatmap_legend_param = heatmap_legend_param)
# mmmmm
# i am confused about why this isn't working 
#-----------why-----you no work----------????------------
?grid::grid.rect()

# ------remove emplty rows-------------
# the default setting is that is samples/genes do not have mutations
# so if they have an NA, they will remain in the data set
# sometimes this is not useful so it is better to remove them
# even though the above visualiation didn't work, let's move a head
# and try and remove some NA's even though I am pretty sure I already did ....


ComplexHeatmap::oncoPrint(mat,
                          alter_fun = alter_fun, col = col, 
                          remove_empty_columns = TRUE, remove_empty_rows = TRUE,
                          column_title = column_title, heatmap_legend_param = heatmap_legend_param)
# mmm so this is not working
# just resulting in a box...
str(mat)
# let's trya dn look just at mat and see what has changed
ComplexHeatmap::Heatmap(mat, name = "mat")
# OKAY WHAT IS HAPPENING
# that is not at all what this is supposed to be happening... but very cool
# eveyrtime I run it there are different colours...
pdf(file = paste(path.figures,"funky.lung.2.pdf", sep="/"))
ComplexHeatmap::Heatmap(mat, name = "mat")
dev.off()


# so at this point, basically everything in this tutorial has stopped working
# I am unclear as to why because all of the logical operations are not different
# BUTTTTTTTTTTT they aren't working 
# or even if i copy and paste a command it is not working













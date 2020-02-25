#
# Data Visualization: An exploration of HEATMAPSSSS
# Created for presentation on Tuesday by Sam Morassutti
# this document is intended to cover the basics of building and analyizing heatmaps in R
# *** V basic! 



#-----notes on matrix-----
# similar to  vector: which is a 1D of data elements
# however, a matrix consists of a 2D array of data elements (rows and columns)
# sort of an extension
# it can only contan one atomic vector type so you can't mutiple
# variable types in the matrix
#SIDE note* you can add these to your heatmap later though
# so no stress



# to start
#----- using matrix function
# need a vector containing values of interest that will be in matrix
# and a dimension (at least one), so you can specify either the
# number of rows or columns

#-----using rbind and cbind
# can also create a matrix by pasting the rows and columns together
# these functions can also be used to add a row or a column to a
# matrix that already exists (cool and handy!)

#----naming matrix---
# can use the functions rownames() and colnames
# rownames(matrix) <("row1", "row2", "row3")

#-----small small matrix-------
# to start, let's create a 3x4 matrix

matrix.small <- matrix(1:12, ncol = 4)
# the input length is 12, adn the number of columns are 4
# check to make sure
matrix.small
# could also use nrow and see how this changes
matrix.small <- matrix(1:12, nrow = 4)
matrix.small
# nice!
# however, for both these cases r fills in the number column by column
# to do this by row use byrow argument and set = TRUE
matrix.small <- matrix(1:12, nrow = 4, byrow = TRUE)
# check again to see how this changes
matrix.small
# this random matrix will have three groups by columns
# and three groups by rows


# rbind and cbind
cbind(1:6, 1:6, 1:6)
rbind(1:6, 1:6, 1:6)

# so if we want to add another column to matrix.small with the values 1-4
matrix.small <- cbind(matrix.small, 1:4)
matrix.small
# great!

#----name row/col-------
# now let's get some names on these bad boys
# we have 4 columns and 4 rows
rownames(matrix.small) <- c("row1", "row2", "row3", "row4")
colnames(matrix.small) <- c("col1", "col2", "col3", "col4")
# always check 
matrix.small
# sweet


#----playing with default plots
# use basic default
heatmap(matrix.small)
# wooo now we see a plot with the default colours and clustering of values 
# clustering is super cool and is a VERY important part of heatmaps
# groups values based on similarity
# this has A LOT of uses 


# okay, back to this little weird small matrix heatmap
# let's get a legend on there and some titles on here
# NOW we gotta get some packages in here
# default heatmap function does not have the same capabilities as other packages
# complexheatmap is used a lot


#---------libraries----------------------------------
# libraries needed (a lot)
source("https://bioconductor.org/biocLite.R")
biocLite("ComplexHeatmap")
library(ComplexHeatmap)
library(BiocManager)
#BiocManager::install("ComplexHeatmap")
library(dplyr)
library(devtools)
# and the package complex heatmap 
install_github("jokergoo/ComplexHeatmap")
library(circlize)
library(cluster)

# still having troubles, but we can do some of this without all the right 
# packages

# check the default setting as it is different in this package
ComplexHeatmap::Heatmap(matrix.small)
# isn't that cute!
# are all positive between the input values are positive
#default clustering has changed the order of the col/rows
# this groups values together

#add a legend title use the name = ""
ComplexHeatmap::Heatmap(matrix.small, name = "matrix.small")
# notice how the plot has been clustered based on similarity of values
# the dendrograms groups these relationships

#---titles-----
# to create a title 
# use the command row_title_, or column_title_
ComplexHeatmap::Heatmap(matrix.small, name = "matrix.small", row_title = "Title of rows", 
                        column_title = "Title of columns")



#-----clustering---
# can also turn clustering off, do for rows and then columns
# depending on if clustering is needed or not
# let's get rid of default clustering 
# used function cluster_rows or cluster_columns = FALSE
ComplexHeatmap::Heatmap(matrix.small, name = "matrix.small", 
                        row_title = "Title of rows", 
                        column_title = "Title of columns", cluster_rows = FALSE)
# first let's remove the clustering of rows
# it's still kinda similar
# let's see what happens when we remove all clustering
ComplexHeatmap::Heatmap(matrix.small, name = "matrix.small", 
                        row_title = "Title of rows", 
                        column_title = "Title of columns", cluster_rows = FALSE, 
                        cluster_columns = FALSE)

# looks pretty different!

#clustering can also been shown on specific sides of heatmap
# for example, t show clustering on the right side of rows, and 
# the bottom of the columns
# to do this use the function row_dend_side = "", column_dend_side
ComplexHeatmap:: Heatmap(matrix.small, name = "matrix.small", 
                         row_title = "Title of rows", 
                         column_title = "Title of columns",
                         row_dend_side = "right", 
                         column_dend_side = "bottom")

#------------Change the colours---------------------------------------
# let's play around with some colours!
# going to use colorRamp2 to do this
# if matrix is continuous, provide a vector of colours 
# the colours will be interpolated linearly
# however, mapping isn't robust against outliers
# becasue it starts with the min and max values
# in the matrix

# therefore, if you set the colour to the max/min
# of the matrix, it will be identical to the
# default plot

# let's change some colours!!!
# create a function and call is col_fun
# input the range of values (we already know these becasue we created them)
# but can also get quickly by the min and max commands
min(matrix.small) #1
max(matrix.small) #12

# because our range is from 1-12, we can select three colours and assign them 
# to the min, mid, max 
col_fun = colorRamp2(c(1, 6, 12), c("mediumvioletred", "white", "mediumspringgreen"))
# now we can input the sequence (how many values in our colour range)
col.fun <- col_fun(seq(1, 12))
ComplexHeatmap::Heatmap(matrix.small, name = "matrix.small", col = col.fun, 
                        row_title = "Title of rows", 
                        column_title = "Title of columns")

# see how this plot show the spread in values much better than the default colours
# the legend is much nicer 
# the clustering also creates groups 







#----------Bigger random matrix------------------------
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
# to get the lendend and legend name
# use function name = "", then put the name of legend wooooo
ComplexHeatmap::Heatmap(mat, name = "mat")

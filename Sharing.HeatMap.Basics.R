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

#----name row/col
# now let's get some names on these bad boys
# we have 4 columns and 4 rows
rownames(matrix.small) <- c("row1", "row2", "row3", "row4")
colnames(matrix.small) <- c("col1", "col2", "col3", "col4")
# always check 
matrix.small
# sweet


#----playing with default plots)
# use basic default
heatmap(matrix.small)
# wooo now we see a plot with the default colours and clustering of values 
# clustering is super cool and is a VERY important part of heatmaps
# groups values based on similarity
# this has A LOT of uses 


#okay, back to this little weird small matrix heatmap
# let's get a legend on there
# NOW we gotta get some packages in here


#---------libraries----------------------------------
# libraries needed (a lot)
source("https://bioconductor.org/biocLite.R")
biocLite("ComplexHeatmap")
library(ComplexHeatmap)

library(BiocManager)
#BiocManager::install("ComplexHeatmap")

# still having troubles, but we can do some of this without all the right 
# packages

# check the default setting as it is different in this package
ComplexHeatmap::Heatmap(matrix.small)
# isn't that cute!
#---------------------------------------------------


# are all positive between the input values were







# if matrix is continuous, provide a vector of colours 
# the colours will be interpolated linearly
# however, mapping isn't robust against outliers
# becasue it starts with the min and max values
# in the matrix

# therefore, if you set the colour to the max/min
# of the matrix, it will be identical to the
# default plot





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

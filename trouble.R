#
# This document is intended to contain the code that is NOT working
# it should be run after the entire MAIN script file of this project has been
# Athored by Sam Morassutti
# February 2020


#-------------------------------------------------------------------
#-----------plotting mutations as colours???--------------------------
#--------------------------------------------------------------------
?oncoPrint()

col = c("MUT" = "#008000", "AMP" = "red", "HOMDEL" = "blue")
column_title = "OncoPrint for TCGA Lung Adenocarcinoma, genes in Ras Raf MEK JNK signalling"
heatmap_legend_param = list(title = "Alternations", at = c("HOMDEL", "AMP", "MUT"), 
                            labels = c("Deep deletion", "Amplification", "Mutation"))
oncoPrint(mat,
          alter_fun = alter_fun, col = col, 
          column_title = column_title, heatmap_legend_param = heatmap_legend_param)

oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "Lung Adenocarcinoma genes in Ras Raf MEK JNK signalling",
          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "HOMDEL", "MUT"), 
                                      labels = c("Amplification", "Deep deletion", "Mutation")))
# want to remove NA values from the matrix 
oncoPrint(mat,
          alter_fun = alter_fun, col = col, 
          remove_empty_columns = TRUE, remove_empty_rows = TRUE,
          column_title = column_title, heatmap_legend_param = heatmap_legend_param)

# I am pretty stuck..... but learning has happened along the way!



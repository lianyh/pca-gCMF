# pca-gCMF
------------------------------

# *Prerequisite*
-------------------------------
Please install R library:<br/>
library(CMF)<br/>
library(pROC)<br/>
library(data.table)<br/>
library(igraph)<br/>
library(sna)<br/>


Command
--------------------------
Rscript pca_gCMF.R matrix_list_pca.txt data/ outputFolder/<br/>

1st input:  the file contains a list of :<br/>
            1st column: input matrices name. <br/>
            2nd column: binary value {1,0} whether the input martrix file required PCA transformation.<br/>
            3rd column: entity type corresponds to the row in the matrix (entity type row index).<br/>
            4th column: entity type corresponds to the column in the matrix (entity type column index).<br/>
2nd input: folder contains a list of files listed in the first argument<br/>
3rd input: your output folder ended with slash / <br/>

All columns are separated by tab delimited.<br/>

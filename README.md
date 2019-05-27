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

matrix_list_pca.txt file contains:<br/>
            -1st column: name list of input matrices. <br/>
            --2nd column: binary value {1,0} whether the input martrix file required PCA transformation.<br/>
            3rd column: entity type corresponds to the row in the matrix (entity type row index).<br/>
            -4th column: entity type corresponds to the column in the matrix (entity type column index).<br/>
data/ folder contains a list of files listed in the first argument<br/>
outputFolder/ is your output folder ended with slash / <br/>

All columns are separated by tab delimited.<br/>

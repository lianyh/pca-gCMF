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
<pre>       1st column: name list of input matrices. </pre><br/>
<pre>       2nd column: binary value {1,0} whether the input martrix file required PCA transformation.</pre>
<pre>       3rd column: entity type corresponds to the row in the matrix (entity type row index).</pre>
<pre>       4th column: entity type corresponds to the column in the matrix (entity type column index).</pre>
data/ folder contains a list of files listed in the first argument<br/>
outputFolder/ is your output folder ended with slash / <br/>

All columns are separated by tab delimited.<br/>

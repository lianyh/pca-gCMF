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
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1st column: name list of input matrices.<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;2nd column: binary value {1,0} whether the input martrix file required PCA transformation.<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;3rd column: entity type (index number) corresponds to the row in the matrix (entity type row index).<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;4th column: entity type (index number) corresponds to the column in the matrix (entity type column index).<br/>
data/ folder contains a list of files listed in the first argument<br/>
outputFolder/ is your output folder ended with slash / <br/>

All columns are separated by tab delimited.<br/>



Citation
--------------------------------
If you find this useful for your research, we would be appreciated if you cite the following papers:

@article{liany2020predicting,
  title={Predicting synthetic lethal interactions using heterogeneous data sources},
  author={Liany, Herty and Jeyasekharan, Anand and Rajan, Vaibhav},
  journal={Bioinformatics},
  volume={36},
  number={7},
  pages={2209--2216},
  year={2020},
  publisher={Oxford University Press}
}


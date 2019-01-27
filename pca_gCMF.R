library(CMF)
library(pROC)
library(data.table)

args<-commandArgs(TRUE)
if(length(args) < 7) {
  args <- c("--help")
}
library(CMF)
library(pROC)
library(data.table)
args<-commandArgs(TRUE)
if(length(args) < 7) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      The R Script
 
      Arguments:
      --arg1 = Known SL Matrix (gene * gene)   - a file
      --arg2 = Test Set - SL Matrix (gene * gene) - a file
      --arg3 = Transformed (PCA) Essentiality Profile (366 * 30)    - a file
      --arg4 = Transformed (PCA) Pairwise Coexpr Profile (366 * 192)    - a file
      --arg5 = Raw RNA Expression Profile (366 * 1100)    - a file
      --arg6 = Transformed (PCA) SCNA Profile (366 * 500)    - a file
      --arg7 = outputFolder/    - folder path ended with /
      --help             
 
      Example:
      ./pca-gCMF.R --arg1=F1_F2_F3_SL_binary_all --arg2=F1_SL_binary_test --arg3=F1_F2_F3_essentiality_pca --arg4=F1_F2_F3_pairwisecoexpr_pca192 --arg5=data_RNA_expression_alltraintest --arg6=data_linear_CNA_alltraintest_pca --arg7=outputFolder/ \n\n")
 
  q(save="no")
}
SL_binary_all=args[1]
test_set=args[2]
essentiality_pca_all=args[3]
pairwisecoexpr_pca_all=args[4]
raw_data_RNA_expression=args[5]
sca_pca_all=args[6]
out=args[7]

if (!(file.exists(out)))
{
	dir.create(out)
} 

triplets=list()
X=list()
F1_F2_F3_SL_binary_all=read.table(SL_binary_all,sep="\t")
F1_F2_F3_SL_binary_all=as.matrix(F1_F2_F3_SL_binary_all)
colcount1=ncol(F1_F2_F3_SL_binary_all)

X[[1]]<-matrix(F1_F2_F3_SL_binary_all,nrow=colcount1,ncol=colcount1)
triplets[[1]]=matrix_to_triplets(X[[1]])

triplets_test=list()
Y=list()

F3_SL_binary_test=read.table(test_set,sep="\t",header=0)
F3_SL_binary_test=as.matrix(F3_SL_binary_test)
Y[[1]]<-matrix(F3_SL_binary_test,nrow=colcount1,ncol=colcount1)
triplets_test[[1]]=matrix_to_triplets(Y[[1]])

train=list()
test=list()
m1=triplets[[1]]
m2=triplets_test[[1]]

m1 <- data.table(m1)
setkey(m1, 'V1', 'V2')
m1[,"index1" := .I]
m2 <- data.table(m2)
setkey(m2, 'V1', 'V2')
m2[,"index2" := .I]

# Join the tables by key #
m3 <- m1[m2]

overlap <- m3[is.na(index1)==FALSE & is.na(index2)==FALSE,]
myIndex=as.matrix(overlap[,4])
train[[1]]=triplets[[1]][-myIndex,]
test[[1]]=triplets[[1]][myIndex,]

myIndex=m2[,V1]

#ADDED ESSENTIALITY
exprfeatures=read.table(essentiality_pca_all,sep="\t",header=0)
colcount2=ncol(exprfeatures)

X[[2]]=matrix(as.matrix(exprfeatures),nrow=colcount2,ncol=colcount2)
triplets[[2]]=matrix_to_triplets(X[[2]])

m2<- data.table(triplets[[2]])
test[[2]]=triplets[[2]][m2[,V1] %in% myIndex,]
 `%not_in%` <- purrr::negate(`%in%`)
train[[2]]=triplets[[2]][m2[,V1] %not_in% myIndex,]

#ADDED PAIRWISE CO-EXPR
exprfeatures=read.table(pairwisecoexpr_pca_all,sep="\t",header=0)
colcount3=ncol(exprfeatures)

X[[3]]=matrix(as.matrix(exprfeatures),nrow=colcount1,ncol=colcount3)
triplets[[3]]=matrix_to_triplets(X[[3]])

m2<- data.table(triplets[[3]])
test[[3]]=triplets[[3]][m2[,V1] %in% myIndex,]
 `%not_in%` <- purrr::negate(`%in%`)
train[[3]]=triplets[[3]][m2[,V1] %not_in% myIndex,]

#ADDED RNA SEQ EXPRESSION PROFILE AND CNA (CBIOPORTAL)
exprfeatures=read.table(raw_data_RNA_expression,sep="\t",header=0)
colcount4=ncol(exprfeatures)

X[[4]]=matrix(as.matrix(exprfeatures),nrow=colcount1,ncol=colcount4)
triplets[[4]]=matrix_to_triplets(X[[4]])

m2<- data.table(triplets[[4]])
test[[4]]=triplets[[4]][m2[,V1] %in% myIndex,]
 `%not_in%` <- purrr::negate(`%in%`)
train[[4]]=triplets[[4]][m2[,V1] %not_in% myIndex,]

cnafeatures=read.table(sca_pca_all,sep="\t",header=0)
colcount5=ncol(cnafeatures)

X[[5]]=matrix(as.matrix(cnafeatures),nrow=colcount1,ncol=colcount5)
triplets[[5]]=matrix_to_triplets(X[[5]])
m2<- data.table(triplets[[5]])
test[[5]]=triplets[[5]][m2[,V1] %in% myIndex,]
 `%not_in%` <- purrr::negate(`%in%`)
train[[5]]=triplets[[5]][m2[,V1] %not_in% myIndex,]

K=5
inds <- matrix(0,nrow=5,ncol=2)
inds[1,]=c(1,1)
inds[2,]=c(1,2)
inds[3,]=c(1,3)
inds[4,]=c(1,4)
inds[5,]=c(1,5)


D=c(colcount1,colcount2,colcount3,colcount4,colcount5)
likelihood <- c("gaussian","gaussian","gaussian","gaussian","gaussian")
opts <- getCMFopts()
opts$iter.max <- 10 # Less iterations for faster computation
model <- CMF(train,inds,K,likelihood,D,test=test,opts=opts)
outm <- predictCMF(test, model)

truth=triplets_test[[1]][,3]

temp_testerr1=model$errors[1,1]
temp_results=outm$out[[1]]

for (i in 1:10) {
	model <- CMF(train,inds,K,likelihood,D,test=test,opts=opts)
	newmodelError=model$errors[1,1]
	#print("newmodelError")
	#print(newmodelError)

	if (newmodelError < temp_testerr1)
	{
		temp_testerr1=newmodelError
		out <- predictCMF(test, model)
		temp_results=out$out[[1]]
	}
}
ROC1 <- roc(as.vector(truth), as.vector(temp_results[,3]))
AUC1 <- auc(ROC1)

print("Test Error:")
print(temp_testerr1)

print("AUC")
print(AUC1)
write.table(temp_results,file=paste(args[7],"gcmf_prediction_results",sep=""),sep="\t",row.names=FALSE,col.names=FALSE)

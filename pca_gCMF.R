library(CMF)
library(pROC)
library(data.table)


args<-commandArgs(TRUE)
if(length(args) < 3) {
  args <- c("--help")
}
## Help section
if("--help" %in% args) {
  cat("
      The R Script
 
      Arguments:
      --arg1 = File contains a list of matrices name - a file
      --arg2 = input folder contains the files listed in argument 1/  - folder path ended with /
      --arg3 = outputFolder/  - folder path ended with /
      --help             
 
      Example:
      ./pca-gCMF.R --arg1=matrix_list_pca.txt --arg2=dat1/ --arg3=outputFolder/ \n\n

      Notes: In argument 1, input a file contains a list of :\n
      1st column: input matrices name. \n
      2nd column: binary value {1,0} whether the input martrix file required PCA transformation.\n
      3rd column: entity type corresponds to the row in the matrix (entity type row index).\n
      4th column: entity type corresponds to the column in the matrix (entity type column index).\n\n

      All columns are separated by tab delimited.\n\n
      ")
  q(save="no")
}
input=args[2]
outFolder=args[3]

if (!(file.exists(outFolder)))
{
	dir.create(outFolder)
} 

fileName <- args[1]
conn <- file(fileName,open="r")
linn <-readLines(conn)
test_set=strsplit(linn[1], '\t') [[1]][1]
train_set=strsplit(linn[2], '\t') [[1]][1]

K=length(linn)-1
size=length(linn) - 2
countList=list()
matrixList=list()
likelihood=vector()
likelihood[1]="gaussian"

inds <- matrix(0,nrow=K,ncol=2)
inds[1,]=c(as.numeric(strsplit(linn[1], '\t') [[1]][3]),as.numeric(strsplit(linn[1], '\t') [[1]][4]))

triplets=list()
X=list()
F1_F2_F3_SL_binary_all=read.table(paste(input,train_set,sep="/"),sep="\t",header=0)
F1_F2_F3_SL_binary_all=as.matrix(F1_F2_F3_SL_binary_all)
countList[[1]]=ncol(F1_F2_F3_SL_binary_all)

X[[1]]<-matrix(F1_F2_F3_SL_binary_all,nrow=countList[[1]],ncol=countList[[1]])
triplets[[1]]=matrix_to_triplets(X[[1]])

triplets_test=list()
Y=list()

F3_SL_binary_test=read.table(paste(input,test_set,sep="/"),sep="\t",header=0)
F3_SL_binary_test=as.matrix(F3_SL_binary_test)

Y[[1]]<-matrix(F3_SL_binary_test,nrow=countList[[1]],ncol=countList[[1]])
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

for (i in 3:length(linn)){

fields =strsplit(linn[i], '\t') [[1]]
if(!is.na(fields[1]))
{
	inds[i-1,]=c(as.numeric(strsplit(linn[i], '\t') [[1]][3]),as.numeric(strsplit(linn[i], '\t') [[1]][4]))
	likelihood[i-1]="gaussian"

	#ADDED ESSENTIALITY
	exprfeatures=read.table(fields[1],sep="\t",header=0)
	colcount2=ncol(exprfeatures)
	countList[[i-1]]=colcount2
	isPCA = fields[2]
		if (isPCA == 1)
		{
			 train1=as.matrix(exprfeatures)
			 pca=prcomp(train1)
			 pcasummary=summary(pca)$importance[3,]
			 pc_keep_len=length(pcasummary[pcasummary<0.99]) #proportion of the variance explained=0.99
			 pca_retain=pca$x[,1:pc_keep_len]
			 countList[[i-1]]=pc_keep_len

			 X[[i-1]]=matrix(as.matrix(pca_retain),nrow=nrow(pca_retain),ncol=ncol(pca_retain))
			 triplets[[i-1]]=matrix_to_triplets(X[[i-1]])
			 m2<- data.table(triplets[[i-1]])
			 test[[i-1]]=triplets[[i-1]][m2[,V1] %in% myIndex,]
			 `%not_in%` <- purrr::negate(`%in%`)
			 train[[i-1]]=triplets[[i-1]][m2[,V1] %not_in% myIndex,]

		}
		else {
			X[[i-1]]=matrix(as.matrix(exprfeatures),nrow=nrow(exprfeatures),ncol=ncol(exprfeatures))
			triplets[[i-1]]=matrix_to_triplets(X[[i-1]])

			m2<- data.table(triplets[[i-1]])
			test[[i-1]]=triplets[[i-1]][m2[,V1] %in% myIndex,]
			 `%not_in%` <- purrr::negate(`%in%`)
			train[[i-1]]=triplets[[i-1]][m2[,V1] %not_in% myIndex,]
		}
	}
}

close(conn)

D=as.numeric(countList)
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
	if (newmodelError < temp_testerr1)
	{
		temp_testerr1=newmodelError
		out <- predictCMF(test, model)
		temp_results=out$out[[1]]
	}

}
ROC1 <- roc(as.vector(truth), as.vector(temp_results[,3]))
AUC1 <- auc(ROC1)

#print("Test Error:")
#print(temp_testerr1)

print("AUC")
print(AUC1)
write.table(temp_results,file=paste(outFolder,"gcmf_pca_based_prediction_results",sep="/"),sep="\t",row.names=FALSE,col.names=FALSE)
cat("\n")

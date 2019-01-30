library(CMF)
library(pROC)
library(data.table)
library('igraph')
library('sna')

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
      ./pca-gCMF.R --arg1=matrix_list_gb.txt --arg2=dat2/ --arg3=outputFolder/ \n\n

      Notes: In argument 1, please list the matrix meant for prediction (test set) in the first row and \n \t \t its related full matrix for training in the second row. \n
      In argument 1, the second column indicates whether a matrix required graph-based features transformation {1,0}, \n \t \t separated by tab.\n
      In argument 1, the third and fourth column indicate entity type, separated by tab, i.e. 1   2.\n\n
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
	print("fields 1")
	print(fields[1])
	#ADDED ESSENTIALITY
	exprfeatures=read.table(fields[1],sep="\t",header=0)
	colcount2=ncol(exprfeatures)
	countList[[i-1]]=colcount2
	isGP = fields[2]
		if (isGP == 1)
		{
			 tableMatrix=as.matrix(exprfeatures)

			 gTable=matrix(tableMatrix,nrow=nrow(tableMatrix),ncol=ncol(tableMatrix))
			
			 deg=degree(gTable, gmode="graph")
			 gDeg=as.matrix(deg,ncol=1)

			 tt=closeness(gTable,gmode="graph")
			 gCloness=as.matrix(tt,ncol=1)
			 
			 tt=betweenness(gTable,gmode="graph")
			 gBetweenness=as.matrix(tt,ncol=1)

			 tt=stresscent(gTable,gmode="graph")
			 gStressCent=as.matrix(tt,ncol=1)

			 tt=infocent(gTable,gmode="graph")
			 gInfoCent=as.matrix(tt,ncol=1)

			 tt=evcent(gTable, gmode="graph")
			 gEvCent=as.matrix(tt,ncol=1)

			 tt=gilschmidt(gTable, gmode="graph")
			 gGilSchmidt=as.matrix(tt,ncol=1)

			 tt=flowbet(gTable, gmode="graph")
			 gFlowBet=as.matrix(tt,ncol=1)

			 graphFeatures=cbind(gDeg,gCloness,gBetweenness,gStressCent,gInfoCent,gEvCent,gGilSchmidt,gFlowBet)

			 pca_retain=ncol(graphFeatures)
			 countList[[i-1]]=pca_retain

			 X[[i-1]]=matrix(as.matrix(graphFeatures),nrow=nrow(tableMatrix),ncol=ncol(graphFeatures))
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
write.table(temp_results,file=paste(outFolder,"gcmf_graph_based_prediction_results",sep="/"),sep="\t",row.names=FALSE,col.names=FALSE)
cat("\n")

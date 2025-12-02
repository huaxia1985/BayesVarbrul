pathname <- "~/Desktop/Shawi" 
setwd(pathname)

#data 
rawdata <- read.csv("rawdata.csv",check.names=F)
library(stringr)
tmp <- str_sub(colnames(rawdata),2)
L <- as.numeric(tmp[length(tmp)])
data <- NULL
data$Xind <- rbind(rawdata$region=="Cahuapanas",rawdata$region=="Sillay",rawdata$gender=="M")*1
data$Xind <- rbind(data$Xind,rawdata$age)
rownames(data$Xind) <- c("Ca","Sa","male","age")
data$Xind <- data$Xind[,-1]
data$lingua <- vector("list",L)
for (i in 1:L) {
    data$lingua[[i]] <- rawdata[-1,which(tmp==as.character(i))]
}
data$type <- vector("list",L)
for (i in 1:L) {
	data$type[[i]] <- as.numeric(rawdata[1,which(tmp==as.character(i))])
}
data$theta0 <- vector("list",L)
for (i in 1:L) {
	data$theta0[[i]] <- numeric(length(data$type[[i]]))
	tmp2 <- which(data$type[[i]]==1)
	data$theta0[[i]][tmp2] <- 1/length(tmp2)
}

##calculate the mean of posterior samples of theta for each variable
nspeaker <- 168
nvar <- 20
burin <- 50
TT <- 2637
theta <- vector("list",nvar)
for (ii in 1:nvar) {
	filename <- paste0("theta",ii,".csv")
	tmp <- read.table(filename,header=F,sep=" ",fill=T,stringsAsFactors=F)
	idx <- which(is.na(tmp[,nspeaker+1]))
	i <- burin+1
	j <- i-burin
	tmp2 <- tmp[(idx[i]+1):(ifelse(i<length(idx),idx[i+1]-1,dim(tmp)[1])),-1]
   	theta[[ii]] <- data.frame(apply(tmp2,2,function (x) as.numeric(x)))
	for (i in (burin+2):TT) {
		j <- i-burin
		tmp2 <- tmp[(idx[i]+1):(ifelse(i<length(idx),idx[i+1]-1,dim(tmp)[1])),-1]
    	theta[[ii]] <- theta[[ii]] + data.frame(apply(tmp2,2,function (x) as.numeric(x)))
    }
    theta[[ii]] <- theta[[ii]]/(TT-burin)
}

#plot PCA
region <- rep("Balsapuerto",nspeaker)
region[data$Xind["Ca",]==1] <- "Cahuapanas"
region[data$Xind["Sa",]==1] <- "Sillay"
color <- rep("black",nspeaker)
color[region=="Cahuapanas"] <- "red"
color[region=="Sillay"] <- "blue"
age <- data$Xind["age",]

library("corrplot")
library("ggplot2")
library("ggfortify")

pca.data <- t(matrix(as.numeric(theta[[1]][-dim(theta[[1]])[1],]),ncol=nspeaker))
for (i in 2:nvar) {
	pca.data <- cbind(pca.data,t(matrix(as.numeric(theta[[i]][-dim(theta[[i]])[1],]),ncol=nspeaker)))
}
var.name <- c("3PL","PROG","because","NMLZ","SEQ","dog","DIM","isu'","SOV","SV","spider","furry catepillar","vowel length","1st person marker","chicken","tasty","jump","all","search for","DESID,REC")
var <- unlist(lapply(1:nvar,function (i) paste0(var.name[i],".",c(1:(dim(theta[[i]])[1]-1)))))
colnames(pca.data) <- var
res.pca <- prcomp(pca.data,scale=T)
autoplot(res.pca, label=TRUE, loadings=TRUE, loadings.label=TRUE, loadings.label.size=5, loadings.colour="dark grey", loadings.label.colour="black",colour=color, size=age/5)+
  theme_bw()

#plot prediction vs observation for each langauge variable
library(klaR)
par(mfrow=c(5,6))
clab <- as.numeric(as.factor(rawdata$region[-1]))
nl <- sapply(data$type,function (i) length(i))
for (i in 1:nvar) {
theta.post <- t(as.matrix(theta[[i]]))
tmp <- rawdata[,which(colnames(rawdata)==paste0("V",i))]
idx <- which(!is.na(rowSums(tmp)))[-1]
tmp <- tmp[idx,]
if (nl[i]>4) {
if (data$type[[i]][nl[i]]==4) {
	tmp <- tmp[,-nl[i]]
	tmp <- tmp/rowSums(tmp)
	tmp <- as.matrix(tmp)
	tmp[is.na(tmp[,1]),] <- 1/(nl[i]-1)
	tmp <- tmp + abs(matrix(rnorm(length(tmp),0,0.1),dim(tmp)))
	tmp <- tmp/rowSums(tmp)
	theta.post2 <- theta.post[,-nl[i]]
	theta.post2 <- theta.post2 + abs(matrix(rnorm(length(theta.post2),0,0.1),dim(theta.post2)))
	theta.post2 <- theta.post2/rowSums(theta.post2)
	quadplot(theta.post2[idx-1,],col=c("red","blue","green","purple","grey")[clab[idx]],pch=16,label=c("v1","v2","v3","v4"),labelpch=17)
	quadplot(tmp,col=c("red","blue","green","purple","grey")[clab[idx]],pch=16,label=c("v1","v2","v3","v4"),labelpch=17)
} else {
	tmp <- tmp/rowSums(tmp)
	tmp <- as.matrix(tmp)
	tmp[is.na(tmp[,1]),] <- 1/nl[i]
	tmp <- tmp + abs(matrix(rnorm(length(tmp),0,0.01),dim(tmp)))
	tmp <- tmp/rowSums(tmp)
	theta.post2 <- theta.post[,1:4]
	theta.post2 <- theta.post2 + abs(matrix(rnorm(length(theta.post2),0,0.01),dim(theta.post2)))
	theta.post2 <- theta.post2/rowSums(theta.post2)
	tmp <- tmp[,1:4]
	tmp <- tmp + abs(matrix(rnorm(length(tmp),0,0.01),dim(tmp)))
	tmp <- tmp/rowSums(tmp)
	quadplot(theta.post2[idx-1,],col=c("red","blue","green","purple","grey")[clab[idx]],pch=16,label=c("v1","v2","v3","v4"),labelpch=17)
	quadplot(tmp,col=c("red","blue","green","purple","grey")[clab[idx]],pch=16,label=c("v1","v2","v3","v4"),labelpch=17)
} 
}
if (nl[i]==4) {
if (data$type[[i]][nl[i]]==4) {
	tmp <- tmp[,-4]
	tmp <- tmp/rowSums(tmp)
	tmp <- as.matrix(tmp)
	tmp[is.na(tmp[,1]),] <- 1/(nl[i]-1)
	tmp <- tmp + abs(matrix(rnorm(length(tmp),0,0.1),dim(tmp)))
	tmp <- tmp/rowSums(tmp)
	theta.post2 <- theta.post[,-4]
	theta.post2 <- theta.post2 + abs(matrix(rnorm(length(theta.post2),0,0.1),dim(theta.post2)))
	theta.post2 <- theta.post2/rowSums(theta.post2)
	triplot(theta.post2[idx-1,],col=c("red","blue","green","purple","grey")[clab[idx]],pch=16,label=c("v1","v2","v3"),grid=F,center=F)
	triplot(tmp,col=c("red","blue","green","purple","grey")[clab[idx]],pch=16,label=c("v1","v2","v3"),grid=F,center=F)
} else {
	tmp <- tmp/rowSums(tmp)
	tmp <- as.matrix(tmp)
	tmp[is.na(tmp[,1]),] <- 1/nl[i]
	tmp <- tmp + abs(matrix(rnorm(length(tmp),0,0.01),dim(tmp)))
	tmp <- tmp/rowSums(tmp)
	theta.post2 <- theta.post
	theta.post2 <- theta.post2 + abs(matrix(rnorm(length(theta.post2),0,0.01),dim(theta.post2)))
	theta.post2 <- theta.post2/rowSums(theta.post2)
	quadplot(theta.post2[idx-1,],col=c("red","blue","green","purple","grey")[clab[idx]],pch=16,label=c("v1","v2","v3","v4"),labelpch=17)
	quadplot(tmp,col=c("red","blue","green","purple","grey")[clab[idx]],pch=16,label=c("v1","v2","v3","v4"),labelpch=17)
}
} 
if (nl[i]==3) {
if (data$type[[i]][nl[i]]==4) {
	tmp <- tmp[,-3]
	tmp <- tmp/rowSums(tmp)
	tmp <- as.matrix(tmp)
	tmp[is.na(tmp[,1]),] <- 1/(nl[i]-1)
	tmp <- tmp + abs(matrix(rnorm(length(tmp),0,0.1),dim(tmp)))
	tmp <- tmp/rowSums(tmp)
	theta.post2 <- theta.post[,-3]
	theta.post2 <- theta.post2 + abs(matrix(rnorm(length(theta.post2),0,0.1),dim(theta.post2)))
	theta.post2 <- theta.post2/rowSums(theta.post2)
	plot(tmp[,1],theta.post2[idx-1,1],col=c("red","blue","green","purple","grey")[clab[idx]],pch=16,xlab="observed frequency of nitun",ylab="predicted frequency of nitun")
} else {
	tmp <- tmp/rowSums(tmp)
	tmp <- as.matrix(tmp)
	tmp[is.na(tmp[,1]),] <- 1/nl[i]
	tmp <- tmp + abs(matrix(rnorm(length(tmp),0,0.01),dim(tmp)))
	tmp <- tmp/rowSums(tmp)
	theta.post2 <- theta.post
	theta.post2 <- theta.post2 + abs(matrix(rnorm(length(theta.post2),0,0.01),dim(theta.post2)))
	theta.post2 <- theta.post2/rowSums(theta.post2)
	triplot(theta.post2[idx-1,],col=c("red","blue","green","purple","grey")[clab[idx]],pch=16,label=c("v1","v2","v3"),grid=F,center=F)
	triplot(tmp,col=c("red","blue","green","purple","grey")[clab[idx]],pch=16,label=c("v1","v2","v3"),grid=F,center=F)
}
}
if (nl[i]==2) {
	tmp <- tmp/rowSums(tmp)
	tmp <- as.matrix(tmp)
	tmp[is.na(tmp[,1]),] <- 1/nl[i]
	tmp <- tmp + abs(matrix(rnorm(length(tmp),0,0.01),dim(tmp)))
	tmp <- tmp/rowSums(tmp)
	theta.post2 <- theta.post
	theta.post2 <- theta.post2 + abs(matrix(rnorm(length(theta.post2),0,0.01),dim(theta.post2)))
	theta.post2 <- theta.post2/rowSums(theta.post2)
	plot(tmp[,1],theta.post2[idx-1,1],col=c("red","blue","green","purple","grey")[clab[idx]],pch=16,xlab="observed frequency of nitun",ylab="predicted frequency of nitun")
}
}
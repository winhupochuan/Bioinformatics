setwd("E:/R/homework1/")
marker <- read.table("marker data 20% missing(6187loci).txt")
marker <- t(marker)
marker1 <- marker[1:3,]
phenotype <- read.table("zhugenchang_16.txt")
pheno <- matrix(0,nrow=321,ncol=2)
pheno[,1:2] <- c(phenotype[,1],phenotype[,17])
marker <- subset(marker,marker[,1]%in%pheno[,1])
pheno <- subset(pheno,pheno[,1]%in%marker[,1])
a <- matrix(0,nrow=3,ncol=2)
pheno <- rbind(a,pheno)
marker[which(marker=="ll")] <- 0
marker[which(marker=="lm")] <- 1
marker[which(marker=="nn")] <- 0
marker[which(marker=="np")] <- 1
marker[which(marker=="hh")] <- 0
marker[which(marker=="hk")] <- 1
marker[which(marker=="kk")] <- 2
marker[which(marker=="--")] <- 9
marker <- rbind(marker1,marker)
rt <- cbind(pheno,marker[,-1])
#---------------------------------L0计算出来放到一个向量里-----------
calc_L0 <- function (pheno,geno){
  miu <- mean(pheno)
  sigma2 <- sum((pheno-miu)^2)/length(pheno)
  fy0 <- dnorm(pheno,miu, sqrt(sigma2), log=TRUE)
  L0 <- sum(fy0)
  return(L0)}
calc_L1 <- function(pheno,geno){
  pheno0 <- pheno[which(geno==0)]
  pheno1 <- pheno[which(geno==1)]
  pheno2 <- pheno[which(geno==2)]
  miu0 <- mean(pheno0)
  miu1 <- mean(pheno1)
  miu2 <- mean(pheno2)
  n0<-length(pheno0)
  n1<-length(pheno1)
  n2<-length(pheno2)
  if (sum(unique(geno))==1){
    sigma2 <- (sum((pheno0-miu0)^2)+sum((pheno1-miu1)^2))/(n0+n1)
    fy0 <- dnorm(pheno0,miu0, sqrt(sigma2), log=T)
    fy1 <- dnorm(pheno1,miu1, sqrt(sigma2), log=T)
    L1 <-  sum(fy0)+sum(fy1)
  }
  if (sum(unique(geno))==2){
    sigma2 <- (sum((pheno0-miu0)^2)+sum((pheno2-miu2)^2))/(n0+n2)
    fy0 <- dnorm(pheno0,miu0, sqrt(sigma2), log=T)
    fy2 <- dnorm(pheno2,miu2, sqrt(sigma2), log=T)
    L1 <- sum(fy0)+sum(fy2)
  }
  if (sum(unique(geno))==3 & length(unique(geno))==2){
    sigma2 <- (sum((pheno2-miu2)^2)+sum((pheno1-miu1)^2))/(n2+n1)
    fy1 <- dnorm(pheno1,miu1, sqrt(sigma2), log=T)
    fy2 <- dnorm(pheno2,miu2, sqrt(sigma2), log=T)
    L1 <- sum(fy2)+sum(fy1)
  }
  if(length(unique(geno))==3){
    sigma2 <- (sum((pheno0-miu0)^2)+sum((pheno1-miu1)^2)+sum((pheno2-miu2)^2))/(n0+n1+n2)
    fy0 <- dnorm(pheno0,miu0, sqrt(sigma2), log=T)
    fy1 <- dnorm(pheno1,miu1, sqrt(sigma2), log=T)
    fy2 <- dnorm(pheno2,miu2, sqrt(sigma2), log=T)
    L1 <- sum(fy0)+sum(fy1)+sum(fy2)} 
  return(L1)     
}
#----------------------------------------------------------------------------

N_snp <- dim(rt)[2]-2
pv <- rep(1,times=N_snp)
cur_geno<-rep(1,times=320)
for (j in 3:dim(rt)[2]){
  cur_geno <- rt[4:dim(rt)[1],j]
  cur_pheno <- rt[4:dim(rt)[1],2]
  misid<-which(cur_geno==9)
  cur_geno <- cur_geno[-misid]
  cur_pheno <- cur_pheno[-misid]
  cur_geno <- as.numeric(cur_geno)
  cur_pheno <- as.numeric(cur_pheno)
  L0 <- calc_L0(cur_pheno,cur_geno)
  L1 <- calc_L1(cur_pheno,cur_geno)
  LR <- (-2)*(L0-L1)
  df <- length(unique(cur_geno))-1
  pv[j-2] <- pchisq(LR,df,lower.tail=F)
}
npv <- -log(pv)
#--------------------------------------------------------------------------
pdf("ANOVA.pdf",25,10)
plot(0,0,xlim=c(0,6187),ylim=c(0,8),xaxt="n",font.lab=1,cex.lab=1,las=1,tck=0,xaxs="i",yaxs="i",xlab="SNPs",ylab="-log(P)") 
#------------------------------------
colors <- c(3,4,7,6,5)
colorset <- c(colors,colors,colors,colors[1:4])
colors <- rep(0,times=6187)
lg_inf <- rt[1,3:dim(rt)[2]]
type_length <- rep(0,times=19)
lg_type<- rep(0,times=19)
for(i in 1:19){
  lg_type[i] <- paste("lg", i, sep="")
  type_id <- which(lg_inf==lg_type[i])
  colors[type_id] <- colorset[i]
  type_length[i] <- length(type_id)
}
#--------------------------------------------------------
type_length <- c(0,cumsum(type_length))
labels_id <- rep(0,times=19)
for(i in 1:19){
  labels_id[i] <- (type_length[i]+type_length[i+1])/2
}
abline(h=seq(2,6,by=2),v=type_length,col="gray",lwd="2")
axis(side=1,at=labels_id,labels=lg_type,pos=0,tcl=0)
#-----------------------------------
i <- c(1:6187)
points(i,npv[i],col=colors[i],pch=16,cex=.5)
abline(-log(.01),0,col="red")
box(col="black")
dev.off()

#rm(list=ls(all=TRUE))
hadoop <- 0;
library(qtl);

if(!hadoop){ 
	source("C:/Users/Behrang/Documents/Revolution/Prunedirect-V4/calc_uboundB.R")
	source("C:/Users/Behrang/Documents/Revolution/Prunedirect-V4/calc_lboundB.R")
	source("C:/Users/Behrang/Documents/Revolution/Prunedirect-V4/CallConstraintsB.R")
	source("C:/Users/Behrang/Documents/Revolution/Prunedirect-V4/CallObjFcnB.R")
	source("C:/Users/Behrang/Documents/Revolution/Prunedirect-V4/DetermineFcnTypeB.R")
	source("C:/Users/Behrang/Documents/Revolution/Prunedirect-V4/DIRdivideB.R")
	source("C:/Users/Behrang/Documents/Revolution/Prunedirect-V4/DIRiniB.R")
	source("C:/Users/Behrang/Documents/Revolution/Prunedirect-V4/find_poB.R")
	source("C:/Users/Behrang/Documents/Revolution/Prunedirect-V4/mydirectB.R")
	source("C:/Users/Behrang/Documents/Revolution/Prunedirect-V4/replaceinfB.R")
	source("C:/Users/Behrang/Documents/Revolution/Prunedirect-V4/ObjectiveFB.R")
	source("C:/Users/Behrang/Documents/Revolution/Prunedirect-V4/probbcB.R")
	source("C:/Users/Behrang/Documents/Revolution/Prunedirect-V4/rssB.R")
	datafile <- read.cross("csv","C:/Users/Behrang/Documents/Revolution/Prunedirect-V4", "hyper2.csv",genotypes=c("BB","BA"));
}else{
	source("/home/hduser/qtlTest/mapreduce-direct-version/Prunedirect-V4/calc_uboundB.R")
	source("/home/hduser/qtlTest/mapreduce-direct-version/Prunedirect-V4/calc_lboundB.R")
	source("/home/hduser/qtlTest/mapreduce-direct-version/Prunedirect-V4/CallConstraintsB.R")
	source("/home/hduser/qtlTest/mapreduce-direct-version/Prunedirect-V4/CallObjFcnB.R")
	source("/home/hduser/qtlTest/mapreduce-direct-version/Prunedirect-V4/DetermineFcnTypeB.R")
	source("/home/hduser/qtlTest/mapreduce-direct-version/Prunedirect-V4/DIRdivideB.R")
	source("/home/hduser/qtlTest/mapreduce-direct-version/Prunedirect-V4/DIRiniB.R")
	source("/home/hduser/qtlTest/mapreduce-direct-version/Prunedirect-V4/find_poB.R")
	source("/home/hduser/qtlTest/mapreduce-direct-version/Prunedirect-V4/mydirectB.R")
	source("/home/hduser/qtlTest/mapreduce-direct-version/Prunedirect-V4/replaceinfB.R")
	source("/home/hduser/qtlTest/mapreduce-direct-version/Prunedirect-V4/ObjectiveFB.R")
	source("/home/hduser/qtlTest/mapreduce-direct-version/Prunedirect-V4/probbcB.R")
	source("/home/hduser/qtlTest/mapreduce-direct-version/Prunedirect-V4/rssB.R")
	datafile <- read.cross("csv","/home/hduser/qtlTest/mapreduce-direct-version/Prunedirect-V4", "hyper2.csv",genotypes=c("BB","BA"));	
}

#--- Extract genotypes + phenotypes
 # choose pheo:datafile$pheno - it both contains genotypes and locations
 # choose only genotypes: datafile$geno[[18]][[1]]
 # choose locations: datafile$geno[[18]][[2]]
 
 #probtemp <- calc.genoprob(datafile, step=step, err=0.001); # 

 chosenChrom1 <- 6#4#6;
 chosenChrom2 <- 15#7#15;

 locations1 <- datafile$geno[[chosenChrom1]][[2]];
 numberOfMarkers1 <- length(locations1);

 locations2 <- datafile$geno[[chosenChrom2]][[2]];
 numberOfMarkers2 <- length(locations2);

 pheno <- datafile$pheno[[1]]; 
 mu <- mean(pheno);
 variance <- var(pheno);

 ideal_step <- 10^(-1);
 step3p_1 <- round(log((locations1[numberOfMarkers1] - locations1[1])/ideal_step,base=3))
 step3p_2 <- round(log((locations2[numberOfMarkers2] - locations2[1])/ideal_step,base=3))

 step1 <- (locations1[numberOfMarkers1] - locations1[1])*3^(-step3p_1);
 step2 <- (locations2[numberOfMarkers2] - locations2[1])*3^(-step3p_2);

 maxdeep <- min(step3p_1,step3p_2)+1;

prob1 <- probbc(datafile, chosenChrom1, step=step1, err=0.001);
prob2 <- probbc(datafile, chosenChrom2, step=step2, err=0.001);

meshSize1 <- (locations1[numberOfMarkers1] - locations1[1])/step1;
mesh1 <- locations1[1] + step1*(0:(meshSize1-1));

meshSize2 <- (locations2[numberOfMarkers2] - locations2[1])/step2; 
mesh2 <- locations2[1] + step2*(0:(meshSize2-1));
 
#location names loc10
j <- step1*(1:(meshSize1-1));
j <- as.character(j);
locnames1 <- 0;
locnames1[1] <- names(datafile$geno[[chosenChrom1]]$map)[1]
for (i in 2:meshSize1)
{
	locnames1[i] <- paste("loc",j[i-1],sep="");
}

j <- step2*(1:(meshSize2-1));
j <- as.character(j);
locnames2 <- 0;
locnames2[1] <- names(datafile$geno[[chosenChrom2]]$map)[1]
for (i in 2:meshSize2)
{
	locnames2[i] <- paste("loc",j[i-1],sep="");
}
#--------------------------------------------------------
probAA1 <- prob1[,locnames1,1];
probAB1 <- prob1[,locnames1,2];

probAA2 <- prob2[,locnames2,1];
probAB2 <- prob2[,locnames2,2];

#----------------------------------------------------------------
bounds=t(data.frame("b1"=c(mesh1[1],mesh1[length(mesh1)]),"b2"=c(mesh2[1],mesh2[length(mesh2)]))) # choóse the boundary of the chromosome : 2d
colnames(bounds)<-c("lower", "upper")
Problem=list("f"="ObjectiveF")

	
#ptm1 <- proc.time();
results <- mydirect(Problem = Problem, pheno ,probAA1, probAA2, variance = variance, mesh1 = mesh1, mesh2 = mesh2, bounds ,maxdeep = maxdeep, showits="none", verbose=T,permutation=0,per_minval=Inf,calculated_rss_saved=0)

per_minval<- 2.182178
#per_minval<- 4.4814
#------------------------
#Pre-calculate the values of RSS and save it for permuation, based on the min value found before
a=bounds[1,] # different index here! Careful
b=bounds[2,]
n <- length(pheno)
	
Thirds<-rep(NA, maxdeep)
Thirds[1] <- 1/3
 for (i in 2:maxdeep){
    Thirds[i]= (1/3)*Thirds[i-1];
 }	 
calculated_rss_saved <- matrix(NA,maxdeep+1,maxdeep+1)
for(l1 in 0:maxdeep)
{
	for(l2 in 0:maxdeep)
	{

	 if(l1==0){ 
		delta1 <- 1
	 }else{
 		delta1 <- Thirds[l1]
	 }
	  
	 if(l2==0){
 		delta2 <- 1
	 }else{
 		delta2 <- Thirds[l2]
	 }	

	 delta1<-(a[2]-a[1])*delta1
	 delta2<-(b[2]-b[1])*delta2

     x <- (delta1+delta2)/4
	 calculated_rss_saved[l1+1,l2+1] <- rss(n,x,per_minval)
	}
}	

#100 permuations 
result_permuted<-matrix(0,100,4)
for (i in 1:100)
{
	set.seed(i)
	pheno_permuted<-sample(pheno)
	print(i)	
	results2 <- mydirect(Problem = Problem, pheno_permuted ,probAA1, probAA2, variance = variance, mesh1 = mesh1, mesh2 = mesh2, bounds ,maxdeep = maxdeep, showits="none", verbose=F,permutation=1,per_minval=per_minval,calculated_rss_saved=calculated_rss_saved)
	result_permuted[i,1] <- results2$minval
	result_permuted[i,2] <- results2$final_point.xatmin[1]
	result_permuted[i,3] <- results2$final_point.xatmin[2]
	result_permuted[i,4] <- length(results2$Hfcncounter)
}

results$final_point.xatmin # final result   66.7 , 17.5 6,15  2.182178
write.table(out1, file = "C:/Users/Behrang/Documents/exh10p5.csv", sep = " ")


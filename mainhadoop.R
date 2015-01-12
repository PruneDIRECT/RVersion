#Main for PruneDIRECT-Hadoop 


#rm(list=ls(all=TRUE))
hadoop <- 1;
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

 chosenChrom1 <- 6;
 chosenChrom2 <- 15;

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
#---- DIRECT 

#bounds=t(data.frame("b1"=c(-2, 2),"b2"=c(-2,2))) # Test function
#bounds=t(data.frame("b1"=c(mesh[1],mesh[length(mesh)]))) # choóse the boundary of the chromosome : 1d

bounds=t(data.frame("b1"=c(mesh1[1],mesh1[length(mesh1)]),"b2"=c(mesh2[1],mesh2[length(mesh2)]))) # choóse the boundary of the chromosome : 2d
colnames(bounds)<-c("lower", "upper")
Problem=list("f"="ObjectiveF2")

#result <- mydirect(Problem = Problem, pheno ,probAA1, probAA2, variance = variance, mesh1 = mesh1, mesh2 = mesh2, bounds ,maxdeep = maxdeep, showits="none", verbose=FALSE)
#result$final_point.xatmin # final result 
result <- mydirect(Problem = Problem, pheno ,probAA1, probAA2, variance = variance, mesh1 = mesh1, mesh2 = mesh2, bounds ,maxdeep = maxdeep, showits="none", verbose=F,permutation=0,per_minval=Inf,calculated_rss_saved=0)

#--- Test with R/QTL

out1<-scantwo(datafile,chr=c(chosenChrom1,chosenChrom2),use=c("all.obs"),incl.markers=FALSE,assumeCondIndep=TRUE,method="hk",n.perm=10^5)
#summary(out1)
per_minval<- 2.182178
#per_minval<- 4.4814
#------------------------
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
#-----------------------------------------
#-----------------------------------------
#-----------------------------------------
#---Permutation testing
n_per <- 10^5# number of permutations

per_sum <- 0;

result_permuted <- 0;

perm <- function (n_per, pheno, Problem, bounds, starting_seeds, probAA1, probAA2, variance, mesh1, mesh2, result,maxdeep) {
per_status <- 0*c(1:n_per);

#set.seed(starting_seed);
min_val_str <- paste("seed", "min_val","xmin_1", "xmin_2", "\n",sep=" ");

for(i in 1:n_per)
{
	set.seed(as.integer(starting_seeds[i]));
	pheno_permuted <- sample(pheno);
	
	#result_permuted  <- mydirect(Problem = Problem, pheno_permuted ,probAA1, probAA2, variance = variance, mesh1 = mesh1, mesh2 = mesh2, bounds ,maxdeep = maxdeep, showits="none", verbose=FALSE); 
	result_permuted <- mydirect(Problem = Problem, pheno_permuted ,probAA1, probAA2, variance = variance, mesh1 = mesh1, mesh2 = mesh2, bounds ,maxdeep = maxdeep, showits="none", verbose=F,permutation=1,per_minval=per_minval,calculated_rss_saved=calculated_rss_saved)
	

	if (result$minval > result_permuted$minval)
	    per_status[i] <- 1;
  
    ss <- paste("still computing -> ",i,sep="");
    min_val_str <- paste(min_val_str, as.integer(starting_seeds[i]),result_permuted$minval,result_permuted$final_point.xatmin[1],result_permuted$final_point.xatmin[2],"\n",sep=" ");
    status(ss);       
}
per_status_list<-c(0:1);
per_status_list[1]<-length(which(per_status==0));
per_status_list[2]<-length(which(per_status==1)); 
##permutation_result <- per_sum/n_per;

Sys.setenv("HADOOP_CMD"="/hadoop/hadoop/bin/hadoop")
library(rhdfs)
hdfs.init()
map_id <- Sys.getenv("mapred_task_id")
map_id_file <- paste("DIRECT/results/f2/",map_id,sep="")
filehandler <- hdfs.file(map_id_file, "w")
hdfs.write(min_val_str, filehandler)
hdfs.close(filehandler)

return(per_status_list);
}

QTL <- function (input, output=NULL) {
  rmr.options(backend="hadoop")
  #rmr.options(backend="local")
  mapreduce(input=input, output=output, input.format="text", map=map, reduce=reduce, backend.parameters = list(hadoop = list(D = "mapred.task.timeout=3600000")) )

}

map <- function(.,lines) {

  seed.list <- strsplit(lines, '\\s')
  seed <- unlist(seed.list)
  w <-c(0:1);
  task <- perm(n_per, pheno, Problem, bounds,seed,probAA1,probAA2,variance, mesh1, mesh2, result,maxdeep);
  keyval(w,task)

}

reduce <- function(w, counts) {

 keyval(w, sum(counts))
}


Sys.setenv("HADOOP_HOME"="/hadoop/hadoop")

library(rmr2)
library(qtl)

options(warn=-1)
hdfs.root <- 'DIRECT'
#hdfs.root = '/home/behrang/mapreduce/R-Test-code'
hdfs.data <- file.path(hdfs.root, 'seeds76/f1')
hdfs.out <- file.path(hdfs.root, 'out-f1')

out <- QTL(hdfs.data, hdfs.out)

results <- from.dfs(out)

results.df <- as.data.frame(results, stringsAsFactors=F)

colnames(results.df) <- c('Opt', 'Rate')

results.df

#keyout <- perm(n_per, pheno, Problem, bounds );
#keyout

#proc.time() - ptm1


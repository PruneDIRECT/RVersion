ObjectiveF <- function(x,pheno_local, probAA1_local, probAA2_local, 
 				 variance_local , mesh1_local , mesh2_local){ # new function
	
	x1 <- x[1];
	x2 <- x[2];

	Px1 <- probAA1_local[,(x1-mesh1_local[1])/step1+1]-0.5;
	Px2 <- probAA2_local[,(x2-mesh2_local[1])/step2+1]-0.5;
	
	Px1x2 <- Px1*Px2; # If they are independet. 
	
	lm.x <- lm(pheno_local~Px1+Px2+Px1x2); 
 
	rss <- summary(lm.x)$sigma;
	sigma2 <- rss^2;
	f <-  -log(abs(-sigma2+variance_local)/variance_local)# make it positive
	return(f);
}

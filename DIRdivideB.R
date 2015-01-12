######################################################################################################
.DIRdivide<-function(a,b,Problem,pheno_local,
									probAA1_local, probAA2_local,variance_local , mesh1_local , mesh2_local,
									index,thirds,
 									lengths, fc, c, con, 
 								    feas_flags, p_fcncounter, szes,
 								    impcons, calltype,showits.flag=TRUE, itctr="", i="", prune,maxdeep,minval_prune,permutation,calculated_rss_saved,...){
 								    
 #%------------------------------------------------------------------%
 #% Function   :  DIRdivide                                          %
 #% Written by :  Dan Finkel                                         %
 #% Created on :  10/19/2002                                         %
 #% Purpose    :  Divides rectangle i that is passed in              %
 #% This function is modified by B. Mahjani on 1/12/2015             %
 #%------------------------------------------------------------------%
 fcncounter <- p_fcncounter
 #%-- 1. Determine which sides are the largest #########################
 li     <- lengths[,index]  
 biggy  <- min(li)
 ls     <- which(li==biggy)
 lssize <- length(ls)
 
 #%-- 2. Evaluate function in directions of biggest size #########################
 #%--    to determine which direction to make divisions
 oldc       <- c[,index,drop=FALSE]
 delta      <- thirds[biggy+1]   

 n <- length(pheno_local)
 

 # add or create new centers? c_old +/- delta*e_i --> 4 new centers(left, right, up, down)!   
 #newc_left  <- oldc[,matrix(1, nrow=1,ncol=lssize)];  ### ??? 
 #newc_right <- oldc(:,ones(1,lssize));                ### ???  
 
 newc_left  <- newc_right <- matrix(rep(oldc, lssize), ncol=lssize  ) 
 
 # initialize
 f_left <- con_left <- fflag_left   <- rep(NA, lssize)
 f_right <- con_right <- fflag_right    <- rep(NA, lssize) 
 
 # for each dimention (parameter) in ls create new centers : left and right from the old center
 for (i in 1:lssize){
     lsi  <- ls[i]
     
     # c_i +/- delta*e_i  , e_i has at the ith position 1, rest=0
     newc_left[lsi,i]  = newc_left[lsi,i] - delta;
     newc_right[lsi,i] = newc_right[lsi,i] + delta;
     
     # f(new_centers_left)
     func.left.list<- .CallObjFcn(Problem,pheno_local, probAA1_local, probAA2_local,variance_local , mesh1_local , mesh2_local,  point.x=newc_left[,i,drop=FALSE],a, b, impcons, calltype, ...)
 		f_left[i]<-  func.left.list$fcn_value
 		con_left[i]<-  func.left.list$con_value
 		fflag_left[i]<-  func.left.list$feas_flag
     
     # f(new_centers_right)
     func.right.list<- .CallObjFcn(Problem,pheno_local,probAA1_local, probAA2_local,variance_local , mesh1_local , mesh2_local, point.x=newc_right[,i,drop=FALSE],a, b, impcons, calltype, ...)
 		f_right[i]<-  func.right.list$fcn_value
 		con_right[i]<-  func.right.list$con_value
 		fflag_right[i]<-  func.right.list$feas_flag
     
 		# counter := add 2 	
     fcncounter <- fcncounter + 2
 }
   
 #%---- 2.1 Calculate w - min of function in new centers #############
 #w = [min(f_left, f_right)' ls]; 
 # like in DIRECTUserGuide:
 # w_i := min( f_right, f_left ) for i in 1:N
 # best function valueS ! in the largest space
 # it means for each dimention find separat min ! 
 w = apply(cbind(f_left, f_right), 1, min)
 
 
 #%-- 3. Sort w for division order #########################
 tmp.sort<-sort(w,index.return=TRUE)
 V<-tmp.sort$x; order<-tmp.sort$ix
 
 #%-- 4. Make divisions in order specified by order #########################
 for (i in 1:length(order) ){

    newleftindex  = p_fcncounter+2*(i-1)+1
    newrightindex = p_fcncounter+2*(i-1)+2
    #%-- 4.1 create new rectangles identical to the old one ########
    oldrect <- lengths[,index, drop= FALSE]
    lengths<- cbind(lengths, oldrect, oldrect)
 
    #%-- old, and new rectangles have been sliced in order(i) direction
    lengths[ls[order[i]],newleftindex ] <- lengths[ls[order[i]],index] + 1
    lengths[ls[order[i]],newrightindex] <-  lengths[ls[order[i]],index]  + 1;
    lengths[ls[order[i]],index]         <-  lengths[ls[order[i]],index]  + 1;
 
    #%-- add new columns to c
    c<-cbind(c, newc_left[,order[i]], newc_right[,order[i]] )
    colnames(c)[(ncol(c)-1):ncol(c)]<- c("left", "right")
	
#Prune for maxdeep

	if ((lengths[ls[order[i]],newleftindex ]+1) < maxdeep)
		{leftprune <- 0;}
	else 
		{leftprune <- 1;}
	
	if ((lengths[ls[order[i]],newrightindex ]+1) < maxdeep)
		{rightprune <- 0;}
	else 
		{rightprune <- 1;}
		
	if ((lengths[ls[order[i]],index ]+1) > maxdeep)
		{prune[index] <- 1;}
    
    #%-- add new values to f
    fc<- c(fc, f_left[order[i]], f_right[order[i]] )
    #Prune for prunedirect
		# Add prune
	l1 <- lengths[1,index]  
    l2 <- lengths[2,index]  

	 if(l1==0) delta1 <- 1
	 else delta1 <- thirds[l1]	
	  
	 if(l2==0) delta2 <- 1
	 else delta2 <- thirds[l2]	

	delta1<-(b[1]-a[1])*delta1
	delta2<-(b[2]-a[2])*delta2
	
   	x <- (delta1+delta2)/4
	#print(minval_prune)
    if(!permutation) calculated_rss <- rss(n,x,minval_prune) 
	else calculated_rss <- calculated_rss_saved[l1+1,l2+1]
	
	if (f_left[order[i]] >  calculated_rss)
		{leftprune <- 1;}
	
	if (f_right[order[i]] >  calculated_rss)
		{rightprune <- 1;}
		
	if (fc[index]  > calculated_rss)
		{prune[index] <- 1;}
		
	prune <- c(prune,leftprune,rightprune)
    #%-- add new values to con
    con<- c(con, con_left[order[i]], con_right[order[i]] )
   
    #%-- add new flag values to feas_flags
    feas_flags<- c(feas_flags, fflag_left[order[i]], fflag_right[order[i]] )
      
    #%-- 01/21/04 Dan Hack
    #%-- store sizes of each rectangle    ### sicher ????  A -vector or matrix
    #n = norm(A), a matrix  returns the largest singular value (!) of A, max(svd(A)$d).
    # if A vector:  norm(A)=sum(abs(A).^p)^(1/p), p=2
    # szes(1,newleftindex)  = 1/2*norm((1/3*ones(size(lengths,1),1)).^(lengths(:,newleftindex)));
    
 		#	   # if matrix
 		#		 tmp<-(1/3*rep(1,nrow(lengths)) )^ ( lengths[,newleftindex, drop=FALSE])
 		#		 max(svd(tmp)$d)
 		#		 # if vector
 		#		 tmp<-(1/3*rep(1,nrow(lengths)) )^ ( lengths[,newleftindex])
 		#	   sum(abs(tmp)^2)^(1/2)
 		#		 ## the same result :-) funny!!!! use matrix!
 		
 		tmp.szes.l<-(1/3*rep(1,nrow(lengths)) )^ ( lengths[,newleftindex, drop=FALSE])
 		tmp.szes.l<- 1/2*max(svd(tmp.szes.l)$d)	
 		
 		tmp.szes.r<-(1/3*rep(1,nrow(lengths)) )^ ( lengths[,newrightindex, drop=FALSE])
 		tmp.szes.r<- 1/2*max(svd(tmp.szes.r)$d)	
 		   
     szes<-c(szes, tmp.szes.l, tmp.szes.r )   ##
    names(szes)[(length(szes)-1):length(szes)]<- c("left", "right")
 } #end of for
   
 
 ## plot old and new centers ####
 if (showits.flag ){
 	my.col<- rep("black", ncol(c))
 	my.col[grep("left",colnames(c))]<- "blue"
 	my.col[grep("right",colnames(c))]<- "red"
 	
 	if (nrow(c)==1){
 		toPlot<- rbind(c, fc) 
 		#plot(as.data.frame(t(toPlot)), pch=20, type="p", xlim=c(0,1), col=my.col, main=paste("Iteration", itctr, "step", i  ) )
		plot(as.data.frame(t(toPlot)), pch=20, type="p")# xlim=c(0,1), col=my.col, main=paste("Iteration", itctr, "step", i  ) )
 		
 		#text(x=toPlot[1,],y=toPlot[2,], round(fc,1), pos=1, cex=0.3)		
 		#text(x=toPlot[1,],y=toPlot[2,], c(1:ncol(c)), pos=3, cex=0.3, col="#008080")	
 		#legend("topright", c("old center", "new left", "new right", "center's number"), fill=c("black", "blue", "red", "#008080"), cex=0.5)
 	} # end for 1 tuning parameter
 	
 	# for 2 tuning parameters
 	if (nrow(c)==2){
 	toPlot<- c
	plot(as.data.frame(t(toPlot)), cex=0.6,pch=20, type="p")#, xlim=c(0,1), ylim=c(0,1), col=my.col, main=paste("Iteration", itctr, "step", i  )  )
 	#plot(as.data.frame(t(toPlot)), pch=20, type="p", xlim=c(0,1), ylim=c(0,1), col=my.col, main=paste("Iteration", itctr, "step", i  )  )
 #	text(x=toPlot[1,],y=toPlot[2,], round(fc,1), pos=1, cex=0.3)		
 #	text(x=toPlot[1,],y=toPlot[2,], c(1:ncol(c)), pos=3, cex=0.3, col="#008080")	
 #	legend("topright", c("old center", "new left", "new right", "center's number"), fill=c("black", "blue", "red", "#008080"), cex=0.5)
 	} # end for 2 tuning parameters
 }
 ## end of plot old and new centers #### 		
 
 
 tmp.szes.ind<-(1/3*rep(1,nrow(lengths)) )^ ( lengths[,index, drop=FALSE])
 tmp.szes.ind<- 1/2*max(svd(tmp.szes.ind)$d)	
 szes[index] <- tmp.szes.ind
 pass = 1
 
# PruneDIRECT
	# Pruning for PruneDIRECT. 
	# minval known - should send it to this function
	#loop over 3 new boxes
		#Calculate rssx
		# save the coeff for prune
	#	if (rss(x) < f(x)) 
	#		save the index to prune permenantly
	#end loop

	#for(i in 1:3)
	#{
		#rss(n,x,fmin)
		
	#}
	#remove the pruned ones

 return(list(lengths=lengths,fc=fc,c=c,con=con, feas_flags=feas_flags,
 			 szes=szes,fcncounter=fcncounter,pass=pass,prune=prune))
 
 }
 
 ####################################################################################################

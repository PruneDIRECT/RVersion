###########################################################################################################      
 .DIRini<-function(Problem,pheno_local,
					probAA1_local, probAA2_local,variance_local , mesh1_local , mesh2_local,
					n,a,b,
 					param.names=c(1:length(a)),
                  #p_lengths,
                  c, fc, con, feas_flags, szes,
                  theglobalmin,
                  maxdeep,testflag,impcons,
				probAA1, probAA2, variance, mesh1, mesh2,
				...){
 
 #%------------------------------------------------------------------%
 #% Function:   DIRini                                               %
 #% Written by: Dan Finkel                                           %
 #% Created on: 10/19/2002                                           %
 #% Purpose   : Initialization of Direct                             %
 #%             to eliminate storing floating points                 %
 #%------------------------------------------------------------------%
 
 # DIRECT begins the optimization by transforming the domain of the problem into the unit
 # hyper-cube. The algorithm works in this normalized space, referring to the original space only when
 # making function calls. The center of this space is c1, and we begin by fnding f(c1).
 
 	#%-- start by calculating the thirds array
 	#%-- here we precalculate (1/3)^i which we will use frequently
 	l_thirds<-rep(NA, maxdeep)
 	l_thirds[1] <- 1/3
 	for (i in 2:maxdeep){
 	   l_thirds[i]= (1/3)*l_thirds[i-1];
 	}
 	#%--lengths=  length array will store # of slices in each dimension for
 	#%-- each rectangle. dimension will be rows; 
 	#%-- each rectangle will be a column
 	
 	#%-- first rectangle is the whole unit hyperrectangle
 	l_lengths <- matrix(0,n,1);
 	rownames(l_lengths)<-  param.names
 	
 	#%01/21/04 HACK
 	#%-- store size of hyperrectangle in vector szes
 	szes = 1
 	names(szes)<- "start"
 	
 	#%-- first element of c is the center of the unit hyperrectangle
 	#l_c(:,1) = matrix(1/2,n,1)
 	# erster Spalte
 	l_c <- matrix(1/2,n,1)
 	colnames(l_c)<-"start"
 	rownames(l_c)<-  param.names
 	
 	
 	#%-- Determine if there are constraints
 	calltype = .DetermineFcnType(Problem,impcons);
 		
 	#%-- first element of f is going to be the function evaluated
 	#%-- at the center of the unit hyper-rectangle.
 	#%om_point   = abs(b - a)*l_c(:,1)+ a;
 	#%l_fc(1)    = feval(f,om_point,varargin{:});
 	func.List<- .CallObjFcn(Problem,pheno_local,probAA1_local, probAA2_local,variance_local , mesh1_local , mesh2_local, point.x=l_c[,1, drop=FALSE],a, b, impcons, calltype, ...)
 	
 	## debugging  
 	# func.List<- .CallObjFcn(Problem, point.x=l_c[,1, drop=FALSE],a, b, impcons, calltype, fmin, fit.gp, muX, muY)
 	## end of debugging
 	
 	l_fc<- func.List$fcn_value
 	l_con<- func.List$con_value
 	l_feas_flags<- func.List$feas_flag
 	fcncounter = 1
 		
 	#%-- initialize minval and point.xatmin to be center of hyper-rectangle !  (NOT in the original intervals!!!!)
 	point.xatmin = l_c[,1, drop=FALSE]      
 	minval   = l_fc[1]
 	if (testflag == 1){
 	    if (theglobalmin != 0){
 	        perror = 100*(minval - theglobalmin)/abs(theglobalmin);
 	    }else{
 	        perror = 100*minval;
 	    }
 	}else {
 	   perror = 2
 	}
 	#%-- initialize History
 	History<-t(matrix(c(0,0,0)))  
   colnames(History)<- c("Iteration Nr", "Function Count", "f_min"  )
   prune <- 0;
 return(list( thirds=l_thirds, lengths=l_lengths,
 						 c=l_c, fc=l_fc, con=l_con, 
 					 	 feas_flags=l_feas_flags,
 		    		 minval=minval,point.xatmin=point.xatmin,perror=perror,
 		         History=History, szes=szes,
 		         fcncounter=fcncounter,calltype=calltype,prune=prune ))
 }
 
 ####################################################################################################
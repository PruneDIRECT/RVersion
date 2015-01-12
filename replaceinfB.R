##################################################################################################################
  .replaceinf<- function(lengths,c,fc,con,flags,pert){
 
 #%------------------------------------------------------------------%
 #% Function   :  replaceinf                                         %
 #% Written by :  Dan Finkel                                         %
 #% Created on :  06/09/2004                                         %
 #% Purpose    :  Assign R. Carter value to given point              %
 #%------------------------------------------------------------------%
 # 
 
 #%-- Initialize fcn_values to original values
 fcn_values <- fc
 
 #%-- Find the infeasible points
 infeas_points <- which(flags == 1)
 
 #%-- Find the feasible points
 feas_points   = which(flags == 0)
 
 #%-- Calculate the max. value found so far
 maxfc<- ifelse( length(feas_points)>0, 
 								max(fc[feas_points] + con[feas_points]),
 								max(fc + con) )
 
 if (length(infeas_points)>0){
 	for (i in 1:length(infeas_points) ){
 	
 	    if (length(feas_points)==0){
 	        #%-- no feasible points found yet
 	        found_points <-found_pointsf <- vector()
 	        index <- infeas_points[i];
 	    } else {
 	        index = infeas_points[i]
 	
 	        #%-- Initialize found points to be entire set
 	        found_points  <- c[,feas_points, drop=FALSE ]
 	        found_pointsf <- fc[feas_points] + con[feas_points]
 	
 	        #%-- Loop through each dimension, and find points who are close enough
 	        for (j in 1:nrow(lengths) ){
 	            neighbors <- which(abs(found_points[j,] - c[j,index]) <=  3^(-lengths[j,index]))
 	            if (length(neighbors)>0 ){
 	                found_points  <- found_points[,neighbors]
 	                found_pointsf <- found_pointsf[neighbors]
 	            } else{
 	                found_points <-found_pointsf <- vector()
 	                break
 	            } 
 	        } 
 	    } # end of if else
 	
 	    #%-- Assign Carter value to the point
 	    if (length(found_pointsf)>0) {
 	        #%-- assign to index the min. value found + a little bit more
 	        fstar <- min(found_pointsf);
 	        fcn_values[index] <- ifelse ((fstar != 0), 
 	        															fstar + pert*abs(fstar),
 	        															fstar + pert*1 )
 	    }else {
 	        fcn_values(index) = maxfc+1
 	        maxfc             = maxfc+1
 	    }
 	} # end of for
 } # end of if  (length(infeas_points)>0)
 return (fcn_values)
 } 

####################################################################################################
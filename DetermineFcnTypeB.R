####################################################################################################
  .DetermineFcnType<-function(Problem,impcons){
 #%------------------------------------------------------------------%
 #% Function   :  DetermineFcnType                                   %
 #% Written by :  Dan Finkel                                         %
 #% Created on :  06/25/2004                                         %
 #% Purpose    :  Determine how constraints are handled              %
 #%------------------------------------------------------------------%
 
 	retval = 0;
 	if ( !("constraint" %in% names(Problem)) & (!impcons) ){
 	    # 1. %-- No constraints at all  and  (no implicit constraint)
 	    retval = 1
 	}
 	if ("constraint" %in% names(Problem)){
 		#%-- There are explicit constraints. Next determine where
     #%-- they are called
     if (length(Problem$constraint)>0){ # slot constraints exists, look at their slots
 		        if (length(Problem$constraint[[1]]$func) == length(Problem$f))
 		            #%-- Constraint values may be returned from objective
 		            #%-- function. Investigate further
 		            if (Problem$constraint[[1]]$func == Problem$f ){  
 		                #%-- f returns constraint values
 		                retval = 2
 		            }else {
 		               # %-- f does not return constraint values
 		               retval = 3;
 		            } # end of There are explicit constraints
 		} else{ 
 			# thre is a EMPTY slot named constraints     
       if (impcons){
 		  	retval = 0;
 		  }else{
 		  	retval = 1;
 		  } 
 		}# end of no constraints
 	} #no slot named constraints    
 	
 	if (impcons){
 		    if (retval== 0 ){
 		        #%-- only implicit constraints
 		        retval = 4;
 		    }else{
 		        #%-- both types of constraints
 		        retval = 5;
 		    }
 	}
 	
 	return(retval)
 }
 
 
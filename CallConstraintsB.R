####################################################################################################
  .CallConstraints<- function(Problem,pheno_local,probAA1_local, probAA2_local,variance_local , mesh1_local , mesh2_local, x,a,b,...){  
 #%------------------------------------------------------------------%
 #% Function   :  CallConstraints                                    %
 #% Written by :  Dan Finkel                                         %    
 #% Created on :  06/07/2004                                         %
 #% Purpose    :  Evaluate Constraints at pointed specified          %
 #%------------------------------------------------------------------%

 
 	#%-- Scale variable back to original space
 	point = abs(b - a)*x+ a;
 	
 	ret_value = 0;
 	if ( ("constraint" %in% names(Problem))  ){
 	    if (length(Problem$constraint)>0){ # slot constraints exists, look at their slots
 	        for (i in 1:Problem$numconstraints){
 	            if (length(Problem$constraint[[i]]$func) == length(Problem$f)){
 	                if ( Problem$constraint[[i]]$func == Problem$f ){                    
 	                    #%-- Dont call constraint; value was returned in obj fcn
 	                    con_value = 0;
 	                }else{
 	                    con_value = eval(parse(text=Problem$constraint[[i]]$func))(point, ...)
 	                } 
 	            }else{
 	               con_value = eval(parse(text=Problem$constraint[[i]]$func))(point, ...)
 	            }
 	            if (con_value > 0) {
 	                #%-- Infeasible, punish with associated pen. param
 	                ret_value = ret_value + con_value * Problem$constraint[[i]]$penalty
 	            } 
 	        } # end of for
 	   } # end of (length(Problem$constraint)>0)
 	}
 	return (ret_value)
 }
 
 ###########################################################################################################
 
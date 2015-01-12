###########################################################################################################
  .CallObjFcn<-function(Problem,pheno_local, probAA1_local, probAA2_local, variance_local , mesh1_local , mesh2_local, point.x,a,b,impcon,calltype,...){
 #%------------------------------------------------------------------%
 #% Function   :  CallObjFcn                                         %
 #% Written by :  Dan Finkel                                         %   
 #% Created on :  06/07/2004                                         %
 #% Purpose    :  Evaluate ObjFcn at pointed specified               %
 #%------------------------------------------------------------------%
 
 ## point.x = vector of values for tuning parametr(s)  
 ## in arguments of Problem function: point - vector of tuning parameters at the first place
 
 	con_value = 0;
 	feas_flag = 0;
 	
 	#%-- Scale variable back to original space
 	point = abs(b - a)*point.x+ a
 	
 	if (calltype == 1){
 	    #%-- No constraints at all
 	    # find the functions value at 'point'
 	    fcn_value = eval(parse(text=Problem$f))(point,pheno_local,probAA1_local, probAA2_local,variance_local , mesh1_local , mesh2_local, ...)
 	    
 	    ## debug
 	    # fcn_value = eval(parse(text=Problem$f))(point, fmin, fit.gp, x.svm,y.svm)
 	    # end of debug
 	}
 	if (calltype == 2){
 	    #%-- f  and   all constraints
 	    tmp.list2<- eval(parse(text=Problem$f))(point,pheno_local,probAA1_local, probAA2_local,variance_local , mesh1_local , mesh2_local, ...)
 	    fcn_value<-tmp.list2$fcn_value
 	    cons<-tmp.list2$cons
 	   # [fcn_value, cons] = feval(Problem.f,point,varargin{:});
 	    for (i in 1:length(cons)){
 	        if (cons > 0){
 	         con_value <- con_value + Problem$constraint[[i]]$penalty * cons(i);
 	        }
 	    }
 	}
 	if (calltype == 3){    
 	    #%-- f returns no constraint values
 	    fcn_value <- eval(parse(text=Problem$f))(point,pheno_local,probAA1_local, probAA2_local,variance_local , mesh1_local , mesh2_local, ...)
 	    con_value <- .CallConstraints(Problem,pheno_local,probAA1_local, probAA2_local,variance_local , mesh1_local , mesh2_local, point.x,a,b,...);
 	}
 	if (calltype == 4){  
 	    #%-- f returns feas flag
 	    tmp.list4<- eval(parse(text=Problem$f))(point,pheno_local,probAA1_local, probAA2_local,variance_local , mesh1_local , mesh2_local, ...)
 	    fcn_value<-tmp.list4$fcn_value
 	    feas_flag<-tmp.list4$feas_flag
 	}
 	if (calltype == 5){
 	    #%-- f returns feas flags, and there are constraints
 	    tmp.list5<- eval(parse(text=Problem$f))(point,pheno_local,probAA1_local, probAA2_local,variance_local , mesh1_local , mesh2_local, ...)
 	    fcn_value<-tmp.list5$fcn_value
 	    feas_flag<-tmp.list5$feas_flag
 	     con_value <- .CallConstraints(Problem,pheno_local,probAA1_local, probAA2_local,variance_local , mesh1_local , mesh2_local, point.x,a,b,...);
 	}
 	
 	if (feas_flag == 1){
 		fcn_value = 10^9
 	  con_value = 0
 	}
 	return(data.frame("fcn_value"=fcn_value, "con_value"=con_value, "feas_flag"=feas_flag))
 }
 
 ##################################################################################################################
####################################################################################################  
 .find_po<-function(fc,lengths,minval,ep,szes, prune){
 
 #%--------------------------------------------------------------------%
 #% Function   :  find_po                                              %
 #% Written by :  Dan Finkel                                           %
 #% Created on :  10/19/2002                                           %
 #% Purpose    :  Return list of potentially optimal hyper-rectangles  %
 #%--------------------------------------------------------------------%
 
 #%-- 1. Find all rects on hub
 # diff_szes = sum(lengths,1);
 diff_szes = colSums(lengths) # col sum or row sums? nicht sicher ??????
 tmp_max = max(diff_szes)
 j=1
 hull<-vector()   ##??? nicht sicher !!!!
 sum_lengths = colSums(lengths)
 
 for (i in 1:(tmp_max+1)){
  if(!prune[1]) #not prune	
   {  tmp_idx <- which(sum_lengths == (i-1))
     if(length(tmp_idx)>0){
 	    tmp_n <- min(fc[tmp_idx])
 	    hullidx <- which.min(fc[tmp_idx]) 
 	    if (length(hullidx) > 0 ){
 	        hull[j] <- tmp_idx[hullidx]
 	        j=j+1;
 	        #%-- 1.5 Check for ties
 	        ties <- which(abs(fc[tmp_idx]- tmp_n) <= 1e-13)
 	        if (length(ties) > 1){
 	            mod_ties <- which(tmp_idx[ties] != hull[j-1])
 	            hull <- c(hull, tmp_idx[ties[mod_ties]])
 	            j <- length(hull)+1;
 	        } # end of the if length(ties) > 1
 	    } # end of the  if length(hullidx) > 0 
     } # end of if length(tmp_idx)>0
	}#	if prune
 } # end of for 
 
 
 #%-- 2. Compute lb and ub for rects on hub
 lbound <- .calc_lbound(lengths,fc,hull,szes)
 ubound <- .calc_ubound(lengths,fc,hull,szes)
 
 #%-- 3. Find indeces of hull who satisfy
 #%--    1st condition
 maybe_po <- which(lbound-ubound <= 0)
 
 #%-- 4. Find indeces of hull who satisfy
 #%--    2nd condition
 t_len  <- length(hull[maybe_po])
 if (minval != 0){
     po = which( ( (minval-fc[hull[maybe_po]])/abs(minval) +
                    szes[hull[maybe_po]] * ubound[maybe_po]/abs(minval)  ) >= ep)
 }else {
     po = which (( fc[hull[maybe_po]] - 
                   szes[hull[maybe_po]] * ubound[maybe_po] ) <= 0)
 } 
 
 final_pos  <- hull[maybe_po[po]]
 
 
 rects <- rbind(final_pos, "szes"=szes[final_pos])  
 return(rects)
 }
 
 ####################################################################################################
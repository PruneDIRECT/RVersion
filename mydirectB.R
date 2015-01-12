 #%------------------------------------------------------------------%
 #% Function   :  mydirect                                           %
 #% Written by :  Dan Finkel                                         %
 #% Created on :  10/19/2002                                         %
 #% Purpose    : 						     %
 #% This function is modified by B. Mahjani on 1/12/2015             %
 #%------------------------------------------------------------------%
 
# c changed to c_value !!! 
mydirect<-function(Problem,pheno_local,
				 probAA1_local, probAA2_local, variance_local , mesh1_local , mesh2_local,
                 bounds, maxdeep,
                 # options
                 maxits =     20000,         #% maximum of iterations
                 maxevals =   500000,       #% maximum # of function evaluations
                 #maxdeep =    7,#1000,        #% maximum number of side divisions
                 testflag =   0,          #%  the optimal value is unknown
                 globalmin =  0,          #% minimum value of function
                 ep =         1e-4,       #% global/local weight parameter.
                 tol =        0.01,       #% allowable relative error if f_reach is set
                 showits =    c("none", "final", "all"), #%  plot iteration stats: none, final iteration, all iterations 
                 verbose = 		TRUE, 			#% print  iteration stats: none, final itertation, all iterations 
                 impcons =    0,          #% flag for using implicit constraint handling
                 pert =       1e-6,       #% pertubation for implicit constraint handling
                 # maxflag =   0,          #% set to 1 for max problems, 0 for min problems
                 # sizeconst = 0.5,        #% constant on rectangle size function
                 # distance =  1,          #% 1/0 for distance/volume measure of size
                 # minlength = 1e-4,       #% stop if best rectangle has all sides 1ess than this
                 #minevals =  0          #% but must evaluate at least this many points 
                 
                 # plot parameter
                 pdf.name=NULL, 
                 pdf.width=12, pdf.height=12,
                 my.mfrow=c(1,1), 
				permutation,
 				per_minval,
				calculated_rss_saved,
                 ...){

##%-- Initialize the variables --------------------------------------%
minval_prune<-per_minval;

lengths<-c_value <- fc <- prune <- vector()
con <- szes <- feas_flags <- vector()
om_lower     <- bounds[,"lower", drop=FALSE]
om_upper     <- bounds[,"upper", drop=FALSE]
fcncounter   <- 0
perror       <- 0
itctr        <- 1
done         <- 0

n            <-nrow(bounds)

#% Determine option values
theglobalmin = globalmin

minc_pruned <- 0
minfc_pruned <- 0


# %-- New 06/08/2004 Pre-allocate memory for storage vectors
if (testflag == 0){
    lengths    = matrix(0,n,c(maxevals + floor(.10*maxevals)))
    c_value          = lengths;
    fc         = matrix(0,1,c(maxevals + floor(.10*maxevals))) 
	prune      = matrix(0,1,c(maxevals + floor(.10*maxevals))) 
    szes       = fc
    con        = fc
    feas_flags = fc
}

#%-- Call DIRini ---------------------------------------------------%
# define the multy dim hypercube with center point c, fc=f(c) 

DIRini.list<-.DIRini(Problem,pheno_local,
		probAA1_local, probAA2_local,variance_local , mesh1_local , mesh2_local,
		n, a=bounds[,"lower"], b=bounds[, "upper"],
        param.names= rownames(bounds),
        c=c_value,fc=fc, con=con, 
        feas_flags=feas_flags,  szes=szes,
        theglobalmin, maxdeep, testflag, impcons
		,...)
          
thirds <- DIRini.list$thirds
#lengths =  length array! will store number of slices in each dimension for 
#each rectangle. dimension will be rows; 
#each rectangle will be a column
lengths <-DIRini.list$lengths
minfc_pruned <- c_value <-DIRini.list$c 
minf_pruned <- fc <-DIRini.list$fc

con <-DIRini.list$con
feas_flags <-DIRini.list$feas_flags
minval <-DIRini.list$minval
point.xatmin <-DIRini.list$point.xatmin
perror <-DIRini.list$perror
History <-DIRini.list$History
szes <-DIRini.list$szes                 # number of regions
fcncounter <-DIRini.list$fcncounter
calltype <-DIRini.list$calltype
prune <- DIRini.list$prune;
 
ret_minval = minval
ret_point.xatmin = point.xatmin

Hfcncounter <- 0

if(showits !="none" & !is.null(pdf.name)) { 
	pdf(pdf.name, pdf.width, pdf.height)
}
if(showits != "none")
	par(mfrow=my.mfrow)
#print("hello")
#%-- MAIN LOOP -----------------------------------------------------%
minval = fc[1] + con[1]
minfc_pruned <- fc[1] + con[1]
mintemp <-minval

#print(paste("---",minval_prune,permutation,sep=" "))
if(!permutation) minval_prune<- mintemp
#print(paste("2---",minval_prune,permutation,sep=" "))
while (perror > tol){
	mintemp <- minval # for prunning
   #%-- Create list S of potentially optimal hyper-rectangles
   S <- .find_po(fc=fc+con,
       				 lengths= lengths,
       				 minval=minval, ep=ep, szes=szes, prune=prune)

	# if we don't find potentially hyper-rectanges --> break!
	if (ncol(S)==0) {print("---> S=0 <---");break;}
	#S_sorted <- sort(S[1,],index.return = TRUE)
	#S <- S[,S_sorted$ix]
	o<-order(S[1,])
	S <- rbind(S[1,o],S[2,o])
	
	#print(ncol(S))
   #%-- Loop through the potentially optimal hyper-rectangles ------%
   #%-- and divide -------------------------------------------------%
   for (i in 1:ncol(S)){
      # plot options: don't plot if not requested
     if ( (showits =="none") ) { 
     	showits.flag <-  FALSE
     } else {
	     # plot  last iteration's step  if requested
	     if ((showits =="final")) showits.flag<- ifelse ( (i == ncol(S)), TRUE, FALSE )
	     # plot all iterations' steps  if requested
	     if ( (showits =="all") )showits.flag <-  TRUE
     }
	#if(max(as.vector(lengths[,S[1,i]]+1))<maxdeep)
	
	
	
    if(!prune[S[1,i]])
	{ 

#print(minval_prune)
	if(!permutation) minval_prune<- mintemp
	 tmp.list.divide<- .DIRdivide(a=bounds[,1], b=bounds[,2],Problem=Problem,pheno_local, probAA1_local, probAA2_local,variance_local , mesh1_local , mesh2_local,
														index=S[1,i], thirds=thirds, lengths=lengths,
   													fc=fc,c=c_value,con=con, feas_flags=feas_flags,
          									p_fcncounter=fcncounter, szes=szes,
          									impcons=impcons, calltype=calltype, showits.flag=showits.flag, itctr=itctr, i=i,prune = prune, maxdeep = maxdeep,minval_prune=minval_prune,permutation,calculated_rss_saved,...)
  
	# if (max(as.vector(tmp.list.divide$lengths)) < maxdeep ){
	  lengths <- tmp.list.divide$lengths
	  fc  <- tmp.list.divide$fc
	  c_value  <- tmp.list.divide$c
	  con <- tmp.list.divide$con
	  feas_flags <- tmp.list.divide$feas_flags
	  szes <-  tmp.list.divide$szes
	  fcncounter <-  tmp.list.divide$fcncounter
	  success <- tmp.list.divide$pass  
   	  prune <- tmp.list.divide$prune
	}
	else{ # pruning the maxdeep boxes 
		 if (verbose){ print("Exceeded Max depth. Pruning step ");}
	     if(1){ # if the min was pruned, we should save it
			if(mintemp >= (fc[S[1,i]] + con[S[1,i]]))# find the min of the pruned boxes
			{
				minfc_pruned <- (fc[S[1,i]] + con[S[1,i]])
				minc_pruned <- (om_upper - om_lower)*c_value[,S[1,i]] + om_lower
				mintemp <- minfc_pruned
			}
 			lengths <- lengths[,-S[1,i]]
			fc  <- fc[-S[1,i]]
			c_value  <- c_value[,-S[1,i]]
			con <- con[-S[1,i]]
			feas_flags <- feas_flags[-S[1,i]]
			szes <-  szes[-S[1,i]]
			fcncounter <-  fcncounter - 1
			prune <- prune[-S[1,i]]
			S[1,] <- S[1,] - 1;
		 }
	}

	
  }
 
 
  
   #%-- update minval, point.xatmin --------------------------------------%
   # [minval,fminindex] =  min(fc(1:fcncounter)+con(1:fcncounter)); sicher ????
   minval<- min(fc + con )
   fminindex<- which.min(fc + con)

 if (minval < minfc_pruned)
 {
   fminindex<- which.min(fc + con)
   penminval = minval + con[fminindex]
   point.xatmin = (om_upper - om_lower)*c_value[,fminindex] + om_lower
 }
else
{
	minval <- minfc_pruned
	penminval <- minval + con[fminindex] # wrong here
	point.xatmin <- minc_pruned
}
   if((minval < per_minval)&&(permutation))
	{
		per_minval <- minval;
		permutation <- 0;
	}
   if ( (con[fminindex] > 0)|(feas_flags[fminindex] != 0) ){ # not sure about this
       #%--- new minval is infeasible, don't do anything
   }else {
       #%--- update return values
       ret_minval <- minval;
       ret_point.xatmin <- point.xatmin;
   } 
 
   #%--see if we are done ------------------------------------------%
   if (testflag == 1){
      #%-- Calculate error if globalmin known
      perror<- ifelse ((theglobalmin != 0), 
      									100*(minval - theglobalmin)/abs(theglobalmin),
      									100*minval )
   }else{
      #%-- Have we exceeded the maxits?
      if (itctr >= maxits){
         if (verbose) print("Exceeded max iterations. Increase maxits")
         done <- 1
      }
      #%-- Have we exceeded the maxevals?
      if (fcncounter > maxevals){
         if (verbose) print("Exceeded max fcn evals. Increase maxevals")
         done <- 1
      }
      if (done == 1)
         perror = -1
		
   } # end of if else
   
   #if (max(as.vector(lengths)) >= maxdeep ){
      #%-- We've exceeded the max depth
      #if (verbose) print("Exceeded Max depth. Increse maxdeep")
      #perror = -1
   #}
      if (min(as.vector(lengths)) >= (0.5*maxdeep) ){ ###############Stoping criteria for DIRECT
      #%-- We've exceeded the max depth
		L_History <- length(History);
		if (abs(History[L_History]-History[L_History-1])<0.000001)
		{
			if (verbose) print("cover most of the space")
			perror = -1
		}
      }
   
   #%-- Store History
   History<-rbind(History, 
                 c(itctr, fcncounter, minval))
  
  #%-- New, 06/09/2004
  #%-- Call replaceinf if impcons flag is set to 1
  if (impcons == 1) {
    fc <- .replaceinf(lengths=lengths,c=c_value,fc=fc,con=con,
                    flags=feas_flags, pert=pert)
  }

  #%-- show iteration stats
  if (verbose)  print(paste("Iter:", itctr, "f_min:", minval, "fn evals:","X1",ret_point.xatmin[1],"X2",ret_point.xatmin[2], fcncounter, sep="   "))
   Hfcncounter <- c(Hfcncounter, fcncounter)
  itctr  = itctr + 1
#if (perror > tol) print(",,,,,");
} # end  of while (perror > tol)

if(showits !="none" & !is.null(pdf.name)) dev.off()

#%-- Return values      #################
#%-- return x*
final_point.xatmin <- ret_point.xatmin;

#%-- chop off (abschneiden) 1st row of History
History<-History[-1,]

return (list(final_point.xatmin=final_point.xatmin,
							minval =minval, 
							c=c_value, fc=fc, 
							History=History, Hfcncounter= Hfcncounter))
}

###########################################################################################################
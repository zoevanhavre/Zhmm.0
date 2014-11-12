#' PostProcessing function univ
#'
#' This function draws samples from a Wishart dist
#' @param v and s
#' @keywords Wishart
#' @export
#' @examples
#' #nope

Zswitch_hmm<-function(GrunK0, PropMin=0.1 ){
			K<-dim(GrunK0$q0)[2]
			
			 # Pick Reference = Max log Likelihood
			wml<-which.max(GrunK0$MAP)
			Zref<-factor(GrunK0$Z[wml,], levels=1:K,ordered=FALSE) 
			FinalOrderChoice<-order(GrunK0$q0[wml,], decreasing=TRUE)		
			non0ref<-FinalOrderChoice[1:sum(table(Zref)>0)]
			refComp<-c(GrunK0$q0[wml,non0ref], GrunK0$Mu[wml,non0ref])  # use only q0 and mu
			Zref<- factor(Zref)

			# RENAME TO neat numbers
			levels(Zref)<- c(1:K)[order(FinalOrderChoice[1:sum(table(Zref)>0)])]
	
			# storage dataframes:
		
			numK0now<-sum(table(Zref )>0) 
			numPar<-length(unlist(lapply(1:dim(GrunK0$Z)[1], rep, numK0now)))
			AllPars<-data.frame(matrix(0, ncol=4+numK0now, nrow=numPar) )
			AllPars[,1]<-unlist(lapply(1:dim(GrunK0$Z)[1], rep, numK0now))
			Zfixed<-GrunK0$Z
			#for each iteration
			for(.iter in 1:dim(GrunK0$q0)[1]){
				
				#Store current states
				Znow<-factor(GrunK0$Z[.iter,])    
				
				#identify potential candidate switches:
				CandiCells<-table(Znow , Zref)/apply(table(Znow , Zref), 1, sum)>PropMin
				getCandi<-function(x) { as.numeric(as.character(row.names(as.matrix(x))[x])) }
				ListCandi<- apply(CandiCells, 1, getCandi)
				
				# R stuff to make sure it deals with inputs correctly
				if(class(ListCandi)=='matrix'){
				ListCandi<-split(ListCandi, rep(1:ncol(ListCandi), each = nrow(ListCandi)))
				Candies<-expand.grid(ListCandi)  # each row is a labelling
				names(Candies)<-row.names(CandiCells)   # RAAAAH
				} else if (class(ListCandi)=='numeric'){
				Candies<-ListCandi
				} else {
				Candies<-expand.grid(ListCandi)  # each row is a labelling
				names(Candies)<-row.names(CandiCells)   # RAAAAH
				}

				namesCandies<-names(Candies)
				# Catch if no appropriate seperation available at all, check all permutations of potential groups. Only occurs in really bad models with bad convergence.				
				done<-0
				if(class(Candies)=='data.frame'){
				if  ( max(sapply(apply(Candies, 1, unique), length))<length(row.names(CandiCells))){
					#Candies<- permutations(K)
					Candies<- matrix( as.numeric( row.names(CandiCells))[permutations(numK0now,numK0now )], ncol=numK0now)	
					 colnames(Candies)<-as.numeric(names(as.data.frame(CandiCells)))
### ISSUE HERE
						#MinusRefPars<-function(x) 	{	flp<- as.numeric(  row.names(CandiCells))[unlist(Candies[x,])]
					#		flp<-na.omit(flp)
				#		if(length(unique(flp))<length(flp)) { Inf
				#		} else {sum(abs( (refComp	-  c(GrunK0$q0[.iter,flp], GrunK0$Mu[.iter,flp]))/refComp))	}}
			
						MinusRefPars_catch<-function(x) 	{ flp<- Candies[x,]
						if(length(unique(flp))<length(flp)) { Inf
						} else {sum(abs( (refComp	-  c(GrunK0$q0[.iter,flp], GrunK0$Mu[.iter,flp]))/refComp))	}}

						BestOne<-which.min( sapply(1:dim(Candies)[1] , MinusRefPars_catch))  # find the best perm out of options
						BestOne<-Candies[BestOne,]
					
						#if(is.null(names(BestOne))) {names(BestOne)<-namesCandies}
						flp<-as.numeric(BestOne)[as.numeric(names(BestOne))]
			
						# Allocations
						Znew<-Znow;				levels(Znew)<-as.numeric(names(BestOne))
						Zfixed[.iter,]<-as.numeric(as.character(Znew))
						# Parameters
						#flp<- as.numeric( row.names(CandiCells)[unlist(BestOne)])
						swQ<-nowQ<-matrix(GrunK0$Q[.iter,], nrow=K, byrow=T)[flp, flp]					
						combinePars<- cbind(.iter, 1:numK0now,  GrunK0$q0[.iter,flp],GrunK0$Mu[.iter,flp],swQ)#[order(as.numeric(BestOne)), decreasing=FALSE),]
						AllPars[AllPars[1]==.iter,]<- combinePars
						done<-1

					}}
					#else{
					#	if  (length(unique(Candies))<length(row.names(CandiCells))){
					#Candies<- permutations(K)}
					#} 
			if (done==0){
			MinusRefPars<-function(x) 	{	flp<- as.numeric(  row.names(CandiCells)[unlist(Candies[x,])])
								flp<-na.omit(flp)
							if(length(unique(flp))<length(flp)) { Inf
							} else {sum(abs( (refComp	-  c(GrunK0$q0[.iter,flp], GrunK0$Mu[.iter,flp]))/refComp))	}}
											
				if( sum( apply(CandiCells, 1, sum)) >  dim(CandiCells)[1] ){
					BestOne<-which.min( sapply(1:dim(Candies)[1] , MinusRefPars))  # find the best perm out of options
					BestOne<-Candies[BestOne,]
					} else {BestOne<- Candies }   # chose this one if no comparing needed
				if(is.null(names(BestOne))) {names(BestOne)<-namesCandies}
	
				# Allocations
				Znew<-Znow; levels(Znew)<-as.numeric(BestOne)
				Zfixed[.iter,]<-as.numeric(as.character(Znew))
				# OK WORKS
					# Parameters
				## FIX BELOW  "  la longueur des donn'ees [21] n'est pas un diviseur ni un multiple du nombre de colonnes [5]"
flp<- as.numeric( row.names(CandiCells)[unlist(BestOne)])
swQ<-nowQ<-matrix(GrunK0$Q[.iter,], nrow=K, byrow=T)[flp, flp]					
combinePars<- cbind(.iter, as.numeric(BestOne),  GrunK0$q0[.iter,flp],GrunK0$Mu[.iter,flp],swQ)#[order(as.numeric(BestOne)), decreasing=FALSE),]
AllPars[AllPars[1]==.iter,]<- combinePars
			}
		}
	
			
			# sumarise Z (find max)
			maxZ<-function (x)  as.numeric(names(which.max(table( x ))))
			Zhat<- factor( apply(Zfixed, 2,maxZ))
			#levels(Zhat)<- levels(Zhat)<-as.character(BestOne)
			return(list('Pars'=AllPars, 'Z'=Zfixed))
			}


		
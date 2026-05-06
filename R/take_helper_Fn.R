dCMP <-	function( x, lambda, mu, nu, log=TRUE, tol=0.01, iter.max=200 ){
  if(!missing(mu) & !missing(lambda)) stop("'mu' and 'lambda' both specified")
  if(missing(mu) & !missing(lambda)) loglike = x*log(lambda) - nu*lfactorial(x) - compute_CMP_constant(Lambda=lambda, Nu=nu, Mu=NA, Tol=tol, Max=iter.max, Log=TRUE, Type="Z")
  if(!missing(mu) & missing(lambda)) loglike = nu*x*log(mu) - nu*lfactorial(x) - compute_CMP_constant(Lambda=NA, Nu=nu, Mu=mu, Tol=tol, Max=iter.max, Log=TRUE, Type="S")
  if(missing(mu) & missing(lambda)) stop("Neither 'mu' or 'lambda' is specified")
  if(log==TRUE) return( loglike )
  if(log==FALSE) return( exp(loglike) )
}
rCMP <-	function( n, lambda, mu, nu, tol=0.01, x_max=200 ){
  loglike_x = rep(NA, x_max+1)
  for( x in 0:x_max ){
    if(!missing(mu) & !missing(lambda)) stop("'mu' and 'lambda' both specified")
    if(missing(mu) & !missing(lambda)) loglike_x[x+1] = dCMP( x=x, lambda=lambda, mu=mu, nu=nu, log=TRUE, tol=tol, iter.max=x_max)
    if(!missing(mu) & missing(lambda)) loglike_x[x+1] = dCMP( x=x, lambda=lambda, mu=mu, nu=nu, log=TRUE, tol=tol, iter.max=x_max)
    if(missing(mu) & missing(lambda)) stop("Neither 'mu' or 'lambda' is specified")
  }
  n = sample( x=0:x_max, size=n, replace=TRUE, prob=exp(loglike_x))
  return(n)
}
compute_CMP_constant <-function(Lambda, Nu, Mu, Tol, Max, Log=TRUE, Type="Z"){       # Mu = Lambda^(1/Nu)
  if( (!is.na(Lambda) & Lambda > 10^Nu) | (!is.na(Mu) & Mu^Nu > 10^Nu) ){
    if(Type=="Z"){
      #Const = exp(Nu*Lambda^(1/Nu)) / ( Lambda^((Nu-1)/(2*Nu)) * (2*pi)^((Nu-1)/2) * sqrt(Nu) )
      ln_Const = Nu*Lambda^(1/Nu) - ((Nu-1)/(2*Nu))*log(Lambda) - ((Nu-1)/2)*log(2*pi) - (1/2)*log(Nu)
    }
    if(Type=="S"){
      #Const = exp(Nu*Mu) / ( Mu^((Nu-1)/(2)) * (2*pi)^((Nu-1)/2) * sqrt(Nu) )
      ln_Const = Nu*Mu - ((Nu-1)/(2))*log(Mu) - ((Nu-1)/2)*log(2*pi) - (1/2)*log(Nu)
    }
  }else{
    Const = rep(0,Max+1)
    Index = 1
    Const[Index] = 1    # Z(0)
    while( Const[Index]/Const[1] > Tol ){
      if(Type=="Z") Const[Index+1] = Const[Index] * ( Lambda / Index^Nu )  # Z(Index)
      if(Type=="S") Const[Index+1] = Const[Index] * ( Mu / Index )^Nu  # Z(Index)
      Index = Index + 1
    }
    ln_Const = log(sum(Const))
  }
  if(Log==TRUE) Return = ln_Const
  if(Log==FALSE) Return = exp(ln_Const)
  return(Return)
}
	# AO EDIT April 2026 
  # In December 2025, we received new data for all months between March 2019 and October 2025 
	# Data have been updated into renamed files to capture the last date of nesting data included, "XX_nests_October2025.csv"
  # Buru logic added!
data.extract <- function(location, year.begin, year.end, season.begin = year.begin, season.end = year.end, file.path = NULL) {
	  if (is.null(season.begin)) season.begin <- year.begin
	  if (is.null(season.end)) season.end <- year.end
	  
	  if (location == "JM") {
	    data.0 <- read.csv(paste0(file.path, "JM_nests_October2025.csv"))
	    data.0 %>% 
	      select(Year_begin, Month_begin, JM_Nests) %>%
	      mutate(Nests = JM_Nests) -> data.0
	    
	  } else if (location == "W") {
	    data.0 <- read.csv(paste0(file.path, "W_nests_October2025.csv"))
	    data.0 %>% 
	      select(Year_begin, Month_begin, W_Nests) %>%
	      mutate(Nests = W_Nests) -> data.0
	    
	  } else if (location == "Buru") {
	    data.0 <- read.csv(paste0(file.path, "Buru_Monthly_Nests.csv"))
	    data.0 %>% 
	      # Use rowSums with na.rm = TRUE for  a total even if one beach is missing
	      # Then replace 0 with a very small number (e.g., 0.1) or NA to avoid log(-Inf)
	      mutate(Nests = Nests_Wamlana + Nests_Waenibe + Nests_Waspait) %>%
	      mutate(Nests = ifelse(Nests <= 0, 0.01, Nests)) %>% 
	      select(Year_begin = Year, Month_begin = Month, Nests) -> data.0
	  }

	# create regularly spaced time series:
	data.2 <- data.frame(Year = rep(min(data.0$Year_begin, na.rm = T):max(data.0$Year_begin, na.rm = T),each = 12),
											 Month_begin = rep(1:12, max(data.0$Year_begin, na.rm = T) - min(data.0$Year_begin, na.rm = T) + 1)) %>%
		mutate(begin_date = as.Date(paste(Year, Month_begin, '01', sep = "-"), format = "%Y-%m-%d"), Frac.Year = Year + (Month_begin-0.5)/12) %>%
		select(Year, Month_begin, begin_date, Frac.Year)
	
	# also make "nesting season" that starts April and ends March
	
	data.0 %>% mutate(begin_date = as.Date(paste(Year_begin, Month_begin, '01', sep = "-"),format = "%Y-%m-%d")) %>%
		mutate(Year = Year_begin,
					 Month = Month_begin,
					 f_month = as.factor(Month),
					 f_year = as.factor(Year),
					 Frac.Year = Year + (Month_begin-0.5)/12) %>%
		select(Year, Month, Frac.Year, begin_date, Nests) %>%
		na.omit() %>%
		right_join(.,data.2, by = "begin_date") %>%
		transmute(Year = Year.y,
							Month = Month_begin,
							Frac.Year = Frac.Year.y,
							Nests = Nests,
							Season = ifelse(Month < 4, Year-1, Year),
							Seq.Month = ifelse(Month < 4, Month + 9, Month - 3)) %>%
		dplyr::arrange(.,"Frac.Year") %>%
		filter(Season >= season.begin & Season <= season.end) -> data.1
	
	data.1 %>% filter(Month > 3 & Month < 10) -> data.summer
	data.1 %>% filter(Month > 9 | Month < 4) %>%
		mutate(Seq.Month = Seq.Month - 6) -> data.winter
	
	jags.data <- list(y = log(data.1$Nests),
										m = data.1$Seq.Month,
										T = nrow(data.1))
	
	y <- matrix(log(data.1$Nests), ncol = 12, byrow = TRUE)
	
	jags.data2 <- list(y = y,
										 m = matrix(data.1$Seq.Month, ncol = 12, byrow = TRUE),
										 n.years = nrow(y))
	
	y.summer <- matrix(log(data.summer$Nests),
										 ncol = 6, byrow = TRUE)
	
	y.winter <- matrix(log(data.winter$Nests),
										 ncol = 6, byrow = TRUE)
	
	jags.data2.summer <- list(y = y.summer,
														m = matrix(data.summer$Seq.Month,ncol = 6, byrow = TRUE),
														n.years = nrow(y.summer))
	
	jags.data2.winter <- list(y = y.winter,
														m = matrix(data.winter$Seq.Month, ncol = 6, byrow = TRUE),
														n.years = nrow(y.winter))
	
	out <- list(jags.data = jags.data,
							jags.data2 = jags.data2,
							jags.data.summer = jags.data2.summer,
							jags.data.winter = jags.data2.winter,
							data.1 = data.1,
							data.summer = data.summer,
							data.winter = data.winter)
	return(out)
}
sum.fn <- function(U){
	meanU <- mean(U)
	varU <- var(U)
	qU <- quantile(U, probs=c(0.5,0.025, 0.975))

	lamb <- exp(U)
	meanL <- mean(lamb)
	varL <- var(lamb)
	qL <- quantile(lamb, probs=c(0.5, 0.025, 0.975))

	summ <- matrix(c(meanU, qU[1], varU, qU[2:3], meanL, qL[1], varL, qL[2:3]), ncol=1)
	rownames(summ) <- c("meanU", "medianU", "varU", "L95U", "U95U", "meanl", "medianl", "varl", "L95l", "U95l")

	return(summ)
}
vbgf_bc <- function(len, prop, VBGF){
	if(missing(len)){
		if(is.null(VBGF$Lknot)){
			with(VBGF, (1/-K)*log(((Linf*prop)-Linf)/(-Linf))+tknot)
		}else{
			with(VBGF, (1/K)*log((Linf - Lknot)/(Linf*(1-prop))))
		}
		
	}else if(missing(prop)){
		if(is.null(VBGF$Lknot)){
			with(VBGF, (1/-K)*log((len-Linf)/(-Linf))+tknot)
		}else{
			with(VBGF, (1/K)*log((Linf - Lknot)/(Linf-len)))
		}
	}
}
maturity <- function(len, casewise=FALSE, p=0.99, VBGF){
	if(casewise){
		ifelse(len>=p*VBGF$Linf,1,1/(1+exp(-(len-VBGF$Lmat)/VBGF$sig_mat)))
	}else{
		1/(1+exp(-(len-VBGF$Lmat)/VBGF$sig_mat))
	}
}
case_ane <- function(draw, mode='det', div_RI){
	draw$ANEj <- (draw$Pj^draw$YatLarge)
	draw$ANEjRI <- draw$ANEj * (1/draw$RI) #ANE

	#Fill in Adult ANE
	draw$ANEa <- 1
	draw$ANEaRI <- draw$ANEa*(1/draw$RI)
	draw$ANEt <- ((draw$ANEj * as.integer(draw$Stage=="J")) + (draw$ANEa * as.integer(draw$Stage=="A")))
	draw$ANEtRI <- ((draw$ANEjRI * as.integer(draw$Stage=="J")) + (draw$ANEaRI * as.integer(draw$Stage=="A")))
	draw$ANEtold <- ((draw$ANEjRI * as.integer(draw$Stage=="J")) + (draw$ANEa * as.integer(draw$Stage=="A")))
	#Calculate the realized ANE
	if(mode=='det'){
		draw$rANE <-  draw$ANEt * (draw$PF * draw$m)
		draw$rANERI <-  draw$ANEtRI * (draw$PF * draw$m)
		draw$rANEold <-  draw$ANEtold * (draw$PF * draw$m)
	}else{
		draw$rANE <-  draw$ANEt * (draw$PF.draw * draw$m.draw)
		draw$rANERI <-  draw$ANEtRI * (draw$PF.draw * draw$m.draw)
		draw$rANEold <-  draw$ANEtold * (draw$PF.draw * draw$m.draw)
	}
	return(draw)
}
case_ane_wrapper <- function(draw, mode='det'){
	b <- rep(0, length(draw))
	for(i in 1:length(draw)){
		b[i] <- sum(ifelse(is.na(draw[[i]]),
		                   0,
		                   case_ane(draw[[i]], 
		                            mode=mode)$rANE))
	}
	return(b)
}
prop_ane <- function(Year,max_year,Age,Pj,Pa,PF,M,RI,VBGF,mat_p){
	y <- seq(Year, max_year)
	l_y <- length(y)
	if(Age>=VBGF$max_age){
		v <- rep(VBGF$Linf*VBGF$mat_p,length(y))
	}else{
		y2mat <- Year+length(seq(Age, VBGF$max_age, by=1))-1
		if(y2mat < max_year){
			v <- c(VBGF$Linf * (1-exp(-VBGF$K*(seq(Age, VBGF$max_age, by=1)-VBGF$tknot))),rep(VBGF$Linf*VBGF$mat_p,max_year-y2mat))
		}else{
			ages <- seq(Age, Age+y2mat, by=1)[1:l_y]
			v <- VBGF$Linf * (1-exp(-VBGF$K*(ages-VBGF$tknot)))
		}
	}

	p <- maturity(v, casewise=TRUE, p=VBGF$mat_p, VBGF=VBGF)
	s <- cumprod((1-p) * Pj + p * Pa)
	p <- rbinom(l_y,1,p)
	if(any(p==1)) p[min(which(p==1)):l_y] <- 1
	sr <- sr2 <- s * p * M * PF * (1/RI)
	if(any(p==1)){
		sr2[pmin(l_y,min(which(p==1))+1):l_y] <- 0
	}else{
		sr2 <- rep(0,l_y)
	}
	xx <- as.data.frame(cbind(y,v,p,s,sr,sr2,
	                    PF=rep(PF,l_y),
	                    M = rep(M, l_y),
	                    RI = rep(RI, l_y)))

	return(xx)
}
draw_propagate <- function(ATL, scenario){
	
	l <- with(scenario,{
		det.TD.draw <- rep(list(NA), n_y) #storage list
		sto.TD.draw <- rep(list(NA), n_y) #storage list
		det.TD.draw.raw <- rep(list(NA),n_y)
		sto.TD.draw.raw <- rep(list(NA),n_y)
		#draw a proportion female in the Take for each year
		PF.sto <- rnorm(n_y, RV$PF, RV$PF_sd)
		PF.det <- RV$PF
		mu.l <- TD_MVN$beta0 + TD_MVN$beta1*ATL #linear model for mean lengths in take (log)
		mu.m <- rep(TD_MVN$mu0, n_y) # mean discard mortality in take (logit)
		for(i in 1:length(ATL)){
			if(ATL[i] > 0){
				#draw from MVN log(length) and logit(discard mortality)
				det.TD.draw.raw[[i]] <- sto.TD.draw.raw[[i]] <- rmvnorm(ATL[i], c(mu.l[i], mu.m[i]), TD_MVN$cov)
				#convert to length and discard mortality (prob. space)
				det.TD.draw[[i]] <- sto.TD.draw[[i]] <- data.frame(l = exp(sto.TD.draw.raw[[i]][,1]), m = inv.logit(sto.TD.draw.raw[[i]][,2]))
				#get an anticpated Age from the VBGF
				#improved with back-calculation
				det.TD.draw[[i]]$Age <- sto.TD.draw[[i]]$Age <- sapply(sto.TD.draw[[i]]$l, function(x) {v <- vbgf_bc(len=x, VBGF=VBGF); ifelse(is.nan(v),vbgf_bc(prop=0.99, VBGF=VBGF),v)})

				#### Stochastic
				#draw a remigration interval for every animal based on a specified distribution (normal or CMP)
				det.TD.draw[[i]]$RI <- RV$RI
				if(RI.dist =="normal"){
					sto.TD.draw[[i]]$RI <- rtruncnorm(nrow(sto.TD.draw[[i]]), a=0, mean=RV$RI, sd=RV$RI_sd) #remigration interval
				}else if(RI.dist =="CMP"){
					sto.TD.draw[[i]]$RI <- rCMP(nrow(sto.TD.draw[[i]]), mu=RV$RI, nu=RV$RI_sd) #remigration interval
					while(any(sto.TD.draw[[i]]$RI==0)){
						sto.TD.draw[[i]]$RI[sto.TD.draw[[i]]$RI==0] <- rCMP(sum(sto.TD.draw[[i]]$RI==0), mu = RV$RI, nu = RV$RI_sd)
					}
				}else{
					sto.TD.draw[[i]]$RI <- RV$RI
					#error handling
					print("Distribution was not specified as normal or CMP")
				}
				#draw a juvenile survival for every animal
				det.TD.draw[[i]]$Pj <- PjV$mu
				sto.TD.draw[[i]]$Pj <- rnorm(nrow(sto.TD.draw[[i]]), PjV$mu, PjV$sd) #juvenile survival
				#draw an adult survival for every animal
				det.TD.draw[[i]]$Pa <- PaV$mu
				sto.TD.draw[[i]]$Pa <- rnorm(nrow(sto.TD.draw[[i]]), PaV$mu, PaV$sd) #juvenile survival


				#calculate the number of years until maturity
				det.TD.draw[[i]]$YatLarge <- sto.TD.draw[[i]]$YatLarge <- VBGF$Amat - sto.TD.draw[[i]]$Age
				
				#determine the animals stage based on maturity
				det.TD.draw[[i]]$Stage <- sto.TD.draw[[i]]$Stage <- ifelse(sto.TD.draw[[i]]$Age > VBGF$Amat, "A", "J")
				#Draw whether the animal is M/F based on yearly proportion
				det.TD.draw[[i]]$PF <- PF.det
				sto.TD.draw[[i]]$PF.draw <- rbinom(nrow(sto.TD.draw[[i]]), 1, 
				                                   PF.sto[i])
				

				#Draw whether the animal lived or died based on discard mortality
				if(grim.reaper){
					det.TD.draw[[i]]$m <- sto.TD.draw[[i]]$m <- 1
				}else{
					sto.TD.draw[[i]]$m.draw <- rbinom(nrow(sto.TD.draw[[i]]), 1, sto.TD.draw[[i]]$m)
				}
			}
		}
		return(list(det=det.TD.draw,
		            sto=sto.TD.draw))
	})
	return(l)
}
draw_atl <- function(scenario){
	with(scenario,{
		if(ATL_scale){
			ATL <- rep(0, n_y)
			for(i in 1:n_y){
				if(ATLp$type=="poisson"){
					ATL[i] <- rpois(1, lambda=ATLp$lambda*cumprod(lambda)[i], nu = ATLp$nu1)
				}else{
					ATL[i] <- rCMP(1, mu = ATLp$mu*cumprod(lambda)[i], nu = ATLp$nu)
				}
			}
		}else{
			ATL <- matrix(0, n_y, 3)
			if(ATLp$type=="poisson"){
				ATL <- rpois(n_y, lambda = ATLp$lambda)
			}else{
				ATL <- rCMP(n_y, mu = ATLp$mu, nu = ATLp$nu)
			}
		}
		return(ATL)
	})
}
draw_ane <- function(ATL, scenario, draw){
	with(scenario,{
		det.TD.draw <- draw$det
		sto.TD.draw <- draw$sto
		ANEdraw.det <- rep(list(NA),n_y) #storage list of MVN draws
		det.Take.old.nRI <- det.Take.old.RI <- det.Take.old <- det.Take <- det.Take.n1 <- rep(0, n_y) #storage vector of take
		ANEdraw.sto <- rep(list(NA),n_y) #storage list of MVN draws
		sto.Take.old.nRI <- sto.Take.old.RI <- sto.Take.old <- sto.Take <- sto.Take.n1 <- rep(0, n_y) #storage vector of take
		for(i in 1:n_y){
			if(ATL[i]>0){
					
				#### Stochastic
				sto.TD.draw[[i]] <- case_ane(sto.TD.draw[[i]],
				                             mode='sto',
				                             div_RI=div_RI)
				
				#### Deterministic
				det.TD.draw[[i]] <- case_ane(det.TD.draw[[i]],
				                             div_RI=div_RI)

				det.Take.old[i] <- sum(det.TD.draw[[i]]$rANEold)
				sto.Take.old[i] <- sum(sto.TD.draw[[i]]$rANEold)
				det.Take.old.nRI[i] <- sum(det.TD.draw[[i]]$rANE)
				sto.Take.old.nRI[i] <- sum(sto.TD.draw[[i]]$rANE)
				det.Take.old.RI[i] <- sum(det.TD.draw[[i]]$rANERI)
				sto.Take.old.RI[i] <- sum(sto.TD.draw[[i]]$rANERI)		
				#Calculate the Adult Nester Equivalency
				
				ANEdraw.det[[i]] <- sapply(1:ATL[i],function(j){
					prop_ane(Year = i, 
                     max_year = n_y,
                     Age = det.TD.draw[[i]]$Age[j],
                     Pj = det.TD.draw[[i]]$Pj[j],
                     Pa = det.TD.draw[[i]]$Pa[j],
                     PF = det.TD.draw[[i]]$PF[j],
                     RI = det.TD.draw[[i]]$RI[j],
                     M = det.TD.draw[[i]]$m[j],
                     VBGF=VBGF)},
                     simplify=FALSE)
				ANEdraw.sto[[i]] <- sapply(1:ATL[i],function(j){
					prop_ane(Year = i, 
                     max_year = n_y,
                     Age = sto.TD.draw[[i]]$Age[j],
                     Pj = sto.TD.draw[[i]]$Pj[j],
                     Pa = sto.TD.draw[[i]]$Pa[j],
                     PF = sto.TD.draw[[i]]$PF.draw[j],
                     RI = sto.TD.draw[[i]]$RI[j],
                     M = sto.TD.draw[[i]]$m.draw[j],
                     VBGF=VBGF)},
                     simplify=FALSE)
				ANEdraw.det.agg <- aggregate(cbind(sr,sr2)~y,
				                             data=do.call('rbind',ANEdraw.det[[i]]),
				                             FUN=sum)
				ANEdraw.sto.agg <- aggregate(cbind(sr,sr2)~y,
				                             data=do.call('rbind',ANEdraw.sto[[i]]),
				                             FUN=sum)
				det.Take[i:n_y] <- det.Take[i:n_y] + ANEdraw.det.agg[,2]
				sto.Take[i:n_y] <- sto.Take[i:n_y] + ANEdraw.sto.agg[,2]
				det.Take.n1[i:n_y] <- det.Take.n1[i:n_y] + ANEdraw.det.agg[,3]
				sto.Take.n1[i:n_y] <- sto.Take.n1[i:n_y] + ANEdraw.sto.agg[,3]
			}
		}
		return(list(old = data.frame(det=det.Take.old, 
										            sto=sto.Take.old,
										            det.nRI=det.Take.old.nRI, 
										            sto.nRI=sto.Take.old.nRI,
										            det.RI=det.Take.old.RI, 
										            sto.RI=sto.Take.old.RI),
		            new = data.frame(det=det.Take, 
										            sto=sto.Take,
										            det.n1=det.Take.n1, 
										            sto.n1=sto.Take.n1)))
	})
}
ane_sim <- function(n, scenario, mc.cores=4){
	fn <- function(scenario){
		ATL <- draw_atl(scenario)
	
		draw <- draw_propagate(ATL, scenario)

		ane_draw <- draw_ane(ATL, scenario, draw)

		res <- cbind(ane_draw$old, ane_draw$new[,3:4])

		scenario$VBGF$mat_p <- 0.975
		ane_draw <- draw_ane(ATL, scenario, draw)$new[,3:4]
		colnames(ane_draw) <- paste0(colnames(ane_draw),'975')
		res <- cbind(res,ane_draw)
		return(as.matrix(res))
	}
	library(parallel)
	reps <- simplify2array(mclapply(1:n,function(i,...)fn(scenario),mc.cores=pmin(mc.cores,detectCores()-2)))
	return(reps)
}
proj_pop <- function(sim, trend, take, scenario, lambda, Q.sd){
	res <- data.frame(Nt_take = rep(0, scenario$n_y),
	                  Nt = rep(0, scenario$n_y),
	                  Nt_take.iQ = rep(0, scenario$n_y),
	                  Nt.iQ = rep(0, scenario$n_y))

	if(scenario$sim_y  > 0){
		res$Nt_take.iQ[1] <- rnorm(1,(trend$N0[sim] - scenario$hist_ANE[1])*lambda[1], Q.sd[1])
		res$Nt.iQ[1] <- rnorm(1, trend$N0[sim] * lambda[1], Q.sd[1])
		res$Nt_take[1] <- (trend$N0[sim] - scenario$hist_ANE[1])*lambda[1]
		res$Nt[1] <-  trend$N0[sim] * lambda[1]
		for(i in 2:(scenario$sim_y-1)){
			#include process variance Q
			res$Nt_take.iQ[i] <- rnorm(1,(res$Nt_take.iQ[i-1] - scenario$hist_ANE[i])*lambda[i], Q.sd[i])
			res$Nt.iQ[i] <- rnorm(1, (res$Nt.iQ[i-1]- scenario$hist_ANE[i]) * lambda[i], Q.sd[i])
			#drop process variance Q
			res$Nt_take[i] <- (res$Nt_take[i-1] - scenario$hist_ANE[i])*lambda[i]
			res$Nt[i] <-  (res$Nt[i-1]- scenario$hist_ANE[i]) * lambda[i]
		}
		for(i in scenario$sim_y:scenario$n_y){
			#include process variance Q
			res$Nt_take.iQ[i] <- rnorm(1,(res$Nt_take.iQ[i-1] - take[(i - (scenario$sim_y-1))])*lambda[i], Q.sd[i])
			res$Nt.iQ[i] <- rnorm(1, res$Nt.iQ[i-1] * lambda[i], Q.sd[i])
			#drop process variance Q
			res$Nt_take[i] <- (res$Nt_take[i-1] - take[(i - (scenario$sim_y-1))])*lambda[i]
			res$Nt[i] <-  res$Nt[i-1] * lambda[i]
		}
	}else{
		res$Nt_take.iQ[1] <- rnorm(1,(trend$N0[sim] - take[1])*lambda[1], Q.sd[1])
		res$Nt.iQ[1] <- rnorm(1, trend$N0[sim] * lambda[1], Q.sd[1])
		res$Nt_take[1] <- (trend$N0[sim] - take[1])*lambda[1]
		res$Nt[1] <-  trend$N0[sim] * lambda[1]
		
		for(i in 2:scenario$n_y){
			#include process variance Q
			res$Nt_take.iQ[i] <- rnorm(1,(res$Nt_take.iQ[i-1] - take[i])*lambda[i], Q.sd[i])
			res$Nt.iQ[i] <- rnorm(1, res$Nt.iQ[i-1] * lambda[i], Q.sd[i])
			#drop process variance Q
			res$Nt_take[i] <- (res$Nt_take[i-1] - take[i])*lambda[i]
			res$Nt[i] <-  res$Nt[i-1] * lambda[i]
		}
	}
	
	return(res)
}
pva.proj <- function(sim, trend, scenario, ane_ver='new', ane_calc_det='det', ane_calc_sto='sto'){
	
	#convert U to lambda
		lambda <- rep(1, scenario$n_y) #exp of U
		Q.sd <- rep(NA, scenario$n_y) #std. dev. of Q (is variance)
		#toggle for dynamic or static U & Q
		if(scenario$dynUQ){
			samp.sim <- sample(1:length(trend$U), scenario$n_y-1, replace = FALSE)
			lambda[1] <- exp(trend$U[sim])
			lambda[2:scenario$n_y] <- exp(trend$U[samp.sim])
			Q.sd[1] <- sqrt(trend$Q[sim])
			Q.sd[2:scenario$n_y] <- sqrt(trend$Q[samp.sim])
		}else{
			lambda[1:scenario$n_y] <- exp(trend$U[sim])
			Q.sd[1:scenario$n_y] <- sqrt(trend$Q[sim])
		}
		
	###-----------------------------------------------------
	#draw a ATL
		#ATL_scale is a logical flag for scaling the ATL at the same rate as the population growth rate

		ATL <- draw_atl(scenario)
	#draw a Size and Mortality based on ATL
		
		draw <- suppressWarnings(draw_propagate(ATL, scenario))

		ane_draw <- draw_ane(ATL, scenario, draw)

			# scenario$VBGF$mat_p <- 0.975
			# ane_draw2 <- draw_ane(ATL, scenario, draw)$new[,3:4]
			# #compare ANE plot
			# 	par(mar=c(2,3,1,1),mfrow=c(1,1), oma=c(2,2,0,0), las=1)
			# 	ylim.r <- c(0,max(sapply(ane_draw,max)))
			# 	# matplot(cbind(ane_draw$old[,1:2], ane_draw$new[,1:2]), 
			# 	#         type='l', lty=c(1,3,1,3), col=c(1,1,2,2), lwd=2,
			# 	#         ylim=ylim.r)
			# 	# legend('topleft', c('Casewise-det',
			# 	#                     'Casewise-sto',
			# 	#                     'p(Mature)-det',
			# 	#                     'p(Mature)-sto'),
			# 	#        lty=c(1,3,1,3), lwd=3, col=c(1,1,2,2), ncol=2)
			# 	# legend('topright',legend=round(colSums(cbind(ane_draw$old[,1:2], ane_draw$new[,1:2])),1),
			# 	#        lty=c(1,3,1,3), lwd=3, col=c(1,1,2,2))
			# 	matplot(cbind(ane_draw$old[,1:2], ane_draw$new[,3:4], ane_draw2), 
			# 	        type='l', lty=c(1,3,1,3,1,3), col=c(1,1,2,2,4,4), lwd=2,
			# 	        ylim=ylim.r)
			# 	legend('topleft', c('Casewise-det',
			# 	                    'Casewise-sto',
			# 	                    'p(Mature)-det-n1',
			# 	                    'p(Mature)-sto-n1',
			# 	                    'p(Mature)-det-n1-97.5%',
			# 	                    'p(Mature)-sto-n1-97.5%'),
			# 	       lty=c(1,3,1,3,1,3), lwd=3, col=c(1,1,2,2,4,4), ncol=1)
			# 	legend('topright',legend=round(colSums(cbind(ane_draw$old[,1:2], ane_draw$new[,3:4], ane_draw2)),1),
			# 	       lty=c(1,3,1,3,1,3), lwd=3, col=c(1,1,2,2,4,4))
			# 	matplot(ane_draw$old[,c(1,3,5)], type='l', lty=1, col=1:3, lwd=2,
			# 	        ylim=ylim.r)
			# 	legend('topleft', c('Martin et al. 2020','no-RI corr.','full RI corr.'),
			# 	       lty=1, lwd=3, col=1:3)
			# 	legend('topright',legend=round(colSums(ane_draw$old[,c(1,3,5)]),1),
			# 	       lty=1, lwd=3, col=1:3)
			# 	matplot(cbind(ane_draw$old[,5:6], ane_draw$new[,3:4]), 
			# 	        type='l', lty=c(1,3,1,3), col=c(1,1,2,2), lwd=2,
			# 	        ylim=ylim.r)
			# 	legend('topleft', c('Casewise-det-fRI',
			# 	                    'Casewise-sto-fRI',
			# 	                    'p(Mature)-det-n1',
			# 	                    'p(Mature)-sto-n1'),
			# 	       lty=c(1,3,1,3), lwd=3, col=c(1,1,2,2), ncol=2)
			# 	legend('topright',legend=round(colSums(cbind(ane_draw$old[,5:6], ane_draw$new[,3:4])),1),
			# 	       lty=c(1,3,1,3), lwd=3, col=c(1,1,2,2))
			# 	mtext('Years',side=1,outer=TRUE)
			# 	mtext(expression(ANE[y]), side=2, outer=TRUE, las=0)
		

		det.proj <- proj_pop(sim=sim, 
		                 trend=trend, 
		                 take=ane_draw[[ane_ver]][,ane_calc_det], 
		                 scenario = scenario,
		                 lambda = lambda,
		                 Q.sd = Q.sd)
		sto.proj <- proj_pop(sim=sim, 
		                 trend=trend, 
		                 take=ane_draw[[ane_ver]][,ane_calc_sto], 
		                 scenario = scenario,
		                 lambda = lambda,
		                 Q.sd = Q.sd)
		colnames(det.proj) <- paste0('det.',colnames(det.proj))
		colnames(sto.proj) <- paste0('sto.',colnames(sto.proj))

		proj <- cbind(det.proj, sto.proj)

	return(as.matrix(proj))
}
curr.abund.fn <- function(dat, RI, round=TRUE){
	curr.abund <- matrix(c(quantile(dat$N_fym0, probs=c(0.5,0.025,0.975)), quantile(dat$N_fym1, probs=c(0.5,0.025,0.975)), quantile(dat$N_fym2, probs=c(0.5,0.025,0.975)), quantile(dat$N_fym3, probs=c(0.5,0.025,0.975)), quantile(rowSums(dat[,c("N_fym0","N_fym1","N_fym2","N_fym3")]),probs=c(0.5,0.025,0.975))*(RI/4)), nrow=5, byrow=T)
	rownames(curr.abund) <- c("N0","N-1","N-2","N-3","Sum")
	colnames(curr.abund) <- c("Median","L95%","U95%")
	if(round) curr.abund <- round(curr.abund)
	return(curr.abund)
}
sim.fn <- function(trend, max.cores, scenario, nsim, fn, ane_ver='new', ane_calc_det='det', ane_calc_sto='sto'){
	pva.proj <- fn
	UseCores <- pmin(max.cores,detectCores() - 2)
	cl <- makeCluster(UseCores)
	registerDoParallel(cl)
	nst <- ceiling(seq(from=1, to = nsim, length.out=(UseCores+1)))[-(UseCores+1)]
	nen <- ceiling(seq(from=1, to = nsim, length.out=(UseCores+1)))[-1] - c(rep(1,(UseCores-1)),0)
	trial <- pva.proj(1, trend, scenario, ane_ver='new', ane_calc_det='det', ane_calc_sto='sto')
	df <- foreach(i = 1:length(nst), .packages=c("mvtnorm","truncnorm"), .export=c('draw_atl','draw_ane','draw_propagate','case_ane','prop_ane','rCMP','dCMP','vbgf_bc','maturity','proj_pop','compute_CMP_constant','inv.logit','logit')) %dopar% {

		iseq <- seq(nst[i],nen[i])
		arr <- array(NA, dim=c(scenario$n_y, ncol(trial), length(iseq)))

		for(j in 1:length(iseq)){
			arr[,,j] <- pva.proj(iseq[j], trend, scenario, ane_ver, ane_calc_det, ane_calc_sto)
		}

		return(arr)
	}
	stopCluster(cl)
	return(df)
}
bootCI <- function(dt, R){
	samps <- replicate(R, sample(1:length(dt), floor(0.9*length(dt))))
	q <- apply(samps, 2, function(v) {sum(dt[v]/length(v))})
	return(q)
}
fig.lab <- function(i, xscale=0.05, yscale, cex=1.4, adj=c(0.5,0.5)){
	text(x = par('usr')[1] + abs(diff(par('usr')[1:2]))*xscale,
	     y = par('usr')[3] + abs(diff(par('usr')[3:4]))*yscale,
	     ifelse(is.numeric(i),LETTERS[i],i), xpd=NA, cex=cex, adj=adj)
}
u.den <- function(notake, notake.old, take){
	d1 <- density(notake, adj=2)
	d2 <- density(take, adj=2)
	d3 <- density(notake.old, adj=2)
	d1$y <- d1$y/max(d1$y)
	d2$y <- d2$y/max(d2$y)
	d3$y <- d3$y/max(d3$y)

	plot(d1$x, d1$y, xlab="r", ylab="Relative Density", ylim=c(0,1.01), las=1, col='chartreuse4', type='l', lwd=3)
	lines(d2$x, d2$y, col='dodgerblue4', lwd=3, lty=2)
	lines(d3$x, d3$y, col='goldenrod4', lwd=3, lty=2)
	abline(v=median(notake), col='chartreuse4')
	abline(v=median(take), col='dodgerblue4')
	abline(v=median(notake.old), col='goldenrod4')
	legend("topleft", legend=c(expression(N),
	                           expression(N[old]), 
	                           expression(N+F)), 
	       lwd=3, col=c('chartreuse4','goldenrod4','dodgerblue4'), bty='n')
}
col2rgbA<-function(color,transparency){
  rgb(t(col2rgb(color))/255,alpha=transparency)
}
proj.summ.fn <- function(sim.l, trend, thres.pop, thres.yr, alpha=0.05, keepers=c(3,4,7,8), spp, mode=NULL, first.yr=5, yr.low=2018, yr.high=2124){

	lci <- alpha/2
	uci <- 1-(alpha/2)

	arr <- abind(sim.l, along=3)
	check.seq <- first.yr:dim(arr)[1]
	check.l <- length(check.seq)

	### Summary arrays
	#dims = length(thres.pop), dim(arr)[2] modes from pva.proj, 4 summary metrics (mean, median, LCI, UCI)
	tab.yr <- array(NA, dim=c(length(thres.pop),dim(arr)[2], 4))
	#dims = length(thres.yr), length(thres.pop), dim(arr)[2] modes from pva.proj, 3 summary metrics (median, LCI, UCI)
	tab.prob <- array(NA, dim=c(length(thres.yr), length(thres.pop), dim(arr)[2], 3))

	### Storage arrays
	#dims = length(thres.pop), dim(arr)[2] modes from pva.proj, nsims
	yr <- array(NA, dim=c(length(thres.pop), dim(arr)[2], dim(arr)[3]))
	#dims = length(thres.yr), length(thres.pop), dim(arr)[2] modes from pva.proj, nsims
	prob <- array(NA, dim=c(length(thres.yr), length(thres.pop), dim(arr)[2], dim(arr)[3]))

	#duplicate array to fill in extinct years with 0's
	arr2 <- arr

	for(i in 1:dim(arr)[3]){
		N0 <- trend$N0[i]
		N0.thres <- N0*thres.pop

		#fills in zeros for extinct runs
		arr2[,,i] <- apply(arr2[,,i], 2, function(x) {if(any(x < 0)){ d <- min(which(x < 0)); x[d:dim(arr)[1]] <- 0}; return(x)})

		#first year to fall below threshold
		yr[,,i] <- apply(arr2[check.seq,,i], 2, FUN = function(x) {sapply(N0.thres, function(v) ifelse(min(which((x < v)==1))==Inf, NA, min(which((x < v)==1))))})

		# prob[,,,i] <- array(apply(arr2[,,i], 2, function(x) {sapply(N0.thres, function(v){ sapply(thres.yr, function(y) sum(x[1:y] < v)/y)})}), dim=c(5,3,dim(arr)[2]))
		prob[,,,i] <- array(apply(arr2[check.seq,,i], 2, function(x) {sapply(N0.thres, function(v){ sapply(thres.yr, function(y) as.integer(any(x[1:y] < v)))})}), dim=c(5,3,dim(arr)[2]))
	}

	#Summarize Years to threshold
	tab.yr[,,1] <- apply(yr, c(1,2), mean, na.rm=T)
	tab.yr[,,2] <- apply(yr, c(1,2), median, na.rm=T)
	tab.yr[,,3] <- apply(yr, c(1,2), quantile, prob = lci, na.rm=T)
	tab.yr[,,4] <- apply(yr, c(1,2), quantile, prob = uci, na.rm=T)

	#Summarize probability of falling below threshold
	tab.prob[,,,1] <- apply(prob, c(1,2,3), function(x) {sum(x)/length(x)})
	bootci.prob <- apply(prob,c(1,2,3), bootCI, R=1000)
	tab.prob[,,,2] <- apply(bootci.prob, c(2,3,4), quantile, prob=lci)
	tab.prob[,,,3] <- apply(bootci.prob, c(2,3,4), quantile, prob=uci)

	# Table 1
		tab1.det.NT <- matrix(c(1-aperm(tab.prob[,,keepers[2],], c(2,1,3))[,5,1], aperm(tab.prob[,,keepers[2],], c(2,1,3))[,5,1], tab.yr[,keepers[2],]), nrow=3)
		tab1.det.T <- matrix(c(1-aperm(tab.prob[,,keepers[1],], c(2,1,3))[,5,1], aperm(tab.prob[,,keepers[1],], c(2,1,3))[,5,1], tab.yr[,keepers[1],]), nrow=3)
		tab1.sto.NT <- matrix(c(1-aperm(tab.prob[,,keepers[4],], c(2,1,3))[,5,1], aperm(tab.prob[,,keepers[4],], c(2,1,3))[,5,1], tab.yr[,keepers[4],]), nrow=3)
		tab1.sto.T <- matrix(c(1-aperm(tab.prob[,,keepers[3],], c(2,1,3))[,5,1], aperm(tab.prob[,,keepers[3],], c(2,1,3))[,5,1], tab.yr[,keepers[3],]), nrow=3)

		rownames(tab1.det.NT) <- rownames(tab1.det.T) <- rownames(tab1.sto.NT) <- rownames(tab1.sto.T) <- paste0(thres.pop*100, "%")
		colnames(tab1.det.NT) <- colnames(tab1.det.T) <- colnames(tab1.sto.NT) <- colnames(tab1.sto.T) <- c("Prob.Above","Prob.Below","MeanYr","MedYr","L95Yr","U95Yr")

		write.csv(round(tab1.det.NT,2), file=paste0(table.path,"Table1.",spp, mode,".det.NT.csv"))
		write.csv(round(tab1.det.T,2), file=paste0(table.path,"Table1.",spp, mode,".det.T.csv"))
		write.csv(round(tab1.sto.NT,2), file=paste0(table.path,"Table1.",spp, mode,".sto.NT.csv"))
		write.csv(round(tab1.sto.T,2), file=paste0(table.path,"Table1.",spp, mode,".sto.T.csv"))

		tab.sum.det <- rbind(tab1.det.NT[1,], tab1.det.T[1,], 
		                     tab1.det.NT[1,] - tab1.det.T[1,],
		                     tab1.det.NT[2,], tab1.det.T[2,], 
		                     tab1.det.NT[2,] - tab1.det.T[2,],
		                     tab1.det.NT[3,], tab1.det.T[3,], 
		                     tab1.det.NT[3,] - tab1.det.T[3,])
		tab.sum.sto <- rbind(tab1.sto.NT[1,], tab1.sto.T[1,], 
		                     tab1.sto.NT[1,] - tab1.sto.T[1,],
		                     tab1.sto.NT[2,], tab1.sto.T[2,], 
		                     tab1.sto.NT[2,] - tab1.sto.T[2,],
		                     tab1.sto.NT[3,], tab1.sto.T[3,], 
		                     tab1.sto.NT[3,] - tab1.sto.T[3,])
		rownames(tab.sum.det) <- rownames(tab.sum.sto) <- paste(rep(paste0(thres.pop*100, "%"), each=3), c("No Take","Take", "∆(NT-T)"))

		write.csv(round(tab.sum.det,2), file=paste0(table.path,"Table1.",spp, mode,".det.SUM.csv"))
		write.csv(round(tab.sum.sto,2), file=paste0(table.path,"Table1.",spp, mode,".sto.SUM.csv"))
	# Table 2
		tp.det.NT <- apply(aperm(tab.prob[,,keepers[2],], c(2,1,3)),2,c)[c(1,4,7,2,5,8,3,6,9),]
		tp.det.T <- apply(aperm(tab.prob[,,keepers[1],], c(2,1,3)),2,c)[c(1,4,7,2,5,8,3,6,9),]
		tp.sto.NT <- apply(aperm(tab.prob[,,keepers[4],], c(2,1,3)),2,c)[c(1,4,7,2,5,8,3,6,9),]
		tp.sto.T <- apply(aperm(tab.prob[,,keepers[3],], c(2,1,3)),2,c)[c(1,4,7,2,5,8,3,6,9),]

		rownames(tp.det.NT) <- rownames(tp.det.T) <- rownames(tp.sto.NT) <- rownames(tp.sto.T) <- paste0(rep(paste0(thres.pop*100,"%"),each=3),  c("", "L95","U95"))
		colnames(tp.det.NT) <- colnames(tp.det.T) <- colnames(tp.sto.NT) <- colnames(tp.sto.T) <- paste0(thres.yr,'yr')

		write.csv(round(tp.det.NT,3), file=paste0(table.path,"Table2.",spp, mode,".det.NT.csv"))
		write.csv(round(tp.det.T,3), file=paste0(table.path,"Table2.",spp, mode,".det.T.csv"))
		write.csv(round(tp.sto.NT,3), file=paste0(table.path,"Table2.",spp, mode,".sto.NT.csv"))
		write.csv(round(tp.sto.T,3), file=paste0(table.path,"Table2.",spp, mode,".sto.T.csv"))

		tab.perc.det <- rbind(tp.det.NT[1,], tp.det.T[1,], 
		                     tp.det.NT[1,] - tp.det.T[1,],
		                     tp.det.NT[2,], tp.det.T[2,], 
		                     tp.det.NT[2,] - tp.det.T[2,],
		                     tp.det.NT[3,], tp.det.T[3,], 
		                     tp.det.NT[3,] - tp.det.T[3,],
		                     tp.det.NT[4,], tp.det.T[4,], 
		                     tp.det.NT[4,] - tp.det.T[4,],
		                     tp.det.NT[5,], tp.det.T[5,], 
		                     tp.det.NT[5,] - tp.det.T[5,],
		                     tp.det.NT[6,], tp.det.T[6,], 
		                     tp.det.NT[6,] - tp.det.T[6,],
		                     tp.det.NT[7,], tp.det.T[7,], 
		                     tp.det.NT[7,] - tp.det.T[7,],
		                     tp.det.NT[8,], tp.det.T[8,], 
		                     tp.det.NT[8,] - tp.det.T[8,],
		                     tp.det.NT[9,], tp.det.T[9,], 
		                     tp.det.NT[9,] - tp.det.T[9,])
		tab.perc.sto <- rbind(tp.sto.NT[1,], tp.sto.T[1,], 
		                     tp.sto.NT[1,] - tp.sto.T[1,],
		                     tp.sto.NT[2,], tp.sto.T[2,], 
		                     tp.sto.NT[2,] - tp.sto.T[2,],
		                     tp.sto.NT[3,], tp.sto.T[3,], 
		                     tp.sto.NT[3,] - tp.sto.T[3,],
		                     tp.sto.NT[4,], tp.sto.T[4,], 
		                     tp.sto.NT[4,] - tp.sto.T[4,],
		                     tp.sto.NT[5,], tp.sto.T[5,], 
		                     tp.sto.NT[5,] - tp.sto.T[5,],
		                     tp.sto.NT[6,], tp.sto.T[6,], 
		                     tp.sto.NT[6,] - tp.sto.T[6,],
		                     tp.sto.NT[7,], tp.sto.T[7,], 
		                     tp.sto.NT[7,] - tp.sto.T[7,],
		                     tp.sto.NT[8,], tp.sto.T[8,], 
		                     tp.sto.NT[8,] - tp.sto.T[8,],
		                     tp.sto.NT[9,], tp.sto.T[9,], 
		                     tp.sto.NT[9,] - tp.sto.T[9,])

		rownames(tab.perc.det) <- rownames(tab.perc.sto) <- paste(rep(paste0(thres.pop*100, "%"), each=9), rep(c("","L95", "U95"),each=3), c("No Take","Take", "∆(NT-T)"))

		write.csv(round(tab.perc.det,3), file=paste0(table.path,"Table2.",spp, mode,".det.SUM.csv"))
		write.csv(round(tab.perc.sto,3), file=paste0(table.path,"Table2.",spp, mode,".sto.SUM.csv"))

	#Median Projections
		proj.med <- apply(arr2, c(1,2), quantile, probs=c(lci, 0.5, uci))
		proj.log.med <- apply(arr2, c(1,2), function(x) {quantile(log(x), probs=c(lci, 0.5, uci))})
		proj.log.med[is.infinite(proj.log.med)] <- min(proj.log.med[!is.infinite(proj.log.med)]) * 5
		diff.det <- arr2[,keepers[1],] - arr2[,keepers[2],]
		diff.sto <- arr2[,keepers[3],] - arr2[,keepers[4],]
		proj.diff.det <- apply(diff.det, 1, quantile, probs=c(lci, 0.5, uci))
		proj.diff.sto <- apply(diff.sto, 1, quantile, probs=c(lci, 0.5, uci))

			png(file=paste0(fig.path,spp,mode,"proj100.png"), units = 'in', res=300, width=7, height=7)
				###PROJECT 100
				layout(matrix(c(1,2,3,3),2,2), width=c(1,0.3))
				par(mar=c(3,4,2,1))
				# Deterministic
				plot(yr.low:yr.high, proj.med[2,,keepers[1]], 
				     xlim=c((yr.low-1),(yr.high+1)), 
				     ylim=c(0,max(proj.med[,,keepers])),
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab="Annual Nesters", 
				     yaxs='i', xaxs='i')
				polygon(x = c(yr.low:yr.high, rev(yr.low:yr.high)), y = c(proj.med[1,,keepers[1]], rev(proj.med[3,,keepers[1]])), border=FALSE, col=col2rgbA("dodgerblue3", 0.5))
				polygon(x = c(yr.low:yr.high, rev(yr.low:yr.high)), y = c(proj.med[1,,keepers[2]], rev(proj.med[3,,keepers[2]])), border=FALSE, col=col2rgbA("chartreuse3", 0.5))
				lines(yr.low:yr.high, proj.med[2,,keepers[1]], col="dodgerblue4", lwd=3)
				lines(yr.low:yr.high, proj.med[2,,keepers[2]], col="chartreuse4", lwd=3)
				yrpretty <- pretty(yr.low:(yr.high+1))[pretty(yr.low:(yr.high+1))>yr.low & pretty(yr.low:(yr.high+1)) < (yr.high+1)]
				axis(1, yrpretty, yrpretty)
				mtext("Deterministic", side=3, font=3)
				fig.lab("W. Pac. Leatherbacks", xscale=0.85, yscale=0.925, cex = 1)
				abline(v=seq(yr.low,yr.high)[first.yr])
				# Stochastic
				par(mar=c(4,4,1,1))
				plot(yr.low:yr.high, proj.med[2,,keepers[3]], 
				     xlim=c((yr.low-1),(yr.high+1)), 
				     ylim=c(0,max(proj.med[,,keepers])), 
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab="Annual Nesters", 
				     yaxs='i', xaxs='i')
				polygon(x = c(yr.low:yr.high, rev(yr.low:yr.high)), y = c(proj.med[1,,keepers[3]], rev(proj.med[3,,keepers[3]])), border=FALSE, col=col2rgbA("dodgerblue3", 0.5))
				polygon(x = c(yr.low:yr.high, rev(yr.low:yr.high)), y = c(proj.med[1,,keepers[4]], rev(proj.med[3,,keepers[4]])), border=FALSE, col=col2rgbA("chartreuse3", 0.5))
				lines(yr.low:yr.high, proj.med[2,,keepers[3]], col="dodgerblue4", lwd=3)
				lines(yr.low:yr.high, proj.med[2,,keepers[4]], col="chartreuse4", lwd=3)
				axis(1, yrpretty, yrpretty)
				mtext("Stochastic", side=3, font=3)
				fig.lab("W. Pac. Leatherbacks", xscale=0.85, yscale=0.925, cex = 1)
				abline(v=seq(yr.low,yr.high)[first.yr])
				par(mar=c(0,0,0,0))
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				legend("center", legend=c(expression(paste(N[t] -F," Median")), expression(paste(N[j] -F," 95%CI")),expression(paste(N[j]," Median")), expression(paste(N[j]," 95%CI"))), lty=c(1,NA,1,NA), lwd=c(3,NA,3,NA), pch=c(NA,15,NA,15), col=c("dodgerblue4",col2rgbA("dodgerblue3",0.5),"chartreuse4", col2rgbA("chartreuse3", 0.5)), bty='n', pt.cex = 2, xpd=NA)

			dev.off()

			png(file=paste0(fig.path,spp,mode,"logproj100.png"), units = 'in', res=300, width=7, height=7)
				###PROJECT 100
				layout(matrix(c(1,2,3,3),2,2), width=c(1,0.3))
				par(mar=c(3,4,2,1))
				# Deterministic
				plot(yr.low:yr.high, proj.log.med[2,,keepers[1]], 
				     xlim=c((yr.low-1),(yr.high+1)), 
				     ylim=c(0.0001,max(proj.log.med[,,keepers])),
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab="Log Annual Nesters", 
				     yaxs='i', xaxs='i')
				polygon(x = c(yr.low:yr.high, rev(yr.low:yr.high)), y = c(proj.log.med[1,,keepers[1]], rev(proj.log.med[3,,keepers[1]])), border=FALSE, col=col2rgbA("dodgerblue3", 0.5))
				polygon(x = c(yr.low:yr.high, rev(yr.low:yr.high)), y = c(proj.log.med[1,,keepers[2]], rev(proj.log.med[3,,keepers[2]])), border=FALSE, col=col2rgbA("chartreuse3", 0.5))
				lines(yr.low:yr.high, proj.log.med[2,,keepers[1]], col="dodgerblue4", lwd=3)
				lines(yr.low:yr.high, proj.log.med[2,,keepers[2]], col="chartreuse4", lwd=3)
				yrpretty <- pretty(yr.low:(yr.high+1))[pretty(yr.low:(yr.high+1))>yr.low & pretty(yr.low:(yr.high+1)) < (yr.high+1)]
				axis(1, yrpretty, yrpretty)
				mtext("Deterministic", side=3, font=3)
				fig.lab("W. Pac. Leatherbacks", xscale=0.85, yscale=0.925, cex = 1)
				abline(v=seq(yr.low,yr.high)[first.yr])
				# Stochastic
				par(mar=c(4,4,1,1))
				plot(yr.low:yr.high, proj.log.med[2,,keepers[3]], 
				     xlim=c((yr.low-1),(yr.high+1)), 
				     ylim=c(0.0001,max(proj.log.med[,,keepers])), 
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab="Log Annual Nesters", 
				     yaxs='i', xaxs='i')
				polygon(x = c(yr.low:yr.high, rev(yr.low:yr.high)), y = c(proj.log.med[1,,keepers[3]], rev(proj.log.med[3,,keepers[3]])), border=FALSE, col=col2rgbA("dodgerblue3", 0.5))
				polygon(x = c(yr.low:yr.high, rev(yr.low:yr.high)), y = c(proj.log.med[1,,keepers[4]], rev(proj.log.med[3,,keepers[4]])), border=FALSE, col=col2rgbA("chartreuse3", 0.5))
				lines(yr.low:yr.high, proj.log.med[2,,keepers[3]], col="dodgerblue4", lwd=3)
				lines(yr.low:yr.high, proj.log.med[2,,keepers[4]], col="chartreuse4", lwd=3)
				axis(1, yrpretty, yrpretty)
				mtext("Stochastic", side=3, font=3)
				fig.lab("W. Pac. Leatherbacks", xscale=0.85, yscale=0.925, cex = 1)
				abline(v=seq(yr.low,yr.high)[first.yr])
				par(mar=c(0,0,0,0))
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				legend("center", legend=c(expression(paste(N[j] -F," Median")), expression(paste(N[j] -F," 95%CI")),expression(paste(N[j]," Median")), expression(paste(N[j]," 95%CI"))), lty=c(1,NA,1,NA), lwd=c(3,NA,3,NA), pch=c(NA,15,NA,15), col=c("dodgerblue4",col2rgbA("dodgerblue3",0.5),"chartreuse4", col2rgbA("chartreuse3", 0.5)), bty='n', pt.cex = 2, xpd=NA)

			dev.off()

			png(file=paste0(fig.path,spp,mode,"proj10.png"), units = 'in', res=300, width=7, height=7)
				#### PROJECT 10
				layout(matrix(c(1,2,3,3),2,2), width=c(1,0.3))
				par(mar=c(3,4,2,1))
				# Deterministic
				plot(yr.low:(yr.low+9+first.yr), proj.med[2,1:(10+first.yr),keepers[1]], 
				     xlim=c((yr.low-0.25),(yr.low+10-0.75+first.yr)), 
				     ylim=c(0,max(proj.med[,1:(10+first.yr),keepers])),
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab="Annual Nesters", 
				     yaxs='i', xaxs='i')
				polygon(x = c(yr.low:(yr.low+9+first.yr), rev(yr.low:(yr.low+9+first.yr))), y = c(proj.med[1,1:(10+first.yr),keepers[1]], rev(proj.med[3,1:(10+first.yr),keepers[1]])), border=FALSE, col=col2rgbA("dodgerblue3", 0.5))
				polygon(x = c(yr.low:(yr.low+9+first.yr), rev(yr.low:(yr.low+9+first.yr))), y = c(proj.med[1,1:(10+first.yr),keepers[2]], rev(proj.med[3,1:(10+first.yr),keepers[2]])), border=FALSE, col=col2rgbA("chartreuse3", 0.5))
				lines(yr.low:(yr.low+9+first.yr), proj.med[2,1:(10+first.yr),keepers[1]], col="dodgerblue4", lwd=3)
				lines(yr.low:(yr.low+9+first.yr), proj.med[2,1:(10+first.yr),keepers[2]], col="chartreuse4", lwd=3)
				yrpretty <- pretty(yr.low:(yr.low+10+first.yr))[pretty(yr.low:(yr.low+10+first.yr))>yr.low & pretty(yr.low:(yr.low+10+first.yr)) < (yr.low+10+first.yr)]
				axis(1, yrpretty, yrpretty)
				mtext("Deterministic", side=3, font=3)
				fig.lab("W. Pac. Leatherbacks", xscale=0.85, yscale=0.925, cex = 1)
				abline(v=seq(yr.low,yr.high)[first.yr])
				# Stochastic
				par(mar=c(4,4,1,1))
				plot(yr.low:(yr.low+9+first.yr), proj.med[2,1:(10+first.yr),keepers[3]], 
				     xlim=c((yr.low-0.25),(yr.low+10-0.75+first.yr)), 
				     ylim=c(0,max(proj.med[,1:(10+first.yr),keepers])),
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab="Annual Nesters", 
				     yaxs='i', xaxs='i')
				polygon(x = c(yr.low:(yr.low+9+first.yr), rev(yr.low:(yr.low+9+first.yr))), y = c(proj.med[1,1:(10+first.yr),keepers[3]], rev(proj.med[3,1:(10+first.yr),keepers[3]])), border=FALSE, col=col2rgbA("dodgerblue3", 0.5))
				polygon(x = c(yr.low:(yr.low+9+first.yr), rev(yr.low:(yr.low+9+first.yr))), y = c(proj.med[1,1:(10+first.yr),keepers[4]], rev(proj.med[3,1:(10+first.yr),keepers[4]])), border=FALSE, col=col2rgbA("chartreuse3", 0.5))
				lines(yr.low:(yr.low+9+first.yr), proj.med[2,1:(10+first.yr),keepers[3]], col="dodgerblue4", lwd=3)
				lines(yr.low:(yr.low+9+first.yr), proj.med[2,1:(10+first.yr),keepers[4]], col="chartreuse4", lwd=3)
				axis(1, yrpretty, yrpretty)
				mtext("Stochastic", side=3, font=3)
				fig.lab("W. Pac. Leatherbacks", xscale=0.85, yscale=0.925, cex = 1)
				abline(v=seq(yr.low,yr.high)[first.yr])
				par(mar=c(0,0,0,0))
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				legend("center", legend=c(expression(paste(N[j] -F," Median")), expression(paste(N[j] -F," 95%CI")),expression(paste(N[j]," Median")), expression(paste(N[j]," 95%CI"))), lty=c(1,NA,1,NA), lwd=c(3,NA,3,NA), pch=c(NA,15,NA,15), col=c("dodgerblue4",col2rgbA("dodgerblue3",0.5),"chartreuse4", col2rgbA("chartreuse3", 0.5)), bty='n', pt.cex = 2)

			dev.off()

			png(file=paste0(fig.path,spp,mode,"diff100.png"), units = 'in', res=300, width=7, height=7)
				###PROJECT 100
				layout(matrix(c(1,2,3,3),2,2), width=c(1,0.3))
				par(mar=c(3,4,2,1))
				# Deterministic
				plot(yr.low:yr.high, proj.diff.det[2,], 
				     xlim=c((yr.low-1),(yr.high+1)), 
				     ylim=range(pretty(range(proj.diff.det))),
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab=expression(paste(Delta," Annual Nesters")), 
				     yaxs='i', xaxs='i')
				polygon(x = c(yr.low:yr.high, rev(yr.low:yr.high)), y = c(proj.diff.det[1,], rev(proj.diff.det[3,])), border=FALSE, col=col2rgbA("darkorchid", 0.5))
				lines(yr.low:yr.high, proj.diff.det[2,], col="darkorchid4", lwd=3)
				yrpretty <- pretty(yr.low:(yr.high+1))[pretty(yr.low:(yr.high+1))>yr.low & pretty(yr.low:(yr.high+1)) < (yr.high+1)]
				axis(1, yrpretty, yrpretty)
				mtext("Deterministic", side=3, font=3)
				fig.lab("W. Pac. Leatherbacks", xscale=0.85, yscale=0.925, cex = 1)
				abline(v=seq(yr.low,yr.high)[first.yr])
				# Stochastic
				par(mar=c(4,4,1,1))
				plot(yr.low:yr.high, proj.diff.sto[2,], 
				     xlim=c((yr.low-1),(yr.high+1)), 
				     ylim=range(pretty(range(proj.diff.sto))),
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab=expression(paste(Delta," Annual Nesters")), 
				     yaxs='i', xaxs='i')
				polygon(x = c(yr.low:yr.high, rev(yr.low:yr.high)), y = c(proj.diff.sto[1,], rev(proj.diff.sto[3,])), border=FALSE, col=col2rgbA("darkorchid", 0.5))
				lines(yr.low:yr.high, proj.diff.sto[2,], col="darkorchid4", lwd=3)
				yrpretty <- pretty(yr.low:(yr.high+1))[pretty(yr.low:(yr.high+1))>yr.low & pretty(yr.low:(yr.high+1)) < (yr.high+1)]
				axis(1, yrpretty, yrpretty)
				mtext("Stochastic", side=3, font=3)
				fig.lab("W. Pac. Leatherbacks", xscale=0.85, yscale=0.925, cex = 1)
				abline(v=seq(yr.low,yr.high)[first.yr])
				par(mar=c(0,0,0,0))
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				legend("center", legend=c(expression(paste(N[j] - (N[j]-F)," Median")), expression(paste(N[j] - (N[j]-F)," 95%CI"))), lty=c(1,NA,1,NA), lwd=c(3,NA,3,NA), pch=c(NA,15,NA,15), col=c("darkorchid4",col2rgbA("darkorchid",0.5)), bty='n', pt.cex = 2, xpd=NA)

			dev.off()

			png(file=paste0(fig.path,spp,mode,"diff10.png"), units = 'in', res=300, width=7, height=7)
				#### PROJECT 10
				layout(matrix(c(1,2,3,3),2,2), width=c(1,0.3))
				par(mar=c(3,4,2,1))
				# Deterministic
				plot(yr.low:(yr.low+9+first.yr), proj.diff.det[2,1:(10+first.yr)], 
				     xlim=c((yr.low-0.25),(yr.low+10-0.75+first.yr)), 
				     ylim=range(pretty(range(proj.diff.det[,1:(10+first.yr)]))),
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab=expression(paste(Delta," Annual Nesters")), 
				     yaxs='i', xaxs='i')
				polygon(x = c(yr.low:(yr.low+9+first.yr), rev(yr.low:(yr.low+9+first.yr))), y = c(proj.diff.det[1,1:(10+first.yr)], rev(proj.diff.det[3,1:(10+first.yr)])), border=FALSE, col=col2rgbA("darkorchid", 0.5))
				lines(yr.low:(yr.low+9+first.yr), proj.diff.det[2,1:(10+first.yr)], col="darkorchid4", lwd=3)
				yrpretty <- pretty(yr.low:(yr.low+10+first.yr))[pretty(yr.low:(yr.low+10+first.yr))>yr.low & pretty(yr.low:(yr.low+10+first.yr)) < (yr.low+10+first.yr)]
				axis(1, yrpretty, yrpretty)
				mtext("Deterministic", side=3, font=3)
				fig.lab("W. Pac. Leatherbacks", xscale=0.85, yscale=0.925, cex = 1)
				abline(v=seq(yr.low,yr.high)[first.yr])
				# Stochastic
				par(mar=c(4,4,1,1))
				plot(yr.low:(yr.low+9+first.yr), proj.diff.sto[2,1:(10+first.yr)], 
				     xlim=c((yr.low-0.25),(yr.low+10-0.75+first.yr)), 
				     ylim=range(pretty(range(proj.diff.sto[,1:(10+first.yr)]))),
				     type='n', las=1, xlab="Years", 
				     xaxt='n', ylab=expression(paste(Delta," Annual Nesters")), 
				     yaxs='i', xaxs='i')
				polygon(x = c(yr.low:(yr.low+9+first.yr), rev(yr.low:(yr.low+9+first.yr))), y = c(proj.diff.sto[1,1:(10+first.yr)], rev(proj.diff.sto[3,1:(10+first.yr)])), border=FALSE, col=col2rgbA("darkorchid", 0.5))
				lines(yr.low:(yr.low+9+first.yr), proj.diff.sto[2,1:(10+first.yr)], col="darkorchid4", lwd=3)
				yrpretty <- pretty(yr.low:(yr.low+10+first.yr))[pretty(yr.low:(yr.low+10+first.yr))>yr.low & pretty(yr.low:(yr.low+10+first.yr)) < (yr.low+10+first.yr)]
				axis(1, yrpretty, yrpretty)
				mtext("Stochastic", side=3, font=3)
				fig.lab("W. Pac. Leatherbacks", xscale=0.85, yscale=0.925, cex = 1)
				abline(v=seq(yr.low,yr.high)[first.yr])
				par(mar=c(0,0,0,0))
				plot.new()
				plot.window(xlim=c(0,1),ylim=c(0,1))
				legend("center", legend=c(expression(paste(N[j] - (N[j]-F)," Median")), expression(paste(N[j] - (N[j]-F)," 95%CI"))), lty=c(1,NA,1,NA), lwd=c(3,NA,3,NA), pch=c(NA,15,NA,15), col=c("darkorchid4",col2rgbA("darkorchid",0.5)), bty='n', pt.cex = 2, xpd=NA)

			dev.off()

	
	return(list(tab.yr = tab.yr,
	            tab.prob = tab.prob))
}
cor.lab <- function(x, y, xscale, yscale, cex=1.2, adj=c(0.5,0.5)){
	xpt <- par('usr')[1]+abs(diff(par('usr')[1:2]))*xscale
	ypt <- par('usr')[3]+abs(diff(par('usr')[3:4]))*yscale

	cor <- round(cor(x,y),2)

	text(xpt, ypt, cor, xpd=NA, cex=cex, adj=adj)
}
kde.plot <- function(x, y, n=100){
	kde <- MASS::kde2d(x, y, n=n)
	kde$z <- kde$z/max(kde$z)
	contour(kde, levels=c(0.05, 0.5), col='black', drawlabels=FALSE, add=TRUE, lwd=c(1,3))
}
plot.joint <- function(trend, spp, graph){
	d.u <- density(trend$U, adj=2)
	d.N0 <- density(trend$N0, adj=2)
	d.Q <- density(trend$Q, adj=2)

	d.u$y <- d.u$y/max(d.u$y)
	d.N0$y <- d.N0$y/max(d.N0$y)
	d.Q$y <- d.Q$y/max(d.Q$y)

	samp <- sample(1:nrow(trend), floor(0.5*nrow(trend)))
	png(file=paste0(fig.path,spp,"joint_post.png"),width=7,height=7,units = 'in', res=300)
		par(mfrow=c(3,3), mar=c(4,4,1,1), cex.axis=1.1, cex.lab=1.1)
	for(i in 1:9){
		if(i ==1) {
			#1
			plot(d.u$x, d.u$y, type='l', xlab='r', ylab='Rel. Density', ylim=c(0,1), yaxt='n', las=1)
			axis(2, at=pretty(c(0,1)), las=1)
		}
		if(i==2){
			#2
			plot(trend$U, trend$N0, pch=".", col='grey80', xlab='r', ylab=expression(N[final]), las=1)
			kde.plot(trend$U[samp], trend$N0[samp], n=100)
		}
		if(i==3){
			#3
			plot(trend$U[samp], trend$Q[samp], pch=".", col='grey80',xlab='r', ylab='Q', las=1)
			kde.plot(trend$U, trend$Q, n=100)
		}
		if(i==4){
			#4
			plot.new()
			box()
			cor.lab(trend$U, trend$N0, xscale=0.5, yscale=0.5, cex=1.8)
			fig.lab(expression(paste("r, ",N[final])), xscale=0.25, yscale=0.1, cex=1.5)
		}
		if(i==5){
			#5
			plot(d.N0$x, d.N0$y, type='l', xlab=expression(N[final]), ylab='Rel. Density', ylim=c(0,1), yaxt='n', las=1)
			axis(2, at=pretty(c(0,1)), las=1)
		}
		if(i==6){
			#6
			plot(trend$N0[samp], trend$Q[samp], pch=".", col='grey80',xlab=expression(N[0]), ylab='Q', las=1)
			kde.plot(trend$N0, trend$Q, n=100)
		}
		if(i==7){
			#7
			plot.new()
			box()
			cor.lab(trend$U, trend$Q, xscale=0.5, yscale=0.5, cex=1.8)
			fig.lab("r, Q", xscale=0.25, yscale=0.1, cex=1.5)
		}
		if(i==8){
			#8
			plot.new()
			box()
			cor.lab(trend$N0, trend$Q, xscale=0.5, yscale=0.5, cex=1.8)
			fig.lab(expression(paste(N[final],", Q")), xscale=0.25, yscale=0.1, cex=1.5)
		}
		if(i==9){
			#9
			plot(d.Q$x, d.Q$y, type='l', xlab='Q', ylab='Rel. Density', ylim=c(0,1), yaxt='n', las=1)
			axis(2, at=pretty(c(0,1)), las=1)
		}
		if(i==graph) legend("topright", legend=c("0.05", "0.5"), lwd=c(1,3), bty='n')
	}
	dev.off()
}
mv.DC.init <- function(chain_id){
	lm <- lm(log(Len)~Year, data=DC.td.df)

	beta0 <- rnorm(1, coef(lm)[1], abs(coef(lm)[1]*0.1))
	beta1 <- rnorm(1, coef(lm)[2], abs(coef(lm)[2]*0.1))

	mu0 <- rnorm(1, 0, 1)

	sigma <- rtruncnorm(2, a=0, mean = apply(DC.td.dat$x, 2, sd), sd = apply(DC.td.dat$x, 2, sd)*0.1)
	rho <- rtruncnorm(1, a=-1, b=1, cor(DC.td.dat$x)[1,2], abs(cor(DC.td.dat$x)[1,2]*0.1))

	return(list(beta0 = beta0,
	            beta1 = beta1,
	            mu0 = mu0,
	            sigma = sigma,
	            rho = rho))
}
mv.DC.hurd.init <- function(chain_id){
	lm <- lm(log(Len)~Year, data=DC.td.df)

	betaL0 <- rnorm(1, coef(lm)[1], abs(coef(lm)[1]*0.1))
	betaL1 <- rnorm(1, coef(lm)[2], abs(coef(lm)[2]*0.1))

	lm <- lm(logit(M.mu)~log(Len), data=DC.td.df[DC.td.df$M.mu != 1,])

	betaD0 <- rnorm(1, coef(lm)[1], abs(coef(lm)[1]*0.1))
	betaD1 <- rnorm(1, coef(lm)[2], abs(coef(lm)[2]*0.1))

	DC.td.df$bin <- as.integer(DC.td.df$M.mu == 1)
	lm <- glm(bin~log(Len), data=DC.td.df, family=binomial)

	betaT0 <- rnorm(1, coef(lm)[1], abs(coef(lm)[1]*0.1))
	betaT1 <- rnorm(1, coef(lm)[2], abs(coef(lm)[2]*0.1))

	sigmaL <- rtruncnorm(1, a=0, 
	                     mean = sd(log(DC.td.df$Len),na.rm=T), 
	                     sd = sd(log(DC.td.df$Len),na.rm=T)*0.1)
	sigmaD <- rtruncnorm(1, a=0, 
                 mean = sd(logit(DC.td.df$M.mu[DC.td.df$M.mu != 1]),na.rm=T), 
                 sd = sd(logit(DC.td.df$M.mu[DC.td.df$M.mu != 1]),na.rm=T)*0.1)

	return(list(betaL0 = betaL0,
	            betaL1 = betaL1,
	            betaD0 = betaD0,
	            betaD1 = betaD1,
	            betaT0 = betaT0,
	            betaT1 = betaT1,
	            sigmaL = sigmaL,
	            sigmaD = sigmaD))
}
mv.CC.init <- function(chain_id){
	lm <- lm(log(Len)~Year, data=CC.td.df)

	beta0 <- rnorm(1, coef(lm)[1], abs(coef(lm)[1]*0.1))
	beta1 <- rnorm(1, coef(lm)[2], abs(coef(lm)[2]*0.1))

	mu0 <- rnorm(1, 0, 1)

	sigma <- rtruncnorm(2, a=0, mean = apply(CC.td.dat$x, 2, sd), sd = apply(CC.td.dat$x, 2, sd)*0.1)
	rho <- rtruncnorm(1, a=-1, b=1, cor(CC.td.dat$x)[1,2], abs(cor(CC.td.dat$x)[1,2]*0.1))

	return(list(beta0 = beta0,
	            beta1 = beta1,
	            mu0 = mu0,
	            sigma = sigma,
	            rho = rho))
}
plot.den <- function(x, adj=2, from=NULL, to=NULL, ...){
	if(is.null(from) & is.null(to)){
		d <- density(x, adj=2)
	}else if(is.null(from)){
		d <- density(x, adj=2, to=to)
	}else if(is.null(to)){
		d <- density(x, adj=2, from=from)
	}else{
		d <- density(x, adj=2, from=from, to=to)
	}
	
	d$y2 <- d$y/max(d$y)
	
	plot(d$x, d$y2, type='l', ...)
}
fig.lab <- function(i, xscale=0.05, yscale, cex=1.4){
	text(x = par('usr')[1] + abs(diff(par('usr')[1:2]))*xscale,
	     y = par('usr')[3] + abs(diff(par('usr')[3:4]))*yscale,
	     ifelse(is.numeric(i),LETTERS[i],i), xpd=NA, cex=cex)
}
poly.den <- function(x, adj=2, from=NULL, to=NULL, col.poly='black', alpha=0.05){
	if(is.null(from) & is.null(to)){
		d <- density(x, adj=2)
	}else if(is.null(from)){
		d <- density(x, adj=2, to=to)
	}else if(is.null(to)){
		d <- density(x, adj=2, from=from)
	}else{
		d <- density(x, adj=2, from=from, to=to)
	}

	d$y2 <- d$y/max(d$y)
	if(is.null(alpha)){
		polygon(c(d$x, 
		          rev(d$x)), 
		        c(rep(0,length(d$x)), 
		          rev(d$y2)), 
		        col=col.poly, border=FALSE)
	}else{
		q <- quantile(x, probs=c(alpha/2, (1-alpha/2)))

		xmin <- sapply(q, function(x) {which.min(abs(d$x - x))})

		polygon(c(d$x[xmin[1]:xmin[2]], 
		          rev(d$x[xmin[1]:xmin[2]])), 
		        c(rep(0,length(d$x[xmin[1]:xmin[2]])), 
		          rev(d$y2[xmin[1]:xmin[2]])), 
		        col=col.poly, border=FALSE)
	}
	
}
poly.den.partial <- function(x, adj=2, col.poly='black', prop=1, from=0,  to=1, xpart=c(0,1), ypart=c(0,1), alpha=NULL){
	d <- density(x, adj=adj, from=from, to=to)
	d$y2 <- d$y/max(d$y) * prop
	d$x2 <- scales::rescale(d$x, from=c(from,to), to=xpart)
	d$y3 <- scales::rescale(d$y2, from=c(0,1), to=ypart)

	if(!is.null(alpha)){
		d$x2[d$y2<alpha] <- NA
		d$y3[d$y2<alpha] <- NA
		d$x2 <- na.omit(d$x2)
		d$y3 <- na.omit(d$y3)
	}
	polygon(c(d$x2, rev(d$x2)), c(rep(min(ypart),length(d$x2)), rev(d$y3)), col=col.poly, border=FALSE)
}
logit <- function(x){
	log(x/(1-x))
}
inv.logit <- function(x){
	exp(x)/(exp(x)+1)
}
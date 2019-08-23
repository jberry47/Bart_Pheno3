################################################
# Author: Dijun Chen (chend@ipk-gatersleben.de)
# Update on June 2, 2014
################################################
cat("loading modeling.models.R \n")

logistic.model <- function(data.y, data.time, pred.time=NULL, 
				logistic.optimize.step=0.01, col="#E4007F", lty=6, ...){
	kmax <- max(data.y, na.rm=TRUE) + 0.01
	##if(kmax<0.4) {kmax=0.4}
	myFunc <- function(data.y, kmax, data.time){
		kmaxn <- kmax
		mlogistic <- lm(log(data.y/(kmax-data.y)) ~ data.time)
		rold <- summary(mlogistic)$adj.r.squared
		for (kmax1 in seq(kmax, 1.5, by=logistic.optimize.step)){
			mlogistic <- lm(log(data.y/(kmax1-data.y)) ~ data.time)
			if (summary(mlogistic)$adj.r.squared > rold){
				rold <- summary(mlogistic)$adj.r.squared
				kmaxn <- kmax1
			}
		}
		kmax <- kmaxn
		mlogistic <- lm(log(data.y/(kmax-data.y)) ~ data.time)
		## summary(mlogistic)
		return(list(kmax=kmax, model=mlogistic))
	}
	
	myTry <- try(myFunc(data.y, kmax, data.time), silent=TRUE)
	if (class(myTry) == "try-error") {
		cat(paste("\t>>> Modeling error\n"))
		cat(myTry)
		return(NULL)
	}else{
		kmax <- myTry$kmax
		mlogistic <- myTry$model
		inflection <- abs(summary(mlogistic)$coefficients[1])/summary(mlogistic)$coefficients[2]
		intrinsic.growth.rate <- summary(mlogistic)$coefficients[2]
		inflection.growth <- kmax/(1+exp(-predict(mlogistic, data.frame(data.time=inflection))))
		inflection.growth.rate <- intrinsic.growth.rate*inflection.growth*(1-inflection.growth/kmax)
		inflection.relative.growth.rate <- inflection.growth.rate / inflection.growth
		mlogistic <- lm(log(data.y/(kmax-data.y)) ~ data.time)
		if(is.null(pred.time)){pred.time <- data.time}
		rsquared <- summary(mlogistic)$adj.r.squared
		prediction <- kmax/(1+exp(-predict(mlogistic, data.frame(data.time=pred.time))))
		legend <- substitute(expression(paste(italic(K*y[0]/(y[0]+(K-y[0])*e^-rt)), 
					", Logistic, ", italic(R)^2 == RR)), list(RR=sprintf("%.4f", rsquared)))[2]
		list(model=mlogistic, kmax=kmax, prediction=prediction, rsquared=rsquared, 
				inflection=inflection, inflection.growth=inflection.growth, inflection.growth.rate=inflection.growth.rate,
				inflection.relative.growth.rate=inflection.relative.growth.rate, intrinsic.growth=intrinsic.growth.rate,
				legend=legend, color=col, lty=lty)
	}
}

gompetz.model <- function(data.y, data.time, pred.time=NULL, kmax=NULL, col="#F39800", lty=5, ...){
	if(is.null(kmax)) kmax <- max(data.y, na.rm=TRUE) + 0.01
	if(is.null(pred.time)){pred.time <- data.time}
	mgompetz <- lm(-log(-log(data.y/kmax)) ~ data.time)
	rsquared <- summary(mgompetz)$adj.r.squared
	prediction <- kmax*exp(-exp(-predict(mgompetz, data.frame(data.time=pred.time))))
	legend <- substitute(expression(paste(italic(K*e^ln(y[0]/K*e^(-rt))), ", Gompetz, ", 
				italic(R)^2 == RR)), list(RR=sprintf("%.4f", rsquared)))[2]
	list(model=mgompetz, kmax=kmax, prediction=prediction, rsquared=rsquared, 
			legend=legend, color=col, lty=lty)
}

monomolecular.model <- function(data.y, data.time, pred.time=NULL, kmax=NULL, col="#7E318E", lty=4, ...){
	if(is.null(kmax)) kmax <- max(data.y, na.rm=TRUE) + 0.01
	if(is.null(pred.time)){pred.time <- data.time}
	mmonomolecular <- lm(log(1/(kmax-data.y)) ~ data.time)
	rsquared <- summary(mmonomolecular)$adj.r.squared
	prediction <- (kmax-exp(-predict(mmonomolecular, data.frame(data.time=pred.time))))
	legend <- substitute(expression(paste(italic(y[0]-(K-y[0])*e^-rt), ", Monomolecular, ", 
				italic(R)^2 == RR)), list(RR=sprintf("%.4f", rsquared)))[2]
	list(model=mmonomolecular, kmax=kmax, prediction=prediction, rsquared=rsquared, 
			legend=legend, color=col, lty=lty)
}

exponential.model <- function(data.y, data.time, pred.time=NULL, col="#B28247", lty=3, ...){
	if(is.null(pred.time)){pred.time <- data.time}
	mexponential <- lm(log(data.y) ~ data.time)
	rsquared <- summary(mexponential)$adj.r.squared
	prediction <- exp(predict(mexponential, data.frame(data.time=pred.time)))
	legend <- substitute(expression(paste(italic(y[0]*e^rt), ", Exponential, ", 
				italic(R)^2 == RR)), list(RR=sprintf("%.4f", rsquared)))[2]
	list(model=mexponential, prediction=prediction, rsquared=rsquared, 
			legend=legend, color=col, lty=lty)
}

linear.model <- function(data.y, data.time, pred.time=NULL, col="#00A29A", lty=2, ...){
	if(is.null(pred.time)){pred.time <- data.time}
	mlinear <- lm(data.y ~ data.time)
	rsquared <- summary(mlinear)$adj.r.squared
	growth.rate <- summary(mlinear)$coefficients[2,1]
	prediction <- predict(mlinear, data.frame(data.time=pred.time))
	legend <- substitute(expression(paste(italic(y[0]+r*t), ", Linear, ", 
				italic(R)^2 == RR)), list(RR=sprintf("%.4f", rsquared)))[2]
	list(model=mlinear, prediction=prediction, rsquared=rsquared, growth.rate=growth.rate, 
			initial.relative.growth.rate=growth.rate/prediction[1], 
			legend=legend, color=col, lty=2)
}

quadratic.model <- function(data.y, data.time, pred.time=NULL, col="#00A29A", lty=2, ...){
	if(is.null(pred.time)){pred.time <- data.time}
	data.time <- data.time+data.time^2
	mquadratic <- lm(data.y ~ data.time)
	rsquared <- summary(mquadratic)$adj.r.squared
	prediction <- predict(mquadratic, data.frame(data.time=pred.time+pred.time^2))
	legend <- substitute(expression(paste(italic(a*t+b*t^2), ", Quadratic, ", 
				italic(R)^2 == RR)), list(RR=sprintf("%.4f", rsquared)))[2]
	list(model=mquadratic, prediction=prediction, rsquared=rsquared, 
			legend=legend, color=col, lty=lty)
}

bellshape1.model <- function(data.y, data.time, pred.time=NULL, col="#B28247", lty=3, ...){
	if(is.null(pred.time)){pred.time <- data.time}
	t.max <- data.time[which(data.y==max(data.y, na.rm=TRUE))][1]
	data.time <- (data.time-t.max)^2
	mbellshape1 <- lm(log(data.y) ~ data.time)
	rsquared <- summary(mbellshape1)$adj.r.squared
	prediction <- exp(predict(mbellshape1, data.frame(data.time=(pred.time-t.max)^2)))
	legend <- substitute(expression(paste(italic(A*e^{a*(t-{t[max]})^2}), ", Bell-shaped 1, ", 
				italic(R)^2 == RR)), list(RR=sprintf("%.4f", rsquared)))[2]
	list(model=mbellshape1, prediction=prediction, rsquared=rsquared, 
			legend=legend, color=col, lty=lty)
}

bellshape2.model <- function(data.y, data.time, pred.time=NULL, col="#7E318E", lty=4, ...){
	if(is.null(pred.time)){pred.time <- data.time}
	data.time.lg <- log(data.time)
	mbellshape2 <- lm(log(data.y) ~ data.time.lg+data.time)
	rsquared <- summary(mbellshape2)$adj.r.squared
	t.max <- abs(summary(mbellshape2)$coefficients[2,1]/summary(mbellshape2)$coefficients[3,1])
	inflection.before <- abs((summary(mbellshape2)$coefficients[2,1]-
				sqrt(summary(mbellshape2)$coefficients[2,1]))/(summary(mbellshape2)$coefficients[3,1]))
	inflection.after <- abs((summary(mbellshape2)$coefficients[2,1]+
				sqrt(summary(mbellshape2)$coefficients[2,1]))/(summary(mbellshape2)$coefficients[3,1]))
	prediction <- exp(predict(mbellshape2, data.frame(data.time=pred.time, data.time.lg=log(pred.time))))
	IPs <- c(t.max, inflection.before, inflection.after)
	pred.infl <- exp(predict(mbellshape2, data.frame(data.time=IPs, data.time.lg=log(IPs))))
	legend <- substitute(expression(paste(italic(A*t^b*e^{-a*t}), ", Bell-shaped 2, ", 
				italic(R)^2 == RR)), list(RR=sprintf("%.4f", rsquared)))[2]
	list(model=mbellshape2, prediction=prediction, rsquared=rsquared, t.max=t.max, 
			inflection.before=inflection.before, inflection.after=inflection.after, 
			inflection.before.pred=pred.infl[2], inflection.after.pred=pred.infl[3], 
			t.max.pred=pred.infl[1], legend=legend, color=col, lty=lty)
}

bellshape3.model <- function(data.y, data.time, pred.time=NULL, col="#F39800", lty=5, ...){
	if(is.null(pred.time)){pred.time <- data.time}
	data.time.sq <- data.time^2
	mbellshape3 <- lm(log(data.y) ~ data.time+data.time.sq)
	rsquared <- summary(mbellshape3)$adj.r.squared
	t.max <- abs(summary(mbellshape3)$coefficients[2,1]/(2*summary(mbellshape3)$coefficients[3,1]))
	inflection.before <- abs((-summary(mbellshape3)$coefficients[2,1]+
				sqrt(0-2*summary(mbellshape3)$coefficients[3,1]))/(2*summary(mbellshape3)$coefficients[3,1]))
	inflection.after <- abs((-summary(mbellshape3)$coefficients[2,1]-
				sqrt(0-2*summary(mbellshape3)$coefficients[3,1]))/(2*summary(mbellshape3)$coefficients[3,1]))
	prediction <- exp(predict(mbellshape3, data.frame(data.time=pred.time, data.time.sq=pred.time^2)))
	IPs <- c(t.max, inflection.before, inflection.after)
	pred.infl <- exp(predict(mbellshape3, data.frame(data.time=IPs, data.time.sq=IPs^2)))
	legend <- substitute(expression(paste(italic(A*e^{b*t-a*t^2}), ", Bell-shaped 3, ", 
				italic(R)^2 == RR)), list(RR=sprintf("%.4f", rsquared)))[2]
	list(model=mbellshape3, prediction=prediction, rsquared=rsquared, t.max=t.max, 
			inflection.before=inflection.before, inflection.after=inflection.after, 
			inflection.before.pred=pred.infl[2], inflection.after.pred=pred.infl[3], 
			t.max.pred=pred.infl[1], legend=legend, color=col, lty=lty)
}

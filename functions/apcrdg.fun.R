apcrdg <-
function(r, nrisk = matrix(1, nrow(r), ncol(r)), c1 = 1, c2 = 1, cc = 1, apc = 1, 
	intercept = T, rr = F, lambda = 0, fam = "pois", gamma = 2, center = F, 
	Plot = F, Pltscale = F, tnull = 0, CI = T, smoothplot = F, expplot = F, 
	panelplot = F)
{
## Compute the intrinsic estimator or an arbitrary estimator of APC model
## SEE FU (1998). 
## 
##  if lambda = 0, an arbitrary estimator with an arbitrary constraint is computed;
##            > 0, the ridge estimator is computed with ridge trace estimator 
##				(good approximation to the intrinsic estimator)
##				for small lambda, such as lambda = 1e-4;
##		  < 0, the intrinsic estimator is computed by removing the null vector 
##				from the solution vector.
##			 The SE is computed via principal component method.
##
## 	c1, c2: column number for the constraint columns.
##	cc : the coef for the constraint, Beta2 = cc * Beta1.
	a <- nrow(r)
	p <- ncol(r)
	y <- c(t(r * nrisk))
	offst <- c(t(nrisk))
	x <- apcmat(a, p)
	sol <- vector("list")
	if(lambda <= 0) {
		if(lambda < 0) {
			c1 <- 1
			c2 <- 2
			cc <- 1
			apc <- 1
		}
		xc <- constr(x, a, p, c1, c2, cc, apc)[, -1]
		if(fam == "pois")
			fit <- glm(y ~ xc + offset(log(offst)), family = poisson(
				link = log), intercept = intercept, control = 
				glm.control(epsilon = 1e-007))
		if(fam == "qlik")
			fit <- glm(y ~ xc + offset(log(offst)), family = quasi(
				link = log, var = "mu"), intercept = intercept, 
				control = glm.control(epsilon = 1e-007))
		sol$dev <- fit$dev
		sol$df <- fit$df
		sol$inter <- fit$coef[1]
		if(apc == 1) {
			aa <- fit$coef[2:(a - 1)]
			aa <- c(aa[1:(c2 - 1)], cc * aa[c1], aa[-1: - (c2 - 1)])
			pp <- fit$coef[a:(a + p - 2)]
			coh <- fit$coef[(a + p - 1):length(fit$coef)]
		}
		if(apc == 2) {
			aa <- fit$coef[2:a]
			pp <- fit$coef[(a + 1):(a + p - 2)]
			pp <- c(pp[1:(c2 - 1)], cc * pp[c1], pp[-1: - (c2 - 1)])
			coh <- fit$coef[(a + p - 1):length(fit$coef)]
		}
		if(apc == 3) {
			aa <- fit$coef[2:a]
			pp <- fit$coef[(a + 1):(a + p - 1)]
			coh <- fit$coef[(a + p):length(fit$coef)]
			coh <- c(coh[1:(c2 - 1)], cc * coh[c1], coh[-1: - (c2 - 1
				)])
		}
		if(lambda < 0) {
###### intercept IS 1st column of X to remove the null eigen vector
#	 yields same estimator as the other case below.
#	eigvec0 <- eigen(t(x)%*%x)$vec[,ncol(x)]
#	fitvec <- c(fit$coef[1:2],cc*fit$coef[2],fit$coef[c(-1,-2)])
#	fitvec0 <- fitvec - eigvec0*c(crossprod(eigvec0,fitvec))
#	aa <- fitvec0[2:(a)]
#	pp <- fitvec0[(a+1):(a+p-1)]
#	coh <- fitvec0[-c(1:(a+p-1))]
###### intercept IS NOT 1st column of X to remove the null eigen vector
			eigvec0 <- eigen(t(x[, -1]) %*% x[, -1])$vec[, ncol(x) - 
				1]
			fitvec <- c(fit$coef[2], cc * fit$coef[2], fit$coef[c(-1, 
				-2
				)])
			fitvec0 <- fitvec - eigvec0 * c(crossprod(eigvec0, fitvec
				))
			aa <- fitvec0[1:(a - 1)]
			pp <- fitvec0[a:(a + p - 2)]
			coh <- fitvec0[ - c(1:(a + p - 2))]
		}
		sol$aa <- c(aa,  - sum(aa))
		sol$pp <- c(pp,  - sum(pp))
		sol$coh <- c(coh,  - sum(coh))
		if(Plot == T) {
			if(Pltscale == F) {
				plot(c(aa,  - sum(aa)), type = "b", xlab = 
				  "age group", ylab = "aa", main = "Age trend")
				plot(c(pp,  - sum(pp)), type = "b", xlab = 
				  "period group", ylab = "pp", main = 
				  "Period trend")
				plot(c(coh,  - sum(coh)), type = "b", xlab = 
				  "cohort group", ylab = "cc", main = 
				  "Cohort trend")
			}
			if(Pltscale == T) {
				pmin <- min(c(aa,  - sum(aa), pp,  - sum(pp), coh,
				   - sum(coh)))
				pmax <- max(c(aa,  - sum(aa), pp,  - sum(pp), coh,
				   - sum(coh)))
				plot(c(aa,  - sum(aa)), type = "b", xlab = 
				  "age group", ylab = "aa", ylim = c(pmin, pmax), 
				  main = "Age trend")
				plot(c(pp,  - sum(pp)), type = "b", xlab = 
				  "period group", ylab = "pp", ylim = c(pmin, 
				  pmax), main = "Period trend")
				plot(c(coh,  - sum(coh)), type = "b", xlab = 
				  "cohort group", ylab = "cc", ylim = c(pmin, 
				  pmax), main = "Cohort trend")
			}
		}
##### Need to correct the var by incorporating
##### the constraint into var/covar matrix
##### SE calculation
		if(lambda == 0) {
			fitstd <- summary.glm(fit)$coef
			p.value <- 2 * (1 - pt(abs(fitstd[, 3]), fit$df))
			std0 <- cbind(fitstd, p.value)
			cor.a <- summary.glm(fit)$corr[2:a, 2:a]	
	####		 matrix(-1/(a-1),a-1,a-1) + (1+1/(a-1))*diag(a-1)
			std.a <- fitstd[2:a, 2]	
	####        cov.aa <- t(as.matrix(rep(1,a-1)))%*%diag(std.a)%*%cor.a%*%diag(std.a)%*%
####			as.matrix(rep(1,a-1))
			aa.std <- sqrt(sum(diag(std.a) %*% cor.a %*% diag(std.a))
				)
			cor.p <- summary.glm(fit)$corr[(a + 1):(a + p - 1), (a + 
				1):(a + p - 1)]	
	####		matrix(-1/(p-1),p-1,p-1) + (1+1/(p-1))*diag(p-1)
			std.p <- fitstd[c((a + 1):(a + p - 1)), 2]	
	####        cov.pp <- t(as.matrix(rep(1,p-1)))%*%diag(std.p)%*%cor.p%*%diag(std.p)%*%
####			as.matrix(rep(1,p-1))
			pp.std <- sqrt(sum(diag(std.p) %*% cor.p %*% diag(std.p))
				)
			cor.c <- summary.glm(fit)$corr[c(-1: - (a + p - 1)), c(-1:
				 - (a + p - 1))]	
	####		matrix(-1/(a+p-2),a+p-2,a+p-2) + (1+1/(a+p-2))*diag(a+p-2)
			std.c <- fitstd[c(-1: - (a + p - 1)), 2]	
	####        cov.cc <- t(as.matrix(rep(1,a+p-2)))%*%diag(std.c)%*%cor.c%*%diag(std.c)%*%
####			as.matrix(rep(1,a+p-2))
			cc.std <- sqrt(sum(diag(std.c) %*% cor.c %*% diag(std.c))
				)
			if(apc == 1) {
				std.a <- std0[2:(a - 1),  ]
				if(c2 <= (a - 2))
				  std.a <- rbind(std.a[1:(c2 - 1),  ], c(std.a[c1,
				    1] * cc, std.a[c1, 2] * abs(cc), sign(cc) * 
				    std.a[c1, 3], std.a[c1, 4]), std.a[c2:(a - 2),
				    ])
				else std.a <- rbind(std.a[1:(c2 - 1),  ], c(std.a[
				    c1, 1] * cc, std.a[c1, 2] * abs(cc), sign(cc) *
				    std.a[c1, 3], std.a[c1, 4]))
				cor.a <- summary.glm(fit)$corr[2:(a - 1), 2:(a - 
				  1)]
				std.a0 <- std.a[, 2][ - c2]
				ccol <- rep(1, a - 2)
				ccol[c1] <- ccol[c1] * (1 + cc)
				aa.std <- sqrt(t(as.matrix(ccol)) %*% diag(std.a0
				  ) %*% cor.a %*% diag(std.a0) %*% as.matrix(ccol
				  ))
				std.a <- rbind(std.a, c( - sum(std.a[, 1]), 
				  aa.std,  - sum(std.a[, 1])/aa.std, 2 * (1 - pt(
				  abs(sum(std.a[, 1])/aa.std), fit$df))))
				std.p <- std0[a:(a + p - 2),  ]
				cor.p <- summary.glm(fit)$corr[a:(a + p - 2), a:(
				  a + p - 2)]
				std.p0 <- std.p[, 2]
				pp.std <- sqrt(t(as.matrix(rep(1, p - 1))) %*% 
				  diag(std.p0) %*% cor.p %*% diag(std.p0) %*% 
				  as.matrix(rep(1, p - 1)))
				std.p <- rbind(std.p, c( - sum(std.p[, 1]), 
				  pp.std,  - sum(std.p[, 1])/pp.std, 2 * (1 - pt(
				  abs(sum(std.p[, 1])/pp.std), fit$df))))
				std.c <- std0[-1: - (a + p - 2),  ]
				cor.c <- summary.glm(fit)$corr[-1: - (a + p - 2), 
				  -1
				  : - (a + p - 2)]
				std.c0 <- std.c[, 2]
				cc.std <- sqrt(t(as.matrix(rep(1, a + p - 2))) %*% 
				  diag(std.c0) %*% cor.c %*% diag(std.c0) %*% 
				  as.matrix(rep(1, a + p - 2)))
				std.c <- rbind(std.c, c( - sum(std.c[, 1]), 
				  cc.std,  - sum(std.c[, 1])/cc.std, 2 * (1 - pt(
				  abs(sum(std.c[, 1])/cc.std), fit$df))))
				sol$std <- rbind(std0[1,  ], std.a, std.p, std.c)
		#		cov.aa0 <- ( t(as.matrix(rep(1,a-1)))%*%diag(std.a)%*%cor.a%*%diag(std.a)%*%
#			as.matrix(rep(1,a-1)) )[1,]
#		cov.aa <- cov.aa + cc*sum(cov.aa0) + cc^2*cov.aa0[1]
#		aa.std<-sqrt(cov.aa)
			}
			if(apc == 2) {
#				cov.pp0 <- (t(as.matrix(rep(1, p - 1))) %*%
#				  diag(std.p) %*% cor.p %*% diag(std.p) %*%
#				  as.matrix(rep(1, p - 1)))[1,  ]
#				cov.pp <- cov.pp + cc * sum(cov.pp0) + 
#				  cc^2 * cov.pp0[1]
#				pp.std <- sqrt(cov.pp)
			}
			if(apc == 3) {
#				cov.cc0 <- (t(as.matrix(rep(1, length(
#				  std.c)))) %*% diag(std.c) %*% cor.c %*% 
#				  diag(std.c) %*% as.matrix(rep(1, 
#				  length(std.c))))[1,  ]
#				cov.cc <- cov.cc + cc * sum(cov.cc0) + 
#				  cc^2 * cov.cc0[1]
#				cc.std <- sqrt(cov.cc)
			}
			sol$aa.std <- c( - sum(aa), aa.std,  - sum(aa)/aa.std, 2 * (
				1 - pt(abs(sum(aa))/aa.std, fit$df)))
			sol$pp.std <- c( - sum(pp), pp.std,  - sum(pp)/pp.std, 2 * (
				1 - pt(abs(sum(pp))/pp.std, fit$df)))
			sol$cc.std <- c( - sum(coh), cc.std,  - sum(coh)/cc.std, 
				2 * (1 - pt(abs(sum(coh))/cc.std, fit$df)))
			if(rr == T) {
				if(apc == 1)
				  rr <- rbind(fitstd[1:(a - 1), 1:2], c( - sum(aa
				    ), aa.std), fitstd[c((a):(a + p - 2)), 1:2], 
				    c( - sum(pp), pp.std), fitstd[c(-1: - (a + p - 
				    2)), 1:2], c( - sum(coh), cc.std))
				if(apc == 2)
				  rr <- rbind(fitstd[1:a, 1:2], c( - sum(aa), 
				    aa.std), fitstd[c((a + 1):(a + p - 2)), 1:2], 
				    c( - sum(pp), pp.std), fitstd[c(-1: - (a + p - 
				    2)), 1:2], c( - sum(coh), cc.std))
				if(apc == 3)
				  rr <- rbind(fitstd[1:a, 1:2], c( - sum(aa), 
				    aa.std), fitstd[c((a + 1):(a + p - 1)), 1:2], 
				    c( - sum(pp), pp.std), fitstd[c(-1: - (a + p - 
				    1)), 1:2], c( - sum(coh), cc.std))
				Expon <- exp(rr[, 1])
				CI.L <- exp(rr[, 1] - 1.96 * rr[, 2])
				CI.R <- exp(rr[, 1] + 1.96 * rr[, 2])
				sol$rr <- cbind(rr, Expon, CI.L, CI.R)
			}
		}
## lambda = 0
	}
	else {
		fit <- as.vector(glmlogcol(x = x[, -1], y = y, lam = lambda, gam
			 = gamma, center = center, tnull = tnull)$coef.bridge)
		sol$intercept <- fit[1]
		aa <- c(fit[2:(a)],  - sum(fit[2:(a)]))
		pp <- c(fit[(a + 1):(a + p - 1)],  - sum(fit[(a + 1):(a + p - 1)]
			))
		coh <- c(fit[ - (1:(a - 1 + p))],  - sum(fit[ - (1:(a - 1 + p))])
			)
		if(Plot == T) {
			if(smoothplot == F) {
				par(mfrow = c(2, 2))
				plot(aa, type = "b", xlab = "age group", ylab = 
				  "aa", main = "Age trend")
				plot(pp, type = "b", xlab = "period group", ylab
				   = "pp", main = "Period trend")
				par(mfrow = c(2, 1), new = F, fig = c(0, 1, 0, 
				  0.5))
				plot(coh, type = "b", xlab = "cohort group", ylab
				   = "cc", main = "Cohort trend")	
	#	par(mfrow=c(1,1))
			}
			if(smoothplot == T) {
# plot age effect
				fit.spl <- smooth.spline(c(1:a), aa, df = 4)	
	# smoothing spline
				res <- (fit.spl$yin - fit.spl$y)/(1 - fit.spl$lev
				  )	# jackknife residuals
				sigma <- sqrt(var(res))	# estimated sd
				upper <- fit.spl$y + 1.96 * sigma * sqrt(fit.spl$
				  lev)	# upper 95% conf. band
				lower <- fit.spl$y - 1.96 * sigma * sqrt(fit.spl$
				  lev)	# lower 95% conf. band
				ymin <- min(aa, fit.spl$y, lower)
				ymax <- max(aa, fit.spl$y, upper)	
	#	matplot(fit$x, cbind(upper, fit.spl$y, lower), type="plp", pch=".")
				plot(aa, type = "b", xlab = "age group", main = 
				  "Age trend", ylim = c(ymin, ymax))
				if(CI == T) {
				  lines(fit.spl$x, upper, lty = 4, col = 3)
				  lines(fit.spl$x, lower, lty = 4, col = 6)
				  lines(fit.spl$x, fit.spl$y, lty = 2, col = 2)
				}
# plot period effect
				fit.spl <- smooth.spline(c(1:p), pp, df = 4)
				res <- (fit.spl$yin - fit.spl$y)/(1 - fit.spl$lev
				  )
				sigma <- sqrt(var(res))
				upper <- fit.spl$y + 1.96 * sigma * sqrt(fit.spl$
				  lev)
				lower <- fit.spl$y - 1.96 * sigma * sqrt(fit.spl$
				  lev)
				ymin <- min(pp, fit.spl$y, lower)
				ymax <- max(pp, fit.spl$y, upper)
				plot(pp, type = "b", xlab = "period group", main
				   = "Period trend", ylim = c(ymin, ymax))
				if(CI == T) {
				  lines(fit.spl$x, fit.spl$y, lty = 2, col = 2)
				  lines(fit.spl$x, upper, lty = 4, col = 3)
				  lines(fit.spl$x, lower, lty = 4, col = 6)
				}
# plot cohort effect
				fit.spl <- smooth.spline(c(1:(a + p - 1)), coh, 
				  df = 4)
				res <- (fit.spl$yin - fit.spl$y)/(1 - fit.spl$lev
				  )
				sigma <- sqrt(var(res))
				upper <- fit.spl$y + 1.96 * sigma * sqrt(fit.spl$
				  lev)
				lower <- fit.spl$y - 1.96 * sigma * sqrt(fit.spl$
				  lev)
				ymin <- min(coh, fit.spl$y, lower)
				ymax <- max(coh, fit.spl$y, upper)
				plot(coh, type = "b", xlab = "cohort group", main
				   = "Cohort trend", ylim = c(ymin, ymax))
				if(CI == T) {
				  lines(fit.spl$x, upper, lty = 4, col = 3)
				  lines(fit.spl$x, lower, lty = 4, col = 6)
				  lines(fit.spl$x, fit.spl$y, lty = 2, col = 2)
				}
				if(expplot == T) {
				  fit.spl <- smooth.spline(c(1:a), exp(aa), df = 
				    4)
				  res <- (fit.spl$yin - fit.spl$y)/(1 - fit.spl$
				    lev)
				  sigma <- sqrt(var(res))
				  upper <- fit.spl$y + 1.96 * sigma * sqrt(
				    fit.spl$lev)
				  lower <- fit.spl$y - 1.96 * sigma * sqrt(
				    fit.spl$lev)
				  ymin <- min(exp(aa), fit.spl$y, lower)
				  ymax <- max(exp(aa), fit.spl$y, upper)
				  plot(exp(aa), type = "b", xlab = "age group", 
				    ylab = "Exp scale", main = "Age trend", ylim
				     = c(ymin, ymax))
				  if(CI == T) {
				    lines(fit.spl$x, upper, lty = 4, col = 3)
				    lines(fit.spl$x, lower, lty = 4, col = 6)
				    lines(fit.spl$x, fit.spl$y, lty = 2, col = 2)
				  }
#	plot(exp(aa),type="b",xlab="age group",ylab="Exp scale", main="Age trend",ylim=c(ymin,ymax))
				  fit.spl <- smooth.spline(c(1:p), exp(pp), df = 
				    4)
				  res <- (fit.spl$yin - fit.spl$y)/(1 - fit.spl$
				    lev)
				  sigma <- sqrt(var(res))
				  upper <- fit.spl$y + 1.96 * sigma * sqrt(
				    fit.spl$lev)
				  lower <- fit.spl$y - 1.96 * sigma * sqrt(
				    fit.spl$lev)
				  ymin <- min(exp(pp), fit.spl$y, lower)
				  ymax <- max(exp(pp), fit.spl$y, upper)
				  plot(exp(pp), type = "b", xlab = "period group",
				    ylab = "Exp scale", main = "Period trend", 
				    ylim = c(ymin, ymax))
				  if(CI == T) {
				    lines(fit.spl$x, upper, lty = 4, col = 3)
				    lines(fit.spl$x, lower, lty = 4, col = 6)
				    lines(fit.spl$x, fit.spl$y, lty = 2, col = 2)
				  }
#	plot(exp(pp),type="b",xlab="period group", ylab="Exp scale", main="Period trend")
				  fit.spl <- smooth.spline(c(1:(a + p - 1)), exp(
				    coh), df = 4)
				  res <- (fit.spl$yin - fit.spl$y)/(1 - fit.spl$
				    lev)
				  sigma <- sqrt(var(res))
				  upper <- fit.spl$y + 1.96 * sigma * sqrt(
				    fit.spl$lev)
				  lower <- fit.spl$y - 1.96 * sigma * sqrt(
				    fit.spl$lev)
				  ymin <- min(exp(coh), fit.spl$y, lower)
				  ymax <- max(exp(coh), fit.spl$y, upper)
				  plot(exp(coh), type = "b", xlab = 
				    "cohort group", ylab = "Exp scale", main = 
				    "Cohort trend", ylim = c(ymin, ymax))
				  if(CI == T) {
				    lines(fit.spl$x, upper, lty = 4, col = 3)
				    lines(fit.spl$x, lower, lty = 4, col = 6)
				    lines(fit.spl$x, fit.spl$y, lty = 2, col = 2)
				  }
				}
			}
#	plot(exp(coh),type="b",xlab="cohort group", ylab="Exp scale", main="Cohort trend")
		}
		sol$aa <- aa
		sol$pp <- pp
		sol$coh <- coh
		sol$lambda <- lambda
		sol$gamma <- gamma
	}
	if(panelplot == T) {
		par(mfrow = c(2, 2))
		plot(c(aa,  - sum(aa)), type = "b", xlab = "age group", main = 
			"Age trend", ylab = "")
		plot(c(pp,  - sum(pp)), type = "b", xlab = "period group", main
			 = "Period trend", ylab = "")
		par(mfrow = c(2, 1), new = F, fig = c(0, 1, 0, 0.5))
		plot(c(coh,  - sum(coh)), type = "b", xlab = "cohort group", main
			 = "Cohort trend", ylab = "")
		par(mfrow = c(1, 1))
	}
## lambda <= 0
	sol
}

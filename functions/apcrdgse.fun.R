apcrdgse <- 
function(r, nrisk = matrix(1, nrow(r), ncol(r)), princomp = 1, lambda = 0, Plot
	 = F, offst = F, stderr = F, CIplot = F, pylim = c(0, 0), ieff = 2, 
	center = F, fam = "pois", amin = 1, pmin = 1, cmin = 1, gapyear = 1, 
	niter = 10)
{
############# Different models controlled by lambda value ###############
## fit APC model to a matrix of rates using log-linear model. 
## if lambda >0, the solution is the ridge penalty solution, which 
## converges to the ridge trace estimator, the intrinsic estimator
## of the APC model. (See FU (1998)).
## if lambda < 0, the estimator is the intrinsic estimator by removing the null
##  			vector from the estimator obtained with a constraint.
################# Estimation of Stderr ##############
## if princomp = 1, principal component method is used to compute the standard errors
## of the intrinsic estimator. 
##
## if princomp = 0, use constraint method as follows.
## ieff is the group number for the extra constraint. It is for the computation 
## of the standard errors. Diff ieff yields diff stderr since the condition 
## number is diff, even though the estimator itself remains the same. S
## For best estimation of the Stderr, set ieff = (a-1), the last available 
## 	age group to obtain the small stderr.
##
##	fam = "pois" : loglinear regression; 
##	fam = "qlik" : quasi-likelihood adjusting for over-(under-) dispersion. 
	fit <- apcrdg(r = r, nrisk = nrisk, c1 = 1, c2 = 1, cc = 1, apc = 1, 
		intercept = T, rr = F, lambda = lambda, gamma = 2, center = 
		center, Plot = Plot, tnull = 0)
	sol <- vector("list")
	if(stderr == F)
		sol <- fit
	else {
		if(princomp == 1) {
			a <- nrow(r)
			p <- ncol(r)
			x <- apcmat(a, p)[, -1]
			v <- eigen(t(x) %*% x)$vector
			xx <- x %*% v[,  - (ncol(v))]
			if(offst == F) {
				if(fam == "pois")
				  fitlog <- glm(as.vector(t(r)) ~ xx, family = 
				    poisson(link = log), control = glm.control(
				    epsilon = 1e-010, maxit = niter))
				if(fam == "qlik")
				  fitlog <- glm(as.vector(t(r)) ~ xx, family = 
				    quasi(link = log, var = "mu"), control = 
				    glm.control(epsilon = 1e-010, maxit = niter))
			}
			else {
				r <- r * nrisk
				if(fam == "pois")
				  fitlog <- glm(as.vector(t(r)) ~ xx + offset(log(
				    as.vector(t(nrisk)))), family = poisson(link
				     = log), control = glm.control(epsilon = 
				    1e-010, maxit = niter))
				if(fam == "qlik")
				  fitlog <- glm(as.vector(t(r)) ~ xx + offset(log(
				    as.vector(t(nrisk)))), family = quasi(link = 
				    log, var = "mu"), control = glm.control(
				    epsilon = 1e-010, maxit = niter))
			}
			b.int <- fitlog$coef[1]
			b <- as.vector(v %*% as.matrix(c(fitlog$coef[-1], 0)))
			aa <- b[1:(a - 1)]
			pp <- b[a:(a + p - 2)]
			coh <- b[ - c(1:(a + p - 2))]
			bvar <- diag(summary.glm(fitlog)$coef[-1,  ][, 2]) %*% 
				summary.glm(fitlog)$corr[-1,  ][, -1] %*% diag(
				summary.glm(fitlog)$coef[-1,  ][, 2])
			sol$dev <- fitlog$dev
			sol$df <- fitlog$df
			sol$int <- rbind(b.int, summary.glm(fitlog)$coef[1, 2])
			bvarall <- v %*% cbind(rbind(bvar, 0), 0) %*% t(v)
			a.se <- sqrt(matrix(-1, 1, a - 1) %*% bvarall[1:(a - 1), 
				1:(a - 1)] %*% matrix(-1, a - 1, 1))
			p.se <- sqrt(matrix(-1, 1, p - 1) %*% bvarall[a:(a + p - 
				2), a:(a + p - 2)] %*% matrix(-1, p - 1, 1))
			c.se <- sqrt(matrix(-1, 1, a + p - 2) %*% bvarall[(a + p - 
				1):(2 * (a + p - 2)), (a + p - 1):(2 * (a + p - 2
				))] %*% matrix(-1, a + p - 2, 1))
			se.all <- sqrt(diag(bvarall))
			a.coef <- rbind(c(aa,  - sum(aa)), c(se.all[1:(a - 1)], 
				a.se))
			p.coef <- rbind(c(pp,  - sum(pp)), c(se.all[a:(a + p - 2)
				], p.se))
			c.coef <- rbind(c(coh,  - sum(coh)), c(se.all[-1: - (a + 
				p - 2)], c.se))
			sol$aa <- round(a.coef, 4)
			sol$pp <- round(p.coef, 4)
			sol$cc <- round(c.coef, 4)
			if(fam == "qlik")
				sol$dispersion <- summary.glm(fitlog)$disp
			if(CIplot == T && Plot == F) {
				par(mfrow = c(2, 2))
				if(pylim[2] == pylim[1])
				  plot((c(1:a) - 1) * gapyear + amin, a.coef[1,  
				    ], xlab = "age", ylab = "age effect", type = 
				    "b")
				else plot((c(1:a) - 1) * gapyear + amin, a.coef[1,
				    ], xlab = "age", ylab = "age effect", type = 
				    "b", ylim = pylim)
				lines((c(1:a) - 1) * gapyear + amin, a.coef[1,  ] +
				  1.96 * a.coef[2,  ], lty = 2)
				lines((c(1:a) - 1) * gapyear + amin, a.coef[1,  ] -
				  1.96 * a.coef[2,  ], lty = 2)
				title(main = "Age trend")	
	#		title (main="Age trend, 95% CI")
				if(pylim[2] == pylim[1])
				  plot((c(1:p) - 1) * gapyear + pmin, p.coef[1,  
				    ], xlab = "period", ylab = "period effect", 
				    type = "b")
				else plot((c(1:p) - 1) * gapyear + pmin, p.coef[1,
				    ], xlab = "period", ylab = "period effect", 
				    type = "b", ylim = pylim)
				lines((c(1:p) - 1) * gapyear + pmin, p.coef[1,  ] +
				  1.96 * p.coef[2,  ], lty = 2)
				lines((c(1:p) - 1) * gapyear + pmin, p.coef[1,  ] -
				  1.96 * p.coef[2,  ], lty = 2)
				title(main = "Period trend")	
	#		title (main="Period trend, 95% CI")
				par(mfrow = c(2, 1), new = F, fig = c(0, 1, 0, 
				  0.5))
				if(pylim[2] == pylim[1])
				  plot((c(1:(a + p - 1)) - 1) * gapyear + cmin, 
				    c.coef[1,  ], xlab = "cohort", ylab = 
				    "cohort effect", type = "b")
				else plot((c(1:(a + p - 1)) - 1) * gapyear + cmin,
				    c.coef[1,  ], xlab = "cohort", ylab = 
				    "cohort effect", type = "b", ylim = pylim)
				lines((c(1:(a + p - 1)) - 1) * gapyear + cmin, 
				  c.coef[1,  ] + 1.96 * c.coef[2,  ], lty = 2)
				lines((c(1:(a + p - 1)) - 1) * gapyear + cmin, 
				  c.coef[1,  ] - 1.96 * c.coef[2,  ], lty = 2)
				title(main = "Cohort trend")	
	#		title (main="Cohort trend, 95% CI")
				par(mfrow = c(1, 1))
			}
		}
		else {
			if(fit$aa[1] != 0) {
				c1 <- 1
				c2 <- ieff
				cc <- fit$aa[ieff]/fit$aa[1]
				apc <- 1
			}
			else if(fit$pp[1] != 0) {
				c1 <- 1
				c2 <- ieff
				cc <- fit$pp[ieff]/fit$pp[1]
				apc <- 2
			}
			else if(fit$coh[1] != 0) {
				c1 <- 1
				c2 <- ieff
				cc <- fit$coh[ieff]/fit$coh[1]
				apc <- 3
			}
			else {
				stop(
				  "first effect = 0 for all age, period and cohort"
				  )
			}
			fit <- apcrdg(r = r, nrisk = nrisk, c1 = c1, c2 = c2, cc
				 = cc, apc = apc, intercept = T, rr = F, lambda
				 = 0, gamma = 2, center = center, Plot = F, tnull
				 = 0)
			sol <- fit
		}
	}
	sol
}

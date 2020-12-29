constr <-
function(x, nr, nc, c1, c2, cc, apc = 1)
{
# APC = 1, age; = 2, period; =3, cohort.
# c1, c2: the column numbers of constrait.
# cc : the coef for the constraint, Beta2 = cc * Beta1.
#nc <- ncol (x)
#nr <- nrow (x)
	if(apc == 1) {
		cc1 <- min(c1, c2) + 1
		cc2 <- max(c1, c2) + 1
	}
	if(apc == 2) {
		cc1 <- min(c1, c2) + nr
		cc2 <- max(c1, c2) + nr
	}
	if(apc == 3) {
		cc1 <- min(c1, c2) + nr + nc - 1
		cc2 <- max(c1, c2) + nr + nc - 1
	}
	c1 <- cc1
	c2 <- cc2
	if(c2 - c1 == 1) {
		xc <- cbind(x[, 1:(c1 - 1)], x[, c1] + cc * x[, c2], x[, c(-1: - 
			c2)])
	}
	if(c2 - c1 > 1) {
		xc <- cbind(x[, 1:(c1 - 1)], x[, c1] + cc * x[, c2], x[, (c1 + 1):
			(c2 - 1)], x[, -1: - c2])
	}
	if(c1 == c2) {
		xc <- x
	}
#endif
	xc
}

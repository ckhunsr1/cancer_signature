apcmat <-
function(a, p)
{
## matrix for APC model, see Fu (1998)
	alpha <- rbind(diag(rep(1, a - 1)), rep(-1, a - 1))
	beta <- rbind(diag(rep(1, p - 1)), rep(-1, p - 1))
	gamma <- rbind(diag(rep(1, a + p - 2)), rep(-1, a + p - 2))
	x <- rep(0, 2 * (a + p) - 3)
	for(i in 1:a) {
		for(j in 1:p) {
			rr <- c(1, alpha[i,  ], beta[j,  ], gamma[a - i + j,  ])
			x <- rbind(x, rr)
		}
	}
	x <- x[-1,  ]
	x
}

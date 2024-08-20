library(microbenchmark)

N <- c(100,1000,10000,100000)
res <- matrix(NA,length(N),3)
names <- c("crossprod","sum","%*%")


for (j in 1:length(N)) {
  n <- N[j]
  v1 <- rnorm(n)
  v1_t = t(v1)
  v2 <- rnorm(n)
  x <- microbenchmark(crossprod = {x <- numeric(1); x <- x + as.numeric(crossprod(v1,v2))},
                      sum = {x <- numeric(1); x <- x + sum(v1*v2)},
                      #"%*%" = {x <- numeric(1); x <- x + v1_t %*% v2},
                      unit = "relative",
                      times = 1000
  )
  for (i in 1:3) {
    res[j,i] <- median(x$time[x$expr == names[i]])
  }
  res[j, ] <- res[j, ] / min(res[j, ])
}

colnames(res) <- names
rownames(res) <- N
t(res)
#           100 1000    10000    1e+05
# crossprod 1.0 1.00 1.000000 1.000000
# sum       1.4 1.00 1.867925 1.946708
# %*%       1.4 1.84 1.971698 1.993368

# xtable::xtable(t(res), digits = 4)

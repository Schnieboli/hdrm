## hdrm1 with matrix
N <- 20
d <- 30
# example data with N = 20 subjects in columns and d = 30 factor levels in rows
M <- matrix(rnorm(d*N), d, N)

# test M for H_0 = mu_1 = mu_2 = mu_3 = 0 time profile
H <- matrix(0,d,d)
H[c(1, d+2, 2*d+3)] <- 1
hdrm1(data = M, hypothesis = H)

## hdrm1 with data.frame in longtable format
head(EEG)
# for 'formula' we are interested in value = value, subject = subject and factor = dimension
formula <- value ~ subject + dimension

# test EEG for equal time profile
hdrm1(data = EEG, formula = formula, hypothesis = "equal")

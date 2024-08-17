## die Idee war, dass zunächst ein hdrm-Object erstellt wir, mit dem dann beide tests durchgeführt werden können
## das sollte das Problem lösen, dass es ein bisschen weird ist mit der Eingabe durch das Formelobjekt
## -> es wäre so aber immer noch weird, da man jan icht immer gruppen braucht... das macht die eingabe auch komisch

create_hdrm <- function(data, group, subject, dimension, ...){
  UseMethod("create_hdrm")
}

create_hdrm_data_frame <- function(data, value, subject, group = rep(1, N), dimension){
  stopifnot(is.data.frame(data))

  N

  f <- value ~ subject + group + dimension
  df <- stats::model.frame(formula = formula, data = data)


}

create_hdrm_matrix <- function(data, group, subject, dimension, ...){
  N <- ncol(data)
  d <- nrow(data)
  a <- length(unique(group))



}

create_hdrm_list <- function(data, group, subject, dimension,...){
  stopifnot(all(sapply(data, is.matrix)))
  d <- sapply(data, nrow)
  stopifnot(length(unique(d)) == 1)

  a <- length(data)
  n <- sapply(L, ncol)
  N <- sum(n)
  d <- d[1]

  group <- numeric(0)
  for (i in 1:a) {
    group <- c(group, rep(i, n[i]))
  }

  group <- as.factor(group)
  if(!is.null(names(data)) & all(is.na(names(data)))) levels(group <- names(data))


  mat <- matrix(0.0, nrow = d, ncol = 0)
  for (i in 1:a) {
    mat <- cbind(mat, L[[i]])
  }

  create_hdrm.matrix(data = X, group = group)

}

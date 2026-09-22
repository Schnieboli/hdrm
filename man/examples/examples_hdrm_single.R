# Long-format data ---------------------------------------------------------

data("EEG")
names(EEG) # already contains columns value, dimension and subject

## select only one diagnostic group for one-group analysis
EEG_single <- EEG[EEG$group == "SCC+", ]

## test whether the time profile is flat
hdrm_single(
  data = EEG_single,
  hypothesis = "flat"
)

## define hypothesis = "flat" via equivalent projection matrix
d <- nlevels(EEG_single$dimension)
flat_projection <- diag(d) -
  matrix(1 / d, nrow = d, ncol = d)

## test whether the time profile is flat via custom hypothesis matrix
hdrm_single(
  data = EEG_single,
  hypothesis = flat_projection
)


# Wide matrix data ---------------------------------------------------------

data("birthrates")

## transform 'birthrates' to matrix and transpose
birthrates_matrix <- t(as.matrix(birthrates))

## test whether the time profile is flat
hdrm_single(
  data = birthrates_matrix,
  hypothesis = "flat"
)

## define hypothesis = "flat" via equivalent projection matrix
d <- ncol(birthrates_matrix)
flat_projection <- diag(d) -
  matrix(
    1 / d,
    nrow = d,
    ncol = d
  )

## test whether the time profile is flat via custom hypothesis matrix
hdrm_single(
  data = birthrates_matrix,
  hypothesis = flat_projection
)

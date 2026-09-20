# Long-format data ---------------------------------------------------------

data("EEG")

# Select one diagnostic group for a one-group analysis.
EEG_single <- droplevels(EEG[EEG$group == "SCC+", ])

hdrm_single(
  data = EEG_single,
  hypothesis = "flat"
)

# A custom projection matrix equivalent to hypothesis = "flat".
d <- nlevels(EEG_single$dimension)
flat_projection <- diag(d) -
  matrix(1 / d, nrow = d, ncol = d)

hdrm_single(
  data = EEG_single,
  hypothesis = flat_projection
)


# Wide matrix data ---------------------------------------------------------

data("birthrates")

birthrates_matrix <- t(as.matrix(birthrates))

hdrm_single(
  data = birthrates_matrix,
  hypothesis = "flat"
)

d <- ncol(birthrates_matrix)
flat_projection <- diag(d) -
  matrix(
    1 / d,
    nrow = d,
    ncol = d
  )

hdrm_single(
  data = birthrates_matrix,
  hypothesis = flat_projection
)

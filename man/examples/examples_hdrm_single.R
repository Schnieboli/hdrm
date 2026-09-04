# Long-format data ---------------------------------------------------------

data("EEG")

# Select one diagnostic group for a one-group analysis.
eeg_single <- droplevels(
  EEG[EEG$group == levels(EEG$group)[1L], ]
)

eeg_single <- eeg_single[
  order(eeg_single$subject, eeg_single$dimension),
]

hdrm_single(
  data = eeg_single$value,
  hypothesis = "flat",
  subject = eeg_single$subject
)

# A custom projection matrix equivalent to hypothesis = "flat".
d <- nlevels(eeg_single$dimension)
flat_projection <- diag(d) -
  matrix(1 / d, nrow = d, ncol = d)

hdrm_single(
  data = eeg_single$value,
  hypothesis = flat_projection,
  subject = eeg_single$subject
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

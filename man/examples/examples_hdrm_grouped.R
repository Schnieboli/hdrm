# Long-format data ---------------------------------------------------------

data("EEG")

# Subject labels need only be unique within groups. The original labels can
# therefore be used directly, even when the same label occurs in another group.
eeg_grouped <- EEG[
  order(EEG$group, EEG$subject, EEG$dimension),
]

# Heterogeneous covariance matrices; the available exact trace estimators are
# used, while the third-trace quantity is still estimated by subsampling.
hdrm_grouped(
  data = eeg_grouped$value,
  hypothesis = "whole",
  group = eeg_grouped$group,
  subject = eeg_grouped$subject,
  cov.equal = FALSE,
  subsampling = FALSE,
  B = "100*N",
  seed = 3141
)

# Replace the remaining available trace estimators by their subsampling
# versions. Here B is used by each subsampling-estimator invocation.
hdrm_grouped(
  data = eeg_grouped$value,
  hypothesis = "sub",
  group = eeg_grouped$group,
  subject = eeg_grouped$subject,
  cov.equal = FALSE,
  subsampling = TRUE,
  B = 10000,
  seed = 3141
)

# A custom projection-matrix hypothesis equivalent to hypothesis = "sub".
custom_eeg_hypothesis <- list(
  TW = matrix(1 / 4, nrow = 4, ncol = 4),
  TS = diag(40) - matrix(1 / 40, nrow = 40, ncol = 40)
)

hdrm_grouped(
  data = eeg_grouped$value,
  hypothesis = custom_eeg_hypothesis,
  group = eeg_grouped$group,
  subject = eeg_grouped$subject,
  B = "100*N",
  seed = 3141
)


# Wide matrix data ---------------------------------------------------------

data("birthrates")

birthrates_matrix <- t(as.matrix(birthrates))

group <- factor(
  c(1, 1, 2, 2, 1, 1, 1, 2, 1, 1, 1, 1, 2, 2, 1, 2),
  labels = c("west", "east")
)

hdrm_grouped(
  data = birthrates_matrix,
  hypothesis = "interaction",
  group = group,
  cov.equal = FALSE,
  subsampling = FALSE,
  B = "100*N",
  seed = 3141
)

# Under equal covariance matrices, B is the total pooled subsampling budget
# across groups and the value of subsampling is ignored.
hdrm_grouped(
  data = birthrates_matrix,
  hypothesis = "whole",
  group = group,
  cov.equal = TRUE,
  B = "100*N",
  seed = 3141
)

d <- ncol(birthrates_matrix)
custom_birthrates_hypothesis <- list(
  TW = matrix(1 / 2, nrow = 2, ncol = 2),
  TS = diag(d) - matrix(1 / d, nrow = d, ncol = d)
)

hdrm_grouped(
  data = birthrates_matrix,
  hypothesis = custom_birthrates_hypothesis,
  group = group,
  B = "100*N",
  seed = 3141
)

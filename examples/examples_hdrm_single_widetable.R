data(birthrates)

# call test
hdrm_single_widetable(data = birthrates,
                       # test if time profiles are flat
                       hypothesis = "flat"
)


# hypothesis as list
hypothesis_matrix <- diag(34) - matrix(1/34, 34, 34)

hdrm_single_widetable(data = birthrates,
                      # equivalent to hypothesis = "flat"
                      hypothesis = "flat"
)

rm(hypothesis_matrix)

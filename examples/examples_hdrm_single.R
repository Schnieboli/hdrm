## Load data set EEG (data frame)
data(EEG)
# ?EEG


## call test
hdrm_single(data = EEG$value,
            ## test if time profiles are flat
            hypothesis = "flat",
            ## if data is a vector, subject must be specified
            subject = EEG$subject
)


## hypothesis as list
hypothesis_matrix <- diag(40) - matrix(1/40, 40, 40)

hdrm_single(data = EEG$value,
            ## test if time profiles are flat
            hypothesis = hypothesis_matrix,
            subject = EEG$subject
)

## Load data set birthrates (matrix)
data(birthrates)
# ?birthrates

## call test

hdrm_single(data = birthrates,
            ## test if time profiles are flat
            hypothesis = "flat"
            ## if data is a matrix, subject does not need to be specified
)


## hypothesis as list
hypothesis_matrix <- diag(34) - matrix(1/34, 34, 34)

hdrm_single(data = birthrates,
            ## equivalent to hypothesis = "flat"
            hypothesis = hypothesis_matrix
)

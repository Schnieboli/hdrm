## load data set EEG (data frame)
data(EEG)

# call test
hdrm_grouped(data = EEG$value,
             # test for no group effect
             hypothesis = "whole",
             group = EEG$group,
             ## if data is a vector, subject has to be specified
             subject = EEG$subject,
             subsampling = FALSE,
             B = "100*N"
)


# test using all subsampling version of the trace estimators
hdrm_grouped(data = EEG$value,
             # test for no time effect
             hypothesis = "sub",
             group = EEG$group,
             subject = EEG$subject,
             subsampling = TRUE,
             # B can also be given as a number
             B = 10000
)


# hypothesis as list
hypothesis_list <- list(TW = matrix(1/4, 4, 4),
                        TS = diag(40) - matrix(1/40, 40, 40)
) # equivalent to hypothesis = "sub"

# call test
hdrm_grouped(data = EEG$value,
             hypothesis = hypothesis_list,
             group = EEG$group,
             subject = EEG$subject,
             subsampling = FALSE,
             B = "100*N"
)
rm(hypothesis_list)


## load data set birthrates (matrix)
data(birthrates)
# ?birthrates

## divide states into (former) east and west with Berlin as east
group <- factor(c(1,1,2,2,1,1,1,2,1,1,1,1,2,2,1,2), labels = c("west", "east"))

## call test for interaction effect
hdrm_grouped(data = birthrates,
             ## test for interaction effect between group and dimension
             hypothesis = "interaction",
             group = group,
             ## if data is a matrix, subject does not have to be specified
             B = "100*N"
)


## test using all subsampling version of the trace estimators
hdrm_grouped(data = birthrates,
             ## test for group effect between group and dimension
             hypothesis = "whole",
             group = group,
             ## use the subsampling version of all trace estimators
             subsampling = TRUE,
             ## B can also be given as a number
             B = 10000
)


## hypothesis as list, equivalent to hypothesis = "sub"
hypothesis_list <- list(TW = matrix(1/2, 2, 2),
                        TS = diag(34) - matrix(1/34, 34, 34)
)

## call test with hypothesis as list
hdrm_grouped(birthrates,
             ## equivalent to hypothesis = "sub"
             hypothesis = hypothesis_list,
             group = group,
             subsampling = TRUE,
             B = "100*N"
)

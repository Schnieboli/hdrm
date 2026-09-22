data(birthrates)

# divide states into (former) east and west
group <- factor(c(1,1,2,2,1,1,1,2,1,1,1,1,2,2,1,2), labels = c("west","east"))
# call test
hdrm_grouped_widetable(data = birthrates,
                       # test for interaction effect between group and dimension
                       hypothesis = "interaction",
                       group = group,
                       bootstrap = FALSE,
                       B = "100*N"
)


# test using all subsampling version of the trace estimators
hdrm_grouped_widetable(data = birthrates,
                       # test for group effect between group and dimension
                       hypothesis = "whole",
                       group = group,
                       subsampling = TRUE,
                       # B can also be given as a number
                       B = 10000
)


# hypothesis as list
hypothesis_list <- list(TW = matrix(1/2, 2, 2),
                        TS = diag(34) - matrix(1/34, 34, 34)
) # equivalent to hypothesis = "sub"

# test call
hdrm_grouped_widetable(birthrates,
                       hypothesis = hypothesis_list,
                       group = group,
                       bootstrap = FALSE,
                       B = "100*N"
)
rm(hypothesis_list)

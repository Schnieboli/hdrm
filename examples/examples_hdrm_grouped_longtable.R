data(EEG)
names(EEG)
# name (index) of columns of interest are:
# group = group (1)
# value = value (2)
# subject = subject (4)
# dimension = dimension (7)

# call test
hdrm_grouped_longtable(data = EEG,
                       # test for no group effect
                       hypothesis = "whole",
                       group = "group",
                       value = "value",
                       # columns can also be specified by their indices
                       subject = 4,
                       dimension = "dimension",
                       subsampling = FALSE,
                       B = "100*N"
)


# test using all subsampling version of the trace estimators
hdrm_grouped_longtable(data = EEG,
                       # test for no time effect
                       hypothesis = "sub",
                       group = "group",
                       value = "value",
                       subject = "subject",
                       dimension = "dimension",
                       subsampling = TRUE,
                       # B can also be given as a number
                       B = 10000
)


# hypothesis as list
hypothesis_list <- list(TW = matrix(1/4, 4, 4),
                        TS = diag(40) - matrix(1/40, 40, 40)
) # equivalent to hypothesis = "sub"

# call test
hdrm_grouped_longtable(data = EEG,
                       hypothesis = hypothesis_list,
                       group = "group",
                       value = "value",
                       subject = "subject",
                       dimension = "dimension",
                       subsampling = FALSE,
                       B = "100*N"
)
rm(hypothesis_list)

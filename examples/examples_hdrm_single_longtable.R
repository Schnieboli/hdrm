data(EEG)
names(EEG)
# name (index) of columns of interest are:
# value = value (2)
# subject = subject (4)
# dimension = dimension (7)

# call test
hdrm_single_longtable(data = EEG,
                      # test if time profiles are flat
                      hypothesis = "flat",
                      value = "value",
                      # columns can also be specified by their indices
                      subject = 4,
                      dimension = "dimension"
)


# hypothesis as list
hypothesis_matrix <- diag(40) - matrix(1/40, 40, 40)

hdrm_single_longtable(data = EEG,
                      # test if time profiles are flat
                      hypothesis = hypothesis_matrix,
                      value = "value",
                      subject = "subject",
                      dimension = "dimension"
)

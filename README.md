# Allocation_CPUE_Covariate_Sharing

Accompanies a paper by Hall, Johnson, Smith (in review).

This simulates fish abundance combining two species (A and B) which have
different habitat preferences. These preferences are used to split the total
abundance into the two species. One species is then standardised, with the
covariates (X1 and X2) used in both models.

The goal is to see if covariate reuse biases the standardised index, which
here has a declining trend. We have two switches: one for whether X2 (~depth)
influences catchability as well as density, and one for whether fishing
effort drifts into deeper water over time.

We see that sharing covariates is generally safe in a well specified model
whereas common situations (like changes in the spatial distribution of
effort) can bias the index if these covariates are not reused.
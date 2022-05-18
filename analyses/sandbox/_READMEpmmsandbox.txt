Started 15 May 2022
By Lizzie

Notes on trying to get a family-level (with species partially pooled below) running with our phylo model
This would take a while to build up, so we should start with test data.


Mike Betancourt provided some very QUICK and incomplete code: oneslopeinterceptcholforsync_withMike.stan


Before that I did some sorting and flailing ... I edited:

* simulate_fit_oneslopeintercept_cholesky_dl.R 

It returns test data values (though I am not convinced it offers a speed up).

I then copied this and started trying to generate test data for the nested model in:

* simulate_fit_oneslopeintercept_cholesky_fam.R 

.... But I am not sure I am doing it right *and* I really should start with code that does not generate intercepts.

I then flailed in the Stan code, with:

* oneslopecholforsyncwfam.stan


So, what are reasonable options for next steps?
1) Look at my general notes on nested models. Organize and improve them if I need to. DONE, see updated repo (and renamed): nestedmodels 
2) Fix the R code so it does the nesting we want and leaves the intercept at 0 or such.
3) Try again on the Stan code. 


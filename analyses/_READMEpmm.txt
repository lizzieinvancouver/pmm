23 August 2020
Updated 15 May 2022
Lizzie finally making some notes on what's in the folder 

Looking for Lemoine PGLS? 
>> git/teaching/gelmanhill/BayesPGLS/pgls.R
>> git/teaching/gelmanhill/BayesPGLS/pgls_lemoine.stan

phlyo_opsree_compact1.R -- Originally, December work using ubermini_2 on OSPREE. Current version fits a one trait (i.e., oneslope) model to the OSPREE data. Model is stan/uber_oneslopeintercept.stan

phlyo_opsree_compact2.R -- Fits a two trait (i.e., twoslope) model to the OSPREE data. Model is stan/uber_twoslopeintercept.stan

phlyo_opsree_compact3.R -- Fits a three trait (i.e., threeslope) model to the OSPREE data. Model is stan/uber_threeslopeintercept.stan

simsmore folder -- work by Geoff L. in 2021, based off ubermini_2.R and related Stan code to build up simulations with multiple slopes and an intercept (simulate_fit_oneslope is a slope-only model, so assuming an intercept of 0)

burstdeltamodels folder -- work started in November 2022 to look at early/late burst models of evolution ... incl. what happens when use sigma/lambda models but set evolution as delta 

*****
ubermini_2.R -- Will's new attempt in late August 2020 with updates from Geoff to set initial conditions ** This one works! **

Geoff Legault added a ton of code in later 2021; much of it I believe duplicating what we already have, but de novo written by him with small tweaks. *We should check and clean this repo someday, as we now have too much code* but for now, here's some notes by Lizzie:


uber_oneslope.stan -- very similar (I think -- need to check) to ubermini_2.stan but with some tweaks (esp. priors) by Geoff

uber_oneslopeintercept.stan -- very similar (I think -- need to check) to ubermini_2.stan but with phylo on slope and intercept (by Geoff)

uber_oneslopeintercept_cholesky.stan -- added cholesky uber_oneslopeintercept.stan (by Geoff)


<><><><><><><><><>
Work by Dinara Mamatova (summer 2022) to compare speed of different formulas
See https://github.com/lizzieinvancouver/pmm/issues/19
speedtests folder
speedtests_viz folder (visualization)

*****

<><><><><><><><><>
wanderings folder (created 29 June 2021, same time Lizzie deleted a bunch of stuff)

max_sims.R from Max Farrell

pmm-pgls-simulations_v2.R -- original code from Simon Joly

pmm-pgls-simulations_v2IMC.R -- edits to Joly code for our purposes

run_MCMCbrms_pglmm.R -- edits to Joly code for our purposes, compares MCMCglmm and brms

uberlessmini.R -- similar to ubermini.R but puts phylo on intercept.

ubermini.R -- my original code to test Will's model, simple test data with intercept of 0, and just slope (with phylo on slope) on linear model

ubertoy_phyr.R -- trying to run model using R package by Ives (phyr)

wanderings/evenolder folder

run_MCMCbrms_pglmm_nointra.R  -- my edits to Joly code to remove intraspecific covariance matrix ... I think superseded by pmm-pgls-simulations_v2IMC.R  (but should check)

notes-from-a-small-island.R -- From Will Pearse, something like PGLS I think ... see pmm_emails.txt

Pearse-code-4brmsdata.R -- Nacho's code trying to run Will's early (spring 2020) models

supplement.R -- code from Will's paper, but has the early method (lambda + alpha) that we think is not working


<><><><><><><><><>
Extra notes ...

PMM in OSPREE:
* Phylo_ospree_reanalyses.R has all the work...
<><><><><><><><><>

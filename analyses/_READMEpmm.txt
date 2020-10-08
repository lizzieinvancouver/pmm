23 August 2020
Lizzie finally making some notes on what's in the folder 

Looking for Lemoine PGLS? 
>> git/teaching/gelmanhill/BayesPGLS/pgls.R
>> git/teaching/gelmanhill/BayesPGLS/pgls_lemoine.stan


max_sims.R from Max Farrell

pmm-pgls-simulations_v2.R -- original code from Simon Joly

pmm-pgls-simulations_v2IMC.R -- edits to Joly code for our purposes

run_MCMCbrms_pglmm.R -- edits to Joly code for our purposes, compares MCMCglmm and brms

uberlessmini.R -- similar to ubermini.R but puts phylo on intercept.

ubermini.R -- my original code to test Will's model, simple test data with intercept of 0, and just slope (with phylo on slope) on linear model

ubermini_2.R -- Will's new attempt in late August 2020 ** This one might work! Has different multi normal formula than the other code!

uberpgls.R -- my attempts to write PGLS code, not successful yet. 

Below are similar to above, but tidied up for Tony Ives: 
ubertoy_nogeiger.R -- this one builds the data exactly as I believe I have coded the model
ubertoy.R -- this one uses geiger to adjust the phylogeny by lambda 



<><><><><><><><><>
notcurrent folder

run_MCMCbrms_pglmm_nointra.R  -- my edits to Joly code to remove intraspecific covariance matrix ... I think superseded by pmm-pgls-simulations_v2IMC.R  (but should check)

notes-from-a-small-island.R -- From Will Pearse, something like PGLS I think ... see pmm_emails.txt

Pearse-code-4brmsdata.R -- Nacho's code trying to run Will's early (spring 2020) models

supplement.R -- code from Will's paper, but has the early method (lambda + alpha) that we think is not working


<><><><><><><><><>
Next steps, various steps etc.
for Lizzie (22 September 2020)


PMM test data:
* Put up an example of ubertoy using Will's method on Stan Discourse
* Work on the code in phyr (get running on OSPREE maybe?)
* organize the damn files here
* follow up with Will, Simone, Tony

PMM in OSPREE:
* Phylo_ospree_reanalyses.R has all the work...
* As of 22 Sep I cannot get the version Will sent in late August (ubermini_2.R) to run (nointer_2levelphyall_2.stan in OSPREE), even with just forcing so I need to work on this more. 
<><><><><><><><><>
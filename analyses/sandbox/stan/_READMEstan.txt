`Started 19 May 2022
By Lizzie

Files checking on pulling out the grandmeans:
	oneslopeinterceptcholforsync_grandmean.stan (pulls out intercept)
	oneslopeinterceptcholforsync_grandmeans.stan (pulls out intercept and slope)

Files with intercept only model (good place to start for working on new additions or changes to model):
	oneinterceptcholforsync.stan -- pulled out grand mean also

Files started 18 May 2022 by Lizzie to try to include phylogeny at higher level (we imagined family) and fit species underneath that without the GP structure. 
	oneinterceptcholforsyncfam.stan -- this is possibly the correct code (Mike B helped a little as he could_, but we need to track down the reason for the divergences 
        oneinterceptcholforsyncfamncp.stan -- quick attempt by Lizzie to non-center the one part of the code she knew how .... But have not even run it. Needs work. 

Not functional code where Mike jotted down some example code of: (1) making our phylogeny model include a nested level and (2) showing how to vectorize code. I added some notes but otherwise have not touched this. 
	oneslopeinterceptcholforsync_withMike.stan

/* From supplement.R */

data{
    int Ntotal;
    int Nspp;
    int presence[Ntotal];
    real env[Ntotal];
    int spp[Ntotal];
    //    
    matrix[Nspp,Nspp]Vphy;     // Give the phylogeny as data
}
parameters{
    vector[Nspp] spp_intercepts;
    vector[Nspp] spp_slopes;
    // Coefficients for the phylogenetically-derived variance in model terms
    real<lower=0,upper=100> lam_intercepts;    // (with priors specified too) 
    real<lower=0,upper=100> lam_slopes;        // 
    // Coefficients for the NON-phylogenetically-derived variance in model terms
    real<lower=0,upper=100> null_intercepts;       // (also with priors)
    real<lower=0,upper=100> null_slopes;           //
}
transformed parameters{
    vector[Ntotal] predictions;
    for (i in 1:Ntotal)
        predictions[i] = spp_intercepts[spp[i]] + spp_slopes[spp[i]]*env[i];
}
model{
    // Now we draw our species coefficients to measure the importance of phylogeny
    spp_intercepts ~ multi_normal(rep_vector(0,Nspp), diag_matrix(rep_vector(null_intercepts,Nspp)) + lam_intercepts*Vphy); 
    spp_slopes ~ multi_normal(rep_vector(0,Nspp), diag_matrix(rep_vector(null_slopes,Nspp)) + lam_slopes*Vphy);
    //
    presence ~ bernoulli_logit(predictions);
}

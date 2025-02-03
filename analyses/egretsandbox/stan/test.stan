data {
  int n; 
  int admit[n];
  int applications[n]; //number of applications for each outcome
  vector[n] gender;
}

parameters {
  real alpha;
  vector[1] beta;
}

model {
  //linear model
  vector[n] p;
  p = alpha + beta[1] * gender;
  
  //likelihood and link function
  admit ~ binomial_logit(applications, p);
}
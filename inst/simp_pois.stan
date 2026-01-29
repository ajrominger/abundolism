// input data:
// samp size: `N`
// response: counts in `y`
// explanatory: reals in `x`
data {
  int<lower=0> N;
  array[N] int<lower=0> y;
  vector[N] x;
}

// parameters:
// a1 is intercept of log mean
// a2 is slope of x for log mean
parameters {
  real a_1;
  real a_2;
}


transformed parameters {
  vector[N] log_mu = a_1 + a_2 * x;
  vector[N] mu = exp(log_mu);
}


model {
  // priors
  a_1 ~ normal(0, 50);
  a_2 ~ normal(0, 50);

  // likelihood
  y ~ poisson(mu);
}


data {
  int<lower=0> N;
  array[N] int<lower=0> y;
  vector[N] x;
}

parameters {
  real a_0;
  real a_1;
  real a_2;
  real a_3;
}

transformed parameters {
  vector[N] logit_p = a_0 + a_1 * x + a_2 * square(x);
  vector[N] p = inv_logit(logit_p);
  vector[N] lambda = p * a_3;
}

model {
  // Priors
  a_0 ~ normal(0, 100);
  a_1 ~ normal(0, 100);
  a_2 ~ normal(0, 100);
  a_3 ~ gamma(10, 0.1);

  // Likelihood
  y ~ poisson(lambda);
}

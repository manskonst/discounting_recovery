data {
  int<lower=0> Ntotal;
  int y[Ntotal];
  vector[Ntotal] r_a;
  vector[Ntotal] r_b;
  vector[Ntotal] d_a;
  vector[Ntotal] d_b;
}  
 
parameters {
  real<lower=0> k;
  real<lower=0,upper=0.5> alpha;
}

transformed parameters {
  vector[Ntotal] theta;
  { //start local variables
  real v_diff;
  for(i in 1:Ntotal){
    v_diff = r_a[i]/(1+k*d_a[i]) - r_b[i]/(1+k*d_b[i]);
    if(v_diff>0){
      theta[i] = (1-alpha);
    } else {
      if(v_diff<0){
        theta[i] = alpha;
      } else {
        if(v_diff==0){
          theta[i] = 0.5;
        }
      }
    }
  }
}
}
  
    
model {
  k~normal(0,5);
  y~bernoulli(theta);
}
  
//generated quantities {
//  vector[Ntotal] log_lik;
//  for (n in 1:Ntotal){
//    log_lik[n] = bernoulli_logit_lpmf(y[n] | p_a_logit[n]);
//  }
//}

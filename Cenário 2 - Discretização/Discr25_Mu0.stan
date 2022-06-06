
  data {
    int<lower=1> k;
    int N[k];
    vector[k] z;
    vector[k] r;
    matrix<lower=0>[k, k] Mdist;
  }
      
  transformed data {
    vector[k] zero;
    for (i in 1:k)
      zero[i] = 0;
  }

  parameters {
    real beta;
    real<lower=0> gamma;
    real<lower=0> tau;
    vector[k] w;
  }
  
  model {
    vector[k] lambda_min;
    vector[k] lambda_mai;
    matrix[k, k] sigmaMat;
    beta ~ normal(0, 100);   
    tau ~ gamma(1, 0.1); 
    gamma ~ gamma(1, 0.06708);
    sigmaMat = (1/tau) * exp(-Mdist/gamma);
    w ~ multi_normal(zero, sigmaMat);
    for (j in 1:k) {
      lambda_min[j] = exp(z[j]*beta + w[j]);
      lambda_mai[j] = r[j]*lambda_min[j];
      N[j] ~ poisson(lambda_mai[j]);
    }
  }
  
  generated quantities {
    real<lower=0> sigma2 = 1/tau;
  }
  
                     

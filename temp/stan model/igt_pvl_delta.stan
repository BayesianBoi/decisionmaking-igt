data {
    int<lower=1> N;
    int<lower=1> T;
    int<lower=1, upper=T> Tsubj[N];
    real rewlos[N, T];
    int ydata[N, T];
}

transformed data {
    vector[4] initV;
    initV  = rep_vector(0.0,4);
}

parameters {
    real mu_A_pr;
    real mu_alpha_pr;
    real mu_cons_pr;
    real mu_lambda_pr;
    
    real<lower=0> sd_A;
    real<lower=0> sd_alpha;
    real<lower=0> sd_cons;
    real<lower=0> sd_lambda;

    // subjective-level raw parameters, declare as vectors for vectorizing
    vector[N] A_pr;
    vector[N] alpha_pr;
    vector[N] cons_pr;
    vector[N] lambda_pr;
}

transformed parameters {
    vector<lower=0,upper=1>[N] A;
    vector<lower=0,upper=2>[N] alpha;
    vector<lower=0,upper=5>[N] cons;
    vector<lower=0,upper=10>[N] lambda;
    
    for (i in 1:N) {
        A[i]      = Phi_approx( mu_A_pr + sd_A * A_pr[i] );
        alpha[i]  = Phi_approx( mu_alpha_pr + sd_alpha * alpha_pr[i] ) * 2;
        cons[i]   = Phi_approx( mu_cons_pr + sd_cons * cons_pr[i] ) * 5;
        lambda[i] = Phi_approx( mu_lambda_pr + sd_lambda * lambda_pr[i] ) * 10;
    }
}

model {
    // hyperparameters
    mu_A_pr      ~ normal(0,1);
    mu_alpha_pr  ~ normal(0,1);
    mu_cons_pr   ~ normal(0,1);
    mu_lambda_pr ~ normal(0,1);

    sd_A      ~ normal(0,.2);
    sd_alpha  ~ normal(0,.2);
    sd_cons   ~ normal(0,.2);
    sd_lambda ~ normal(0,.2);

    // individual parameters w/ Matt trick
    A_pr      ~ normal(0,1);
    alpha_pr  ~ normal(0,1);
    cons_pr   ~ normal(0,1);
    lambda_pr ~ normal(0,1);
    
    for (i in 1:N) {
        vector[4] ev;
        real curUtil;     // utility of curFb
        real theta;       // theta = 3^c - 1  

        theta = pow(3, cons[i]) -1;
        ev = initV; // initial ev values
        
        for (t in 1:(Tsubj[i]-1)) {
            if ( rewlos[i,t] >= 0) {  // x(t) >= 0
                curUtil = pow(rewlos[i,t], alpha[i]);
            } else {                  // x(t) < 0
                curUtil = -1 * lambda[i] * pow( -1*rewlos[i,t], alpha[i]);
            }
            
            // delta
            ev[ ydata[i, t] ] = ev[ ydata[i, t] ] + A[i]*(curUtil - ev[ ydata[i, t] ]);

            // softmax choice
            ydata[i, t+1] ~ categorical_logit( theta * ev );
        }
    }
}

generated quantities {
    real<lower=0,upper=1> mu_A;
    real<lower=0,upper=2> mu_alpha;
    real<lower=0,upper=5> mu_cons;
    real<lower=0,upper=10> mu_lambda;
    real log_lik[N];

    mu_A      = Phi_approx(mu_A_pr);
    mu_alpha  = Phi_approx(mu_alpha_pr) * 2;
    mu_cons   = Phi_approx(mu_cons_pr) * 5;
    mu_lambda = Phi_approx(mu_lambda_pr) * 10;

    { // local section, this saves time and space
        for (i in 1:N) {
            vector[4] ev;
            real curUtil;     // utility of curFb
            real theta;       // theta = 3^c - 1  
           
            log_lik[i] = 0;
            theta = pow(3, cons[i]) -1;
            ev = initV; // initial ev values
            
            for (t in 1:(Tsubj[i]-1)) {
                if ( rewlos[i,t] >= 0) {  // x(t) >= 0
                    curUtil = pow(rewlos[i,t], alpha[i]);
                } else {                  // x(t) < 0
                    curUtil = -1 * lambda[i] * pow( -1*rewlos[i,t], alpha[i]);
                }
                
                // delta
                ev[ ydata[i, t] ] = ev[ ydata[i, t] ] + A[i]*(curUtil - ev[ ydata[i, t] ]);

                // softmax choice
                log_lik[i] = log_lik[i] + categorical_logit_log( ydata[i, t+1], theta * ev );
            }
        }
    }
}

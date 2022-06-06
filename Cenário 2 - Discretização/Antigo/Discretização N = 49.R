############## Estudo de Simulação ##############
########### Cenário 2 - Discretização ###########
#################### Caso 2 #####################
########## Resolução Média | N = 49 #############  

# Carregando os pacotes necessários -------------
library(MASS)
library(spatstat)
library(matrixcalc)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(beepr)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Geração dos Dados -----------------------------

# 1. Definindo a região S: Região regular [0,10] x [0,10].
# 2. Definindo a discretização em sub-regiões, grade 7x7 com 49 sub-regiões de área 2.040816.
# 3. Calculando os centróides de cada subregião para a matriz de distâncias
x <- rep(seq(5/7, 65/7, by = 10/7), 7)
y <- rep(5/7, 7); for (i in 1:6) y <- c(y, rep(5/7 + (i*10/7), 7))
dados.sim <- as.data.frame(cbind(x, y))

# 4. Calculando a matriz de distâncias entre os centróides
Mdist <- as.matrix(dist(dados.sim))

# 5. Fixando o valor do parâmetro de alcance
max_dist <- max(Mdist)
corr <- 0.04
gamma.sim <- -(max_dist/(2*log(corr)))

# 6. Fixando a média e precisão
mu <- 0
tau.sim <- 1

# 7. Matriz de variância-covariância para o processo Gaussiano w
Sigma.sim <- (1/tau.sim)*exp(-Mdist/gamma.sim)
matrixcalc::is.positive.semi.definite(Sigma.sim)

# 8. Gerando o processo gaussiano w com média constante igual a 0
w <- mvrnorm(n = 1, mu = mu * rep(1, 49), Sigma = Sigma.sim)
dados.sim$w <- w

# 9. Gerando a covariável z
z <- mvrnorm(n = 1, mu = rep(0, 49), Sigma = diag(49))
dados.sim$z <- z

# 10. Fixando o valor do coeficiente beta
beta <- 1.5

# 11. Gerando o valor de offset r
r <- rep(1, 49) # Considerando r(s) = 1, para todo s
dados.sim$r <- r

# 12. Definindo lambda(s) e Lambda(s)
lambda <- exp(z*beta + w)
Lambda <- r * lambda
dados.sim$lambda <- lambda
dados.sim$Lambda <- Lambda

# 13. Gerando as contagens Poisson em cada subregião
N <- rpois(n = 49, lambda = Lambda)
dados.sim$N <- N

k <- length(N)

# Inferência - Stan -----------------------------
0.3/median(Mdist)

cat("
  data {
    int<lower=1> k;
    int N[k];
    vector[k] z;
    vector[k] r;
    matrix<lower=0>[k, k] Mdist;
  }
      
  transformed data {
    vector[k] um;
    for (i in 1:k)
      um[i] = 1;
  }

  parameters {
    real beta;
    real mu;
    real<lower=0> gamma;
    real<lower=0> tau;
    vector[k] w;
  }
  
  model {
    vector[k] lambda_min;
    vector[k] lambda_mai;
    vector[k] mediaVet;
    matrix[k, k] sigmaMat;
    beta ~ normal(0, 100);   
    mu ~ normal(0, 100);
    tau ~ gamma(1, 0.1); 
    gamma ~ gamma(1, 0.05824);
    mediaVet = mu * um;
    sigmaMat = (1/tau) * exp(-Mdist/gamma);
    w ~ multi_normal(mediaVet, sigmaMat);
    for (j in 1:k) {
      lambda_min[j] = exp(z[j]*beta + w[j]);
      lambda_mai[j] = r[j]*lambda_min[j];
      N[j] ~ poisson(lambda_mai[j]);
    }
  }
  
  generated quantities {
    real<lower=0> sigma2 = 1/tau;
  }
  
                     
", file = "Cenário 2 - Discretização/Modelo1Cen2Caso2.stan")

data_1_2_2 <- list(k = k, N = N, z = z, r = r, Mdist = Mdist)

initf2 <- function(chain_id = 1) {
  list(mu = chain_id, beta = chain_id, gamma = chain_id, tau = chain_id, w = rep(chain_id,49))
} 

Mod1Cen2Caso2 <- stan(file = "Cenário 2 - Discretização/Modelo1Cen2Caso2.stan", 
                      data = data_1_2_2, 
                      iter = 20000,
                      warmup = 5000, 
                      thin = 20,
                      chains = 2,
                      init = initf2)
beepr::beep(3)

print(Mod1Cen2Caso2)
res_1_2_2 <- rstan::extract(Mod1Cen2Caso2)
#save.image("Cenário 2 - Discretização/Mod1Cen2Caso2.RData")

# Resultados --------------------
estimativas <- data.frame(cbind(beta = c(beta, mean(res_1_2_2$beta), median(res_1_2_2$beta)),
                                tau = c(tau.sim, mean(res_1_2_2$tau), median(res_1_2_2$tau)),
                                gamma = c(gamma.sim, mean(res_1_2_2$gamma), median(res_1_2_2$gamma)),
                                mu = c(mu, mean(res_1_2_2$mu), median(res_1_2_2$mu))), row.names = c("valor real","média","mediana"))
estimativas %<>% round(digits = 3)
estimativas

DFplots <- data.frame(iter = 1:1500,
                      beta = res_1_2_2$beta,
                      mu = res_1_2_2$mu,
                      tau = res_1_2_2$tau,
                      gamma = res_1_2_2$gamma)

# Beta
rstan::plot(Mod1Cen2Caso2, plotfun = "stan_hist", pars = c("beta"))
rstan::plot(Mod1Cen2Caso2, plotfun = "stan_trace", pars = c("beta"))

ggplot() +
  geom_path(aes(x = iter, y = beta), colour = "blue", data = DFplots) +
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "red", size = 1) +
  xlab("Iterações") +
  ylab(bquote(paste('Parâmetro '*beta*''))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = rel(2)),
        legend.text = element_text(size = rel(2)),
        axis.title.y = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)))

ggplot(aes(x = beta), data = DFplots) +
  geom_histogram(aes(y = ..density..), color = "black", fill = "grey", bins = 30) +
  geom_vline(xintercept = 1.5, linetype = "dashed", color = "red", size = 2) +
  geom_density(alpha = .2, fill="blue") +
  xlab("") +
  ylab("Densidade") +
  theme_bw() +
  theme(axis.title.x = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        axis.title.y = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)))

# Mu
rstan::plot(Mod1Cen2Caso2, plotfun = "stan_hist", pars = c("mu"))
rstan::plot(Mod1Cen2Caso2, plotfun = "stan_trace", pars = c("mu"))

ggplot() +
  geom_path(aes(x = iter, y = mu), colour = "blue", data = DFplots) +
  geom_hline(yintercept=0, linetype = "dashed", color = "red", size = 1) +
  xlab("Iterações") +
  ylab(bquote(paste('Parâmetro '*mu*''))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = rel(2)),
        legend.text = element_text(size = rel(2)),
        axis.title.y = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)))

ggplot(aes(x = mu), data = DFplots) +
  geom_histogram(aes(y = ..density..), color = "black", fill = "grey", bins = 30) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 2) +
  geom_density(alpha = .2, fill="blue") +
  xlab("") +
  ylab("Densidade") +
  theme_bw() +
  theme(axis.title.x = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        axis.title.y = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)))

# Tau
rstan::plot(Mod1Cen2Caso2, plotfun = "stan_hist", pars = c("tau"))
rstan::plot(Mod1Cen2Caso2, plotfun = "stan_trace", pars = c("tau"))

ggplot() +
  geom_path(aes(x = iter, y = tau), colour = "blue", data = DFplots) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) +
  xlab("Iterações") +
  ylab(bquote(paste('Parâmetro '*tau*''))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = rel(2)),
        legend.text = element_text(size = rel(2)),
        axis.title.y = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)))

ggplot(aes(x = tau), data = DFplots) +
  geom_histogram(aes(y = ..density..), color = "black", fill = "grey", bins = 30) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 2) +
  geom_density(alpha = .2, fill="blue") +
  xlab("") +
  ylab("Densidade") +
  theme_bw() +
  theme(axis.title.x = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        axis.title.y = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)))

# Gamma
rstan::plot(Mod1Cen2Caso2, plotfun = "stan_hist", pars = c("gamma"))
rstan::plot(Mod1Cen2Caso2, plotfun = "stan_trace", pars = c("gamma"))

ggplot() +
  geom_path(aes(x = iter, y = gamma), colour = "blue", data = DFplots) +
  geom_hline(yintercept = gamma.sim, linetype = "dashed", color = "red", size = 1) +
  xlab("Iterações") +
  ylab(bquote(paste('Parâmetro '*gamma*''))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = rel(2)),
        legend.text = element_text(size = rel(2)),
        axis.title.y = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)))

ggplot(aes(x = gamma), data = DFplots) +
  geom_histogram(aes(y = ..density..), color = "black", fill = "grey", bins = 30) +
  geom_vline(xintercept = gamma.sim, linetype = "dashed", color = "red", size = 2) +
  geom_density(alpha = .2, fill="blue") +
  xlab("") +
  ylab("Densidade") +
  theme_bw() +
  theme(axis.title.x = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        axis.title.y = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)))

# w
g1w_1_2_2 <- data.frame(x = seq(1,49,1),
                        media = apply(res_1_2_2$w, 2, mean)[seq(1,49,1)],
                        mediana = apply(res_1_2_2$w, 2, median)[seq(1,49,1)],
                        real = dados.sim$w[seq(1,49,1)],
                        inf = apply(res_1_2_2$w, 2, function(x) quantile(x, 0.025))[seq(1,49,1)],
                        sup = apply(res_1_2_2$w, 2, function(x) quantile(x, 0.975))[seq(1,49,1)])


Val.g1w_1_2_2 <- g1w_1_2_2 %>% pivot_longer(cols = 2:6, names_to = "ind", values_to = "valor") %>% filter(ind %in% c("media","real"))

ggplot() +
  geom_point(aes(x = x, y = valor, group = ind, shape = ind, color = ind), size = 2, data = Val.g1w_1_2_2) +
  geom_errorbar(aes(x = x, ymax = sup, ymin = inf), data = g1w_1_2_2) +
  xlab("Sub-regiões") +
  ylab("Intervalo de credibilidade de 95% para w") +
  labs(color = "", shape = "") +
  scale_shape_manual(values=c(17, 16), labels = c("Média", "Valor real")) +
  scale_color_manual(values=c('red','blue'), labels = c("Média", "Valor real")) +
  theme_bw() +
  theme(legend.position = "top",
        axis.title.x = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        axis.title.y = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)))

g2w_1_2_2 <- data.frame(valor.real = dados.sim$w,
                        valor.media = apply(res_1_2_2$w, 2, mean))
ggplot() +
  geom_point(aes(x = valor.real, y = valor.media), size = 2, data = g2w_1_2_2) +
  geom_abline(col = 2, size = 1) +
  xlab("Valor original") +
  ylab("Média") +
  theme_bw() +
  theme(axis.title.x = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        axis.title.y = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)))

# rstan::plot(Mod1Cen2Caso2, plotfun = "stan_dens")
# rstan::plot(Mod1Cen2Caso2, plotfun = "stan_trace")
# rstan::plot(Mod1Cen2Caso2, plotfun = "stan_hist")
# rstan::plot(Mod1Cen2Caso2, plotfun = "stan_diag")
# rstan::plot(Mod1Cen2Caso2, plotfun = "stan_scat", pars = c("gamma","tau"))
# rstan::plot(Mod1Cen2Caso2, plotfun = "stan_rhat")
# rstan::plot(Mod1Cen2Caso2, plotfun = "stan_ess")
# rstan::plot(Mod1Cen2Caso2, plotfun = "stan_mcse")

#Gráfico do padrão de pontos resultante
Sk <- list(NULL)
S <- data.frame(NULL)
for (i in 1:7) {
  for (j in 1:7) {
    Sk[[7*(i-1) + j]] <- runifpoint(dados.sim$N[7*(i-1) + j], win = owin(c((j-1)*10/7,j*10/7), c((i-1)*10/7,i*10/7)))
    S <- rbind(S, as.data.frame(Sk[[7*(i-1) + j]])) 
  }
}
S <- as.ppp(S, W = owin(c(0,10), c(0,10)))
plot(S, pch = 20, main="")
qS <- quadratcount(S, 7, 7)
plot(qS, add = T, col = "", cex = 1.5, lty = 2)
dados.sim$N

#load("Cenário 2 - Discretização/Mod1Cen2Caso2.RData")
#save.image("Cenário 2 - Discretização/Mod1Cen2Caso2.RData")

#Salvando objetos para comparação entre simulações
dados.sim_1_2_2 <- dados.sim
gamma.sim_1_2_2 <- gamma.sim
rm(list = ls()[!(ls() %in% c("res_1_2_2","dados.sim_1_2_2","gamma.sim_1_2_2"))])

#load("Cenário 2 - Discretização/Comp.RData")
#save.image("Cenário 2 - Discretização/Comp.RData")

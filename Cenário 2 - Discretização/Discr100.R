#################################### Estudo de Simulação #####################################
################################# Cenário 2 - Discretização ##################################
########################################## Caso 1 ############################################
################################# Grade mais fina N = 100 ####################################

# Carregando os pacotes necessários ----------------------------------------------------------
library(MASS)
library(spatstat)
library(matrixcalc)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtext)
library(beepr)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Geração dos Dados -----------------------------

# 1. Definindo a região S: Região regular [0,10] x [0,10].
# 2. Definindo a discretização em sub-regiões, grade 10x10 com 100 sub-regiões de área 1.
# 3. Calculando os centroides de cada sub-região para a matriz de distâncias
x <- rep(seq(0.5, 9.5, by = 1), 10)
y <- rep(0.5, 10); for (i in 1:9) y <- c(y, rep(0.5 + i, 10))
dados.sim <- as.data.frame(cbind(x, y))

# 4. Calculando a matriz de distâncias entre os centroides
Mdist <- as.matrix(dist(dados.sim))

# 5. Fixando o valor do parâmetro de alcance
max_dist <- max(Mdist)
corr <- 0.04
gamma.sim <- -(max_dist/(2*log(corr)))

# 6. Fixando a média e precisão
mu <- 0
tau.sim <- 2

# 7. Matriz de variância-covariância para o processo Gaussiano w
Sigma.sim <- (1/tau.sim)*exp(-Mdist/gamma.sim)
matrixcalc::is.positive.semi.definite(Sigma.sim)

# 8. Gerando o processo gaussiano w com média constante igual a 0
set.seed(NULL)
w <- mvrnorm(n = 1, mu = mu * rep(1, 100), Sigma = Sigma.sim)
dados.sim$w <- w
plot(dados.sim$w, type = "l")
abline(h = 0, col = "red")

# 9. Gerando a covariável z
set.seed(2021)
z <- mvrnorm(n = 1, mu = rep(0, 100), Sigma = diag(100))
dados.sim$z <- z

# 10. Fixando o valor do coeficiente beta
beta <- 1.65

# 11. Gerando o valor de offset r
r <- rep(1, 100) # Considerando r(s) = 1, para todo s
dados.sim$r <- r

# 12. Definindo lambda(s) e Lambda(s)
lambda <- exp(z*beta + w)
Lambda <- r * lambda
dados.sim$lambda <- lambda
dados.sim$Lambda <- Lambda

# 13. Gerando as contagens Poisson em cada sub-região
set.seed(NULL)
N <- rpois(n = 100, lambda = Lambda)
dados.sim$N <- N

k <- length(N)
sum(dados.sim$N)
dados.sim$N
mean(dados.sim$w)

# Inferência - Stan --------------------------------------------------------------------------
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
    gamma ~ gamma(1, 0.05883);
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
  
                     
", file = "Cenário 1 - Alcance/Modelo1Cen1Caso2.stan")

data_1_1_2 <- list(k = k, N = N, z = z, r = r, Mdist = Mdist)

initf2 <- function(chain_id = 1) {
  list(mu = chain_id, beta = chain_id, gamma = chain_id, tau = chain_id, w = rep(chain_id,100))
} 

Mod1Cen1Caso2 <- stan(file = "Cenário 1 - Alcance/Modelo1Cen1Caso2.stan", 
                      data = data_1_1_2, 
                      iter = 20000,
                      warmup = 5000, 
                      thin = 20,
                      chains = 2,
                      init = initf2)
beepr::beep(3)

print(Mod1Cen1Caso2)
res_1_1_2 <- rstan::extract(Mod1Cen1Caso2)
#save.image("Cenário 1 - Alcance/Mod1Cen1Caso2.RData")

# Resultados ---------------------------------------------------------------------------------
estimativas <- data.frame(cbind(beta = c(beta, mean(res_1_1_2$beta), median(res_1_1_2$beta)),
                                tau = c(tau.sim, mean(res_1_1_2$tau), median(res_1_1_2$tau)),
                                gamma = c(gamma.sim, mean(res_1_1_2$gamma), 
                                          median(res_1_1_2$gamma)),
                                mu = c(mu, mean(res_1_1_2$mu), median(res_1_1_2$mu))), 
                          row.names = c("valor real","média","mediana"))
estimativas %<>% round(digits = 3)
estimativas

DFplots <- data.frame(iter = 1:1500,
                      beta = res_1_1_2$beta,
                      mu = res_1_1_2$mu,
                      tau = res_1_1_2$tau,
                      gamma = res_1_1_2$gamma)

# Beta
rstan::plot(Mod1Cen1Caso2, plotfun = "stan_hist", pars = c("beta"))
rstan::plot(Mod1Cen1Caso2, plotfun = "stan_trace", pars = c("beta"))

ggplot() +
  geom_path(aes(x = iter, y = beta), colour = "blue", data = DFplots) +
  geom_hline(yintercept=1.65, linetype = "dashed", color = "red", size = 2) +
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
  geom_vline(xintercept = 1.65, linetype = "dashed", color = "red", size = 2) +
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
rstan::plot(Mod1Cen1Caso2, plotfun = "stan_hist", pars = c("mu"))
rstan::plot(Mod1Cen1Caso2, plotfun = "stan_trace", pars = c("mu"))

ggplot() +
  geom_path(aes(x = iter, y = mu), colour = "blue", data = DFplots) +
  geom_hline(yintercept=0, linetype = "dashed", color = "red", size = 2) +
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
rstan::plot(Mod1Cen1Caso2, plotfun = "stan_hist", pars = c("tau"))
rstan::plot(Mod1Cen1Caso2, plotfun = "stan_trace", pars = c("tau"))

ggplot() +
  geom_path(aes(x = iter, y = tau), colour = "blue", data = DFplots) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red", size = 2) +
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
  geom_vline(xintercept = 2, linetype = "dashed", color = "red", size = 2) +
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
rstan::plot(Mod1Cen1Caso2, plotfun = "stan_hist", pars = c("gamma"))
rstan::plot(Mod1Cen1Caso2, plotfun = "stan_trace", pars = c("gamma"))

ggplot() +
  geom_path(aes(x = iter, y = gamma), colour = "blue", data = DFplots) +
  geom_hline(yintercept = gamma.sim, linetype = "dashed", color = "red", size = 2) +
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
g1w_1_1_2 <- data.frame(x = seq(1,100,1),
                        media = apply(res_1_1_2$w, 2, mean)[seq(1,100,1)],
                        mediana = apply(res_1_1_2$w, 2, median)[seq(1,100,1)],
                        real = dados.sim$w[seq(1,100,1)],
                        inf = apply(res_1_1_2$w, 2, 
                                    function(x) quantile(x, 0.025))[seq(1,100,1)],
                        sup = apply(res_1_1_2$w, 2, 
                                    function(x) quantile(x, 0.975))[seq(1,100,1)])


Val.g1w_1_1_2 <- g1w_1_1_2 %>% 
                  pivot_longer(cols = 2:6, names_to = "ind", values_to = "valor") %>% 
                  filter(ind %in% c("media","real"))

ggplot() +
  geom_point(aes(x = x, y = valor, group = ind, shape = ind, color = ind), 
             size = 2, data = Val.g1w_1_1_2) +
  geom_errorbar(aes(x = x, ymax = sup, ymin = inf), width = 0, data = g1w_1_1_2) +
  xlab("Sub-regiões") +
  ylab("Intervalo de credibilidade de 95% para w") +
  labs(color = "", shape = "") +
  scale_shape_manual(values=c(17, 16), labels = c("Média*a posteriori*", "Valor real")) +
  scale_color_manual(values=c('red','blue'), labels = c("Média*a posteriori*", "Valor real")) +
  coord_cartesian(ylim = c(-3, 3)) +
  theme_bw() +
  theme(legend.position = "top",
        axis.title.x = element_text(size = rel(1.5)),
        legend.text = element_markdown(size = rel(1.5)),
        legend.text.align = 0,
        axis.title.y = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)))

ggplot() +
  geom_point(aes(x = real, y = media), size = 2, data = g1w_1_1_2) +
  geom_abline(col = 2, size = 1) +
  xlab("Valor original") +
  ylab("Média") +
  coord_cartesian(ylim = c(-2, 2), xlim = c(-2,2)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        axis.title.y = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)))

# rstan::plot(Mod1Cen1Caso2, plotfun = "stan_dens")
# rstan::plot(Mod1Cen1Caso2, plotfun = "stan_trace")
# rstan::plot(Mod1Cen1Caso2, plotfun = "stan_hist")
# rstan::plot(Mod1Cen1Caso2, plotfun = "stan_diag")
# rstan::plot(Mod1Cen1Caso2, plotfun = "stan_scat", pars = c("gamma","tau"))
# rstan::plot(Mod1Cen1Caso2, plotfun = "stan_rhat")
# rstan::plot(Mod1Cen1Caso2, plotfun = "stan_ess")
# rstan::plot(Mod1Cen1Caso2, plotfun = "stan_mcse")

#Gráfico do padrão de pontos resultante
Sk <- list(NULL)
S <- data.frame(NULL)
for (i in 1:10) {
  for (j in 1:10) {
    Sk[[10*(i-1) + j]] <- runifpoint(dados.sim$N[10*(i-1) + j], win = owin(c(j-1,j), c(i-1,i)))
    S <- rbind(S, as.data.frame(Sk[[10*(i-1) + j]])) 
  }
}
S <- as.ppp(S, W = owin(c(0,10), c(0,10)))
plot(S, pch = 20, main="")
qS <- quadratcount(S, 10, 10)
plot(qS, add = T, col = "", cex = 1.5, lty = 2)
dados.sim$N

#load("Cenário 1 - Alcance/Mod1Cen1Caso2.RData")
#save.image("Cenário 1 - Alcance/Mod1Cen1Caso2.RData")

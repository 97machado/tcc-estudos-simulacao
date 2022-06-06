#################################### Estudo de Simulação #####################################
################################# Cenário 2 - Discretização ##################################
########################################## Caso 1 ############################################
################################# Grade mais grossa N = 25 # (mu = 0) fixado #################

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

# Carregando os dados gerados para grade com 100 sub-regiões ---------------------------------
#load("Cenário 1 - Alcance/Mod1Cen1Caso2.RData")

# Salvando os dados
#save.image("Cenário 2 - Discretização/Discr100.RData")

# Geração das sub-regiões --------------------------------------------------------------------

# 1. Definindo a região S: Região regular [0,10] x [0,10].
# 2. Definindo a discretização em sub-regiões, grade 5x5 com 25 sub-regiões de área 4.
# 3. Calculando os centroides de cada sub-região para a matriz de distâncias
x_25 <- rep(seq(1, 9, by = 2), 5)
y_25 <- rep(1, 5); for (i in 1:4) y_25 <- c(y_25, rep(1 + (i*2), 5))
dados.sim_25 <- as.data.frame(cbind(x_25, y_25))

# 4. Calculando a matriz de distâncias entre os novos centroides
Mdist_25 <- as.matrix(dist(dados.sim_25))

# 5. Gerando nova covariável z de dimensão 25
# set.seed(2021)
# z_25 <- mvrnorm(n = 1, mu = rep(0, 25), Sigma = diag(25))
# set.seed(NULL)
# dados.sim_25$z_25 <- z_25

z_25 <- c(NULL)
for (i in 1:25) {
  z_25[i] <- (dados.sim$z[(i + (i-1)) + 10*floor((i-1)/5)] + 
                dados.sim$z[(i+1 + (i-1)) + 10*floor((i-1)/5)] + 
                dados.sim$z[(i+10 + (i-1)) + 10*floor((i-1)/5)] + 
                dados.sim$z[(i+11 + (i-1)) + 10*floor((i-1)/5)])/4
}
dados.sim_25$z_25 <- z_25



# 6. Gerando o valor de offset r de tamanho 25
r_25 <- rep(1, 25) # Considerando r(s) = 1, para todo s
dados.sim_25$r_25 <- r_25

# 7. Recalculando o número de pontos em cada sub-região
N_25 <- c(NULL)
for (i in 1:25) {
  N_25[i] <- dados.sim$N[(i + (i-1)) + 10*floor((i-1)/5)] + 
    dados.sim$N[(i+1 + (i-1)) + 10*floor((i-1)/5)] + 
    dados.sim$N[(i+10 + (i-1)) + 10*floor((i-1)/5)] + 
    dados.sim$N[(i+11 + (i-1)) + 10*floor((i-1)/5)]
}
dados.sim$N_25 <- N_25

k_25 <- length(N_25)

# Inferência - Stan --------------------------------------------------------------------------
0.3/median(Mdist_25)

cat("
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
  
                     
", file = "Cenário 2 - Discretização/Discr25_Mu0.stan")

data_Discr25_Mu0 <- list(k = k_25, N = N_25, z = z_25, r = r_25, Mdist = Mdist_25)

initf2_25Mu0 <- function(chain_id = 1) {
  list(beta = chain_id, gamma = chain_id, tau = chain_id, w = rep(chain_id,25))
} 

ModDiscr25_Mu0 <- stan(file = "Cenário 2 - Discretização/Discr25_Mu0.stan", 
                   data = data_Discr25_Mu0, 
                   iter = 20000,
                   warmup = 5000, 
                   thin = 20,
                   chains = 2,
                   init = initf2_25Mu0)
beepr::beep(3)

print(ModDiscr25_Mu0)
res_Discr25_Mu0 <- rstan::extract(ModDiscr25_Mu0)
save.image("Cenário 2 - Discretização/Discr25_Mu0.RData")

# Resultados ---------------------------------------------------------------------------------
estimativas_25 <- data.frame(cbind(beta = c(beta, mean(res_Discr25_Mu0$beta), 
                                            median(res_Discr25_Mu0$beta)),
                                   tau = c(tau.sim, mean(res_Discr25_Mu0$tau), 
                                           median(res_Discr25_Mu0$tau)),
                                   gamma = c(gamma.sim, mean(res_Discr25_Mu0$gamma), 
                                             median(res_Discr25_Mu0$gamma)),
                                   mu = c(mu, mean(res_Discr25_Mu0$mu), 
                                          median(res_Discr25_Mu0$mu))), 
                             row.names = c("valor real","média","mediana"))
estimativas_25 %<>% round(digits = 3)
estimativas_25

DFplots_25 <- data.frame(iter = 1:1500,
                         beta = res_Discr25_Mu0$beta,
                         mu = res_Discr25_Mu0$mu,
                         tau = res_Discr25_Mu0$tau,
                         gamma = res_Discr25_Mu0$gamma)

# Beta
rstan::plot(ModDiscr25_Mu0, plotfun = "stan_hist", pars = c("beta"))
rstan::plot(ModDiscr25_Mu0, plotfun = "stan_trace", pars = c("beta"))

ggplot() +
  geom_path(aes(x = iter, y = beta), colour = "blue", data = DFplots_25) +
  geom_hline(yintercept=1.65, linetype = "dashed", color = "red", size = 1) +
  xlab("Iterações") +
  ylab(bquote(paste('Parâmetro '*beta*''))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = rel(2)),
        legend.text = element_text(size = rel(2)),
        axis.title.y = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)))

ggplot(aes(x = beta), data = DFplots_25) +
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
rstan::plot(ModDiscr25_Mu0, plotfun = "stan_hist", pars = c("mu"))
rstan::plot(ModDiscr25_Mu0, plotfun = "stan_trace", pars = c("mu"))

ggplot() +
  geom_path(aes(x = iter, y = mu), colour = "blue", data = DFplots_25) +
  geom_hline(yintercept=0, linetype = "dashed", color = "red", size = 1) +
  xlab("Iterações") +
  ylab(bquote(paste('Parâmetro '*mu*''))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = rel(2)),
        legend.text = element_text(size = rel(2)),
        axis.title.y = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)))

ggplot(aes(x = mu), data = DFplots_25) +
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
rstan::plot(ModDiscr25_Mu0, plotfun = "stan_hist", pars = c("tau"))
rstan::plot(ModDiscr25_Mu0, plotfun = "stan_trace", pars = c("tau"))

ggplot() +
  geom_path(aes(x = iter, y = tau), colour = "blue", data = DFplots_25) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red", size = 1) +
  xlab("Iterações") +
  ylab(bquote(paste('Parâmetro '*tau*''))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = rel(2)),
        legend.text = element_text(size = rel(2)),
        axis.title.y = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)))

ggplot(aes(x = tau), data = DFplots_25) +
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
rstan::plot(ModDiscr25_Mu0, plotfun = "stan_hist", pars = c("gamma"))
rstan::plot(ModDiscr25_Mu0, plotfun = "stan_trace", pars = c("gamma"))

ggplot() +
  geom_path(aes(x = iter, y = gamma), colour = "blue", data = DFplots_25) +
  geom_hline(yintercept = gamma.sim, linetype = "dashed", color = "red", size = 1) +
  xlab("Iterações") +
  ylab(bquote(paste('Parâmetro '*gamma*''))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = rel(2)),
        legend.text = element_text(size = rel(2)),
        axis.title.y = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)))

ggplot(aes(x = gamma), data = DFplots_25) +
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
DFw25 <- data.frame(x = 1:25,
                    media = apply(res_Discr25_Mu0$w, 2, mean),
                    inf = apply(res_Discr25_Mu0$w, 2, function(x) 
                      quantile(x, 0.025)),
                    sup = apply(res_Discr25_Mu0$w, 2, function(x) 
                      quantile(x, 0.975)))
#LngDFw25 <- DFw25 %>% pivot_longer(cols = 2:4, names_to = "ind", values_to = "valor")

ordem <- c(NULL)
for (i in 1:25) {
  ordem <- c(ordem, 
             (i + (i-1)) + 10*floor((i-1)/5), 
             (i+1 + (i-1)) + 10*floor((i-1)/5), 
             (i+10 + (i-1)) + 10*floor((i-1)/5), 
             (i+11 + (i-1)) + 10*floor((i-1)/5))
}

DFw100 <- data.frame(x = seq(0.125,24.875, length.out = 100),
                     media = apply(res_1_1_2$w, 2, mean)[ordem],
                     inf = apply(res_1_1_2$w, 2, function(x) 
                       quantile(x, 0.025))[ordem],
                     sup = apply(res_1_1_2$w, 2, function(x) 
                       quantile(x, 0.975))[ordem])

#LngDFw100 <- DFw100 %>% pivot_longer(cols = 2:4, names_to = "ind", values_to = "valor")

ggplot() +
  geom_rect(aes(ymin = inf, ymax = sup, xmin = x-1, xmax = x), alpha = 0.4, 
            data = DFw25) +  
  geom_errorbar(aes(x = x, ymax = sup, ymin = inf), width = 0, data = DFw100) +
  geom_segment(aes(x = x-1, xend = x, y = media, yend = media), data = DFw25,
               linetype = "dashed", col = "blue", size = 0.9) +
  geom_point(aes(x = x, y = media), size = 2, 
             data = DFw100, color = "red") +
  geom_vline(aes(xintercept = 0:25), linetype = "dashed") +
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


# g1w_Discr25_Mu0 <- data.frame(x = seq(1,25,1),
#                         media = apply(res_Discr25_Mu0$w, 2, mean)[seq(1,25,1)],
#                         mediana = apply(res_Discr25_Mu0$w, 2, median)[seq(1,25,1)],
#                         #real = dados.sim$w[seq(1,25,1)],
#                         inf = apply(res_Discr25_Mu0$w, 2, function(x) 
#                                       quantile(x, 0.025))[seq(1,25,1)],
#                         sup = apply(res_Discr25_Mu0$w, 2, function(x) 
#                                       quantile(x, 0.975))[seq(1,25,1)])
# 
# 
# Val.g1w_Discr25_Mu0 <- g1w_Discr25_Mu0 %>% 
#                       pivot_longer(cols = 2:6, names_to = "ind", values_to = "valor") %>% 
#                       filter(ind %in% c("media","real"))
# 
# ggplot() +
#   geom_point(aes(x = x, y = valor, group = ind, shape = ind, color = ind), size = 2, 
#              data = Val.g1w_Discr25_Mu0) +
#   geom_errorbar(aes(x = x, ymax = sup, ymin = inf), data = g1w_Discr25_Mu0) +
#   xlab("Sub-regiões") +
#   ylab("Intervalo de credibilidade de 95% para w") +
#   labs(color = "", shape = "") +
#   scale_shape_manual(values=c(17, 16), labels = c("Média", "Valor real")) +
#   scale_color_manual(values=c('red','blue'), labels = c("Média", "Valor real")) +
#   theme_bw() +
#   theme(legend.position = "top",
#         axis.title.x = element_text(size = rel(1.5)),
#         legend.text = element_text(size = rel(1.5)),
#         axis.title.y = element_text(size = rel(1.5)),
#         axis.text.y = element_text(size = rel(2)),
#         axis.text.x = element_text(size = rel(2)))
# 
# g2w_Discr25_Mu0 <- data.frame(valor.real = dados.sim$w,
#                         valor.media = apply(res_Discr25_Mu0$w, 2, mean))
# ggplot() +
#   geom_point(aes(x = valor.real, y = valor.media), size = 2, data = g2w_Discr25_Mu0) +
#   geom_abline(col = 2, size = 1) +
#   xlab("Valor original") +
#   ylab("Média") +
#   theme_bw() +
#   theme(axis.title.x = element_text(size = rel(1.5)),
#         legend.text = element_text(size = rel(1.5)),
#         axis.title.y = element_text(size = rel(1.5)),
#         axis.text.y = element_text(size = rel(2)),
#         axis.text.x = element_text(size = rel(2)))


#Gráfico do padrão de pontos resultante
Sk <- list(NULL)
S <- data.frame(NULL)
for (i in 1:5) {
  for (j in 1:5) {
    Sk[[5*(i-1) + j]] <- runifpoint(dados.sim$N[5*(i-1) + j], win = owin(c((j-1)*2,j*2), c((i-1)*2,i*2)))
    S <- rbind(S, as.data.frame(Sk[[5*(i-1) + j]])) 
  }
}
S <- as.ppp(S, W = owin(c(0,10), c(0,10)))
qS_25 <- quadratcount(S, 5, 5)
plot(S, pch = 20, main="")
plot(qS, add = T, col = "2", cex = 1.5, lty = 2)
plot(S, pch = 20, main="")
plot(qS_25, add = T, col = "2", cex = 1.5, lty = 2)
dados.sim$N
N_25
#load("Cenário 2 - Discretização/Discr25_Mu0.RData")
#save.image("Cenário 2 - Discretização/Discr25_Mu0.RData")

#Salvando objetos para comparação entre simulações
dados.sim_Discr25_Mu0 <- dados.sim
gamma.sim_Discr25_Mu0 <- gamma.sim
rm(list = ls()[!(ls() %in% c("res_Discr25_Mu0","dados.sim_Discr25_Mu0","gamma.sim_Discr25_Mu0"))])

#save.image("Cenário 2 - Discretização/Comp.RData")

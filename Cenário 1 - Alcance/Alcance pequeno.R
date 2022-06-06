#################################### Estudo de Simulação #####################################
################################### Cenário 1 - Alcance ######################################
########################################## Caso 1 ############################################
################################# Alcance Pequeno - 0.01 #####################################

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

# Geração dos Dados --------------------------------------------------------------------------

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
corr <- 0.01
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
beta <- 1.8

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
  
                     
", file = "Cenário 1 - Alcance/Modelo1Cen1Caso1.stan")

data_1_1_1 <- list(k = k, N = N, z = z, r = r, Mdist = Mdist)

initf2 <- function(chain_id = 1) {
  list(mu = chain_id, beta = chain_id, gamma = chain_id, tau = chain_id, w = rep(chain_id,100))
} 

Mod1Cen1Caso1 <- stan(file = "Cenário 1 - Alcance/Modelo1Cen1Caso1.stan", 
                       data = data_1_1_1, 
                       iter = 20000,
                       warmup = 5000, 
                       thin = 20,
                       chains = 2,
                       init = initf2)
beepr::beep(3)

print(Mod1Cen1Caso1)
res_1_1_1 <- rstan::extract(Mod1Cen1Caso1)
save.image("Cenário 1 - Alcance/Mod1Cen1Caso1.RData")

# Resultados ---------------------------------------------------------------------------------
estimativas <- data.frame(cbind(beta = c(beta, mean(res_1_1_1$beta), median(res_1_1_1$beta)),
                                tau = c(tau.sim, mean(res_1_1_1$tau), median(res_1_1_1$tau)),
                                gamma = c(gamma.sim, mean(res_1_1_1$gamma), 
                                          median(res_1_1_1$gamma)),
                                mu = c(mu, mean(res_1_1_1$mu), median(res_1_1_1$mu))),
                          row.names = c("valor real","média","mediana"))
estimativas %<>% round(digits = 3)
estimativas

DFplots <- data.frame(iter = 1:1500,
                      beta = res_1_1_1$beta,
                      mu = res_1_1_1$mu,
                      tau = res_1_1_1$tau,
                      gamma = res_1_1_1$gamma)

# Beta
rstan::plot(Mod1Cen1Caso1, plotfun = "stan_hist", pars = c("beta"))
rstan::plot(Mod1Cen1Caso1, plotfun = "stan_trace", pars = c("beta"))

ggplot() +
  geom_path(aes(x = iter, y = beta), colour = "blue", data = DFplots) +
  geom_hline(yintercept=1.8, linetype = "dashed", color = "red", size = 2) +
  xlab("Iterações") +
  ylab(bquote(paste('Parâmetro'))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = rel(2)),
        legend.text = element_text(size = rel(2)),
        axis.title.y = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)))

ggplot(aes(x = beta), data = DFplots) +
  geom_histogram(aes(y = ..density..), color = "black", fill = "grey", bins = 30) +
  geom_vline(xintercept = 1.8, linetype = "dashed", color = "red", size = 2) +
  geom_density(alpha = .2, fill="blue") +
  xlab("") +
  ylab("Densidade") +
  coord_cartesian(ylim = c(0, 4.5)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = rel(2.5)),
        legend.text = element_text(size = rel(1.5)),
        axis.title.y = element_text(size = rel(2.5)),
        axis.text.y = element_text(size = rel(3)),
        axis.text.x = element_text(size = rel(3)))

# Mu
rstan::plot(Mod1Cen1Caso1, plotfun = "stan_hist", pars = c("mu"))
rstan::plot(Mod1Cen1Caso1, plotfun = "stan_trace", pars = c("mu"))

ggplot() +
  geom_path(aes(x = iter, y = mu), colour = "blue", data = DFplots) +
  geom_hline(yintercept=0, linetype = "dashed", color = "red", size = 2) +
  xlab("Iterações") +
  ylab(bquote(paste('Parâmetro'))) +
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
  coord_cartesian(ylim = c(0, 1.5)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = rel(2.5)),
        legend.text = element_text(size = rel(1.5)),
        axis.title.y = element_text(size = rel(2.5)),
        axis.text.y = element_text(size = rel(3)),
        axis.text.x = element_text(size = rel(3)))

# Tau
rstan::plot(Mod1Cen1Caso1, plotfun = "stan_hist", pars = c("tau"))
rstan::plot(Mod1Cen1Caso1, plotfun = "stan_trace", pars = c("tau"))

ggplot() +
  geom_path(aes(x = iter, y = tau), colour = "blue", data = DFplots) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red", size = 2) +
  xlab("Iterações") +
  ylab(bquote(paste('Parâmetro'))) +
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
  coord_cartesian(ylim = c(0, 0.45)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = rel(2.5)),
        legend.text = element_text(size = rel(1.5)),
        axis.title.y = element_text(size = rel(2.5)),
        axis.text.y = element_text(size = rel(3)),
        axis.text.x = element_text(size = rel(3)))

# Gamma
rstan::plot(Mod1Cen1Caso1, plotfun = "stan_hist", pars = c("gamma"))
rstan::plot(Mod1Cen1Caso1, plotfun = "stan_trace", pars = c("gamma"))

ggplot() +
  geom_path(aes(x = iter, y = gamma), colour = "blue", data = DFplots) +
  geom_hline(yintercept = gamma.sim, linetype = "dashed", color = "red", size = 2) +
  xlab("Iterações") +
  ylab(bquote(paste('Parâmetro'))) +
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
  coord_cartesian(ylim = c(0, 0.35)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = rel(2.5)),
        legend.text = element_text(size = rel(1.5)),
        axis.title.y = element_text(size = rel(2.5)),
        axis.text.y = element_text(size = rel(3)),
        axis.text.x = element_text(size = rel(3)))

# Sorteando 4 w's para exibir as cadeias 

# Cadeia Fraca
set.seed(1); sample(1:100, size = 4) # w[1], w[34], w[39] e w[68]
DFplots_w <- data.frame(iter = 1:1500,
                        w1 = res_1_1_1$w[,1],
                        w34 = res_1_1_1$w[,34],
                        w39 = res_1_1_1$w[,39],
                        w68 = res_1_1_1$w[,68])

# w[1]
ggplot() +
  geom_path(aes(x = iter, y = w1), colour = "blue", data = DFplots_w) +
  geom_hline(yintercept = dados.sim$w[1], linetype = "dashed", color = "red", size = 2) +
  xlab("Iterações") +
  ylab(bquote(paste('Processo Gaussiano w'['1']*''))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = rel(2)),
        legend.text = element_text(size = rel(2)),
        axis.title.y = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)))

# w[34]
ggplot() +
  geom_path(aes(x = iter, y = w34), colour = "blue", data = DFplots_w) +
  geom_hline(yintercept = dados.sim$w[34], linetype = "dashed", color = "red", size = 2) +
  xlab("Iterações") +
  ylab(bquote(paste('Processo Gaussiano w'['34']*''))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = rel(2)),
        legend.text = element_text(size = rel(2)),
        axis.title.y = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)))

# w[39]
ggplot() +
  geom_path(aes(x = iter, y = w39), colour = "blue", data = DFplots_w) +
  geom_hline(yintercept = dados.sim$w[39], linetype = "dashed", color = "red", size = 2) +
  xlab("Iterações") +
  ylab(bquote(paste('Processo Gaussiano w'['39']*''))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = rel(2)),
        legend.text = element_text(size = rel(2)),
        axis.title.y = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)))

# w[68]
ggplot() +
  geom_path(aes(x = iter, y = w68), colour = "blue", data = DFplots_w) +
  geom_hline(yintercept = dados.sim$w[68], linetype = "dashed", color = "red", size = 2) +
  xlab("Iterações") +
  ylab(bquote(paste('Processo Gaussiano w'['68']*''))) +
  theme_bw() +
  theme(axis.title.x = element_text(size = rel(2)),
        legend.text = element_text(size = rel(2)),
        axis.title.y = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)))

# w
g1w_1_1_1 <- data.frame(x = seq(1,100,1),
                        media = apply(res_1_1_1$w, 2, mean)[seq(1,100,1)],
                        mediana = apply(res_1_1_1$w, 2, median)[seq(1,100,1)],
                        real = dados.sim$w[seq(1,100,1)],
                        inf = apply(res_1_1_1$w, 2, 
                                    function(x) quantile(x, 0.025))[seq(1,100,1)],
                        sup = apply(res_1_1_1$w, 2, 
                                    function(x) quantile(x, 0.975))[seq(1,100,1)])

Val.g1w_1_1_1 <- g1w_1_1_1 %>% 
                  pivot_longer(cols = 2:6, names_to = "ind", values_to = "valor") %>% 
                  filter(ind %in% c("media","real"))

ggplot() +
  geom_point(aes(x = x, y = valor, group = ind, shape = ind, color = ind), 
             size = 2, data = Val.g1w_1_1_1) +
  geom_errorbar(aes(x = x, ymax = sup, ymin = inf), width = 0, data = g1w_1_1_1) +
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
  geom_point(aes(x = real, y = media), size = 2, data = g1w_1_1_1) +
  geom_abline(col = 2, size = 1) +
  xlab("Valor original") +
  ylab("Média") +
  coord_cartesian(ylim = c(-2, 2), xlim = c(-2,2)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = rel(2.5)),
        legend.text = element_text(size = rel(1.5)),
        axis.title.y = element_text(size = rel(2.5)),
        axis.text.y = element_text(size = rel(3)),
        axis.text.x = element_text(size = rel(3)))

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

# Mapa das covariáveis e PG
DFMapa <- data.frame(N = dados.sim$N, z = dados.sim$z, w = dados.sim$w, 
                     cx = dados.sim$x, cy = dados.sim$y)

# ggplot() +
#   geom_point(aes(x = x, y = y), data = as.data.frame(S)) +
#   geom_segment(aes(x = 0, xend = 10, y = 1:9, yend = 1:9),
#                linetype = "dashed", col = "black", size = 0.7) +
#   geom_segment(aes(x = 1:9, xend = 1:9, y = 0, yend = 10),
#                linetype = "dashed", col = "black", size = 0.7) +
#   geom_segment(aes(x = 0, xend = 10, y = c(0,10), yend = c(0,10)),
#                col = "black", size = 0.7) +
#   geom_segment(aes(x = c(0,10), xend = c(0,10), y = 0, yend = 10),
#                col = "black", size = 0.7) +
#   geom_point(aes(x = cx, y = cy, col = z), data = DFMapa, size = 10,
#              alpha = .2) +
#   scale_size_continuous(range = c(1, 10)) +
#   coord_cartesian(ylim = c(0, 10), xlim = c(0,10)) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         panel.border = element_blank(),
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank())

im_z <- DFMapa %>% select(cx, cy, z) %>% as.im()
im_w <- DFMapa %>% select(cx, cy, w) %>% as.im()

g1w_1_1_1$cx <- dados.sim$x
g1w_1_1_1$cy <- dados.sim$y
g1w_1_1_1$contido <- ifelse((g1w_1_1_1$real >= g1w_1_1_1$inf) & 
                            (g1w_1_1_1$real <= g1w_1_1_1$sup), 1, 0)
table(g1w_1_1_1$contido)
ForaICFraca <- g1w_1_1_1 %>% filter(contido == 0) %>% select(cx, cy)  %>% 
                             as.ppp(W = owin(c(0,10), c(0,10)))

# par(mfrow=c(1,2))
# Mapa de z
plot(im_z, main = "", col = heat.colors(128, rev = TRUE))
plot(qS, add = T, col = "", cex = 1.5, lty = 2, lwd = 2)
plot(S, add = T, cex = 1.5, lty = 2, pch = 20, lwd = 1)

# Mapa de w
plot(im_w, main = "", col = terrain.colors(128))
plot(qS, add = T, col = "", cex = 1.5, lty = 2, lwd = 2)
plot(S, add = T, cex = 1.5, lty = 2, pch = 20, lwd = 1)
plot(ForaICFraca, add = T, pch = "x", col = "red", cex = 1.5)

#load("Cenário 1 - Alcance/Mod1Cen1Caso1.RData")
save.image("Cenário 1 - Alcance/Mod1Cen1Caso1.RData")

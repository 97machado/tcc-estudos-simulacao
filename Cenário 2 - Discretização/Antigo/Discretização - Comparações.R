#################################### Estudo de Simulação #######################################
################################# Cenário 2 - Discretização ####################################
####################################### Comparações ############################################

# Carregando os pacotes necessários ------------------------------------------------------------
#library(MASS)
# library(spatstat)
# library(matrixcalc)
#library(magrittr)
#library(dplyr)
#library(tidyr)
library(ggplot2)
library(gridExtra)
#library(rstan)

#load("Cenário 2 - Discretização/Comp.RData")

# Comparações ----------------------------------------------------------------------------------

# w
g2w_discr_1 <- data.frame(valor.real = dados.sim_1_2_1$w,
                        valor.media = apply(res_1_2_1$w, 2, mean),
                        grade = as.factor(rep(100,100)))

g2w_discr_2 <- data.frame(valor.real = dados.sim_1_2_2$w,
                          valor.media = apply(res_1_2_2$w, 2, mean),
                          grade = as.factor(rep(49,49)))

g2w_discr_3 <- data.frame(valor.real = dados.sim_1_2_3$w,
                          valor.media = apply(res_1_2_3$w, 2, mean),
                          grade = as.factor(rep(25,25)))

g2w_discr <- rbind(g2w_discr_1, g2w_discr_2, g2w_discr_3)

ggplot() +
  geom_point(aes(x = valor.real, y = valor.media, group = grade, shape = grade, color = grade), 
             size = 2.5, data = g2w_discr) +
  geom_abline(col = 1, size = 1) +
  xlab("Valor original") +
  ylab("Média") +
  labs(color = "", shape = "") +
  scale_shape_manual(values=c(17, 16, 15), labels = c("N = 100", "N = 49", "N = 25")) +
  scale_color_manual(values=c('red','blue', "darkgreen"), 
                     labels = c("N = 100", "N = 49", "N = 25")) +
  theme_bw() +
  theme(legend.position = "top",
        axis.title.x = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        axis.title.y = element_text(size = rel(1.5)),
        axis.text.y = element_text(size = rel(2)),
        axis.text.x = element_text(size = rel(2)))

# Parâmetros

DFplots_25 <- data.frame(iter = 1:1500,
                      beta = res_1_2_3$beta,
                      mu = res_1_2_3$mu,
                      tau = res_1_2_3$tau,
                      gamma = res_1_2_3$gamma)
DFplots_49 <- data.frame(iter = 1:1500,
                      beta = res_1_2_2$beta,
                      mu = res_1_2_2$mu,
                      tau = res_1_2_2$tau,
                      gamma = res_1_2_2$gamma)

beta_25 <- ggplot(aes(x = beta), data = DFplots_25) +
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

beta_49 <- ggplot(aes(x = beta), data = DFplots_49) +
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

mu_25 <- ggplot(aes(x = mu), data = DFplots_25) +
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
        
mu_49 <- ggplot(aes(x = mu), data = DFplots_49) +
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

tau_25 <- ggplot(aes(x = tau), data = DFplots_25) +
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

tau_49 <- ggplot(aes(x = tau), data = DFplots_49) +
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

gamma_25 <- ggplot(aes(x = gamma), data = DFplots_25) +
              geom_histogram(aes(y = ..density..), color = "black", fill = "grey", bins = 30) +
              geom_vline(xintercept = gamma.sim_1_2_3, linetype = "dashed", 
                         color = "red", size = 2) +
              geom_density(alpha = .2, fill="blue") +
              xlab("") +
              ylab("Densidade") +
              theme_bw() +
              theme(axis.title.x = element_text(size = rel(1.5)),
                    legend.text = element_text(size = rel(1.5)),
                    axis.title.y = element_text(size = rel(1.5)),
                    axis.text.y = element_text(size = rel(2)),
                    axis.text.x = element_text(size = rel(2)))

gamma_49 <- ggplot(aes(x = gamma), data = DFplots_49) +
              geom_histogram(aes(y = ..density..), color = "black", fill = "grey", bins = 30) +
              geom_vline(xintercept = gamma.sim_1_2_2, linetype = "dashed", 
                         color = "red", size = 2) +
              geom_density(alpha = .2, fill="blue") +
              xlab("") +
              ylab("Densidade") +
              theme_bw() +
              theme(axis.title.x = element_text(size = rel(1.5)),
                    legend.text = element_text(size = rel(1.5)),
                    axis.title.y = element_text(size = rel(1.5)),
                    axis.text.y = element_text(size = rel(2)),
                    axis.text.x = element_text(size = rel(2)))

grid.arrange(beta_25, beta_49, mu_25, mu_49, tau_25, tau_49, gamma_25, gamma_49, ncol = 2)



#save.image("Cenário 2 - Discretização/Comp.RData")
library(cli)

# Data generating process for IV misspecification study.
# Model:
#   U ~ N(0, 1)           instrument
#   Z ~ N(0, 1)           unobserved confounder (not in returned data)
#   X ~ N(a*U + b*Z, 1)   endogenous regressor
#   Y ~ N(c*X + d*Z, 1)   outcome

simulate_data <- function(n, a, b, c, d) {
  U <- rnorm(n)
  Z <- rnorm(n)
  X <- rnorm(n, mean = a * U + b * Z)
  Y <- rnorm(n, mean = c * X + d * Z)
  data.frame(U = U, X = X, Y = Y)
}

iv_estimator <- function(data) {
  cov(data$U, data$Y) / cov(data$U, data$X)
}


# Bayesian model (misspecified):
#   X ~ N(a*U + b*Z, 1),  Y ~ N(c*X + d*Z^p, 1)
#   Priors: a, b, c, d ~ N(0,1) independently; Z unobserved
#
# Conditionals for a, b, c, d are Gaussian (derived below).
# Z_i's conditional is non-Gaussian due to Z_i^p in Y, so we use an
# MH step: propose Z_i* from p(Z_i | X_i, a, b) = N(b*r_i/(1+b^2), 1/(1+b^2))
# where r_i = X_i - a*U_i, then accept/reject on the Y likelihood ratio.
gibbs_sampler <- function(data, n_iter = 5000, n_warmup = 1000, p=2) {
  U <- data$U
  X <- data$X
  Y <- data$Y
  n <- nrow(data)
  
  a <- 0; b <- 0; coef_c <- 0; d <- 0
  Z <- rnorm(n)
  
  c_samples <- numeric(n_iter)
  
  for (iter in seq_len(n_iter + n_warmup)) {
    # a | rest ~ N(sum((X - b*Z)*U) / (1 + sum(U^2)),  1/(1 + sum(U^2)))
    prec_a <- 1 + sum(U^2)
    a <- rnorm(1, sum((X - b * Z) * U) / prec_a, 1 / sqrt(prec_a))
    
    # b | rest ~ N(sum((X - a*U)*Z) / (1 + sum(Z^2)),  1/(1 + sum(Z^2)))
    prec_b <- 1 + sum(Z^2)
    b <- rnorm(1, sum((X - a * U) * Z) / prec_b, 1 / sqrt(prec_b))
    
    # c | rest ~ N(sum((Y - d*Z^p)*X) / (1 + sum(X^2)),  1/(1 + sum(X^2)))
    prec_c <- 1 + sum(X^2)
    coef_c <- rnorm(1, sum((Y - d * (Z*(1-p)+p*Z^2)) * X) / prec_c, 1 / sqrt(prec_c))
    
    # d | rest ~ N(sum((Y - c*X)*Z^p) / (1 + sum(Z^(2*p))),  1/(1 + sum(Z^(2*p))))
    prec_d <- 1 + sum((Z*(1-p)+p*Z^2)^2)
    d <- rnorm(1, sum((Y - coef_c * X) * (Z*(1-p)+p*Z^2)) / prec_d, 1 / sqrt(prec_d))
    
    # Z_i | rest: MH with proposal q(Z_i) = p(Z_i | X_i, a, b)
    #   acceptance log-ratio = log p(Y_i | Z_i*) - log p(Y_i | Z_i)
    prec_z <- 1 + b^2
    mu_z   <- b * (X - a * U) / prec_z
    Z_prop <- rnorm(n, mu_z, 1 / sqrt(prec_z))
    s <- Y - coef_c * X
    log_alpha <- -(s - d * (Z_prop*(1-p)+p*Z_prop^2))^2 / 2 + (s - d * (Z*(1-p)+p*Z^2))^2 / 2
    accept <- log(runif(n)) < log_alpha
    Z[accept] <- Z_prop[accept]
    
    if (iter > n_warmup) {
      c_samples[iter - n_warmup] <- coef_c
    }
  }
  
  mean(c_samples)
}


my_sims <- function(replicas, sample_sizes, n, n_iter, n_warm, seed = NULL){
  
  reps <- expand.grid(1:replicas, p, sample_sizes)
  
  df <- data.frame(coef_c = rep(NA, dim(reps)[1]), reps)
  names(df) <- c("coef_c", "n_chain", "p", "n")
  
  coef_c <- numeric()
  
  set.seed(seed)
  
  cli_progress_bar("Simulations done", total = dim(reps)[1])
  for(iter in 1:dim(reps)[1]){
    p <- reps[iter,2]
    n <- reps[iter,3]
    dat <- simulate_data(n, a = 1, b = 2, c = 3, d = 4)
    coef_c <- append(coef_c,
                     gibbs_sampler(dat[1:n,], p = p,
                                   n_iter = n_iter, n_warmup = n_warm))
    cli_progress_update()
  }
  cli_progress_done()
  df$coef_c <- coef_c
  
  df
}


################################################################################


n_iter = 5000
n_warm = 1000

replicas <- 200

sample_sizes <- c(100, 250, 500, 1000, 2500, 5000)

p <- c(0, .25, .5, .75, 1)

df <- my_sims(replicas, sample_sizes, n, n_iter, n_warm, seed = 42)

#saveRDS(df, file = "iv_miss.rds")


################################################################################

library(tidyverse)
library(ggplot2)

my_theme <- theme(
  legend.justification = c("right", "top"),
  legend.box.just = "right",
  legend.margin = margin(6, 6, 6, 6),
  legend.box.background = element_rect(colour = "black", size = .5),
  legend.text.align = 0,
  legend.title = element_text(size = 30, face = "bold"),
  legend.text = element_text(size = 25),
  axis.text = element_text(size = 25),
  axis.title = element_text(size = 35),
  strip.text = element_text(size = 35, face = "bold"),
  strip.background = element_blank()
)



df <- readRDS("iv_miss.rds")

df %>%
  filter(n >= 500) %>%
  ggplot(aes(x = factor(n), y = coef_c, fill = factor(p))) +
  geom_violin() + 
  geom_hline(yintercept = 3, linetype = "dashed") +
  ylab("Média a posteriori") +
  coord_transform(y = "log10") +
  xlab("Tamanho da amostra") +
  theme_bw() +
  scale_fill_discrete(name = "p =") +
  my_theme +
  theme(legend.position = "inside",
        legend.position.inside = c(0.99,0.98))

df %>%
  filter(n >= 500) %>%
  group_by(n, p) %>%
  summarise(mean(coef_c))


new_grid <- expand.grid(1:replicas, sample_sizes)

iv_est <- numeric(dim(new_grid)[1])

for(iter in 1:dim(new_grid)[1]){
  n <- new_grid[iter,2]
  dat <- simulate_data(n, a = 1, b = 2, c = 3, d = 4)
  iv_est[iter] <- iv_estimator(dat)
}

df2 <- data.frame(coef_c = iv_est, n_chain=1:replicas, p = -1, n = new_grid[,2])

df12 <- bind_rows(df, df2)

df12 <- df12 %>%
  group_by(n, p) %>%
  summarise(mean = mean(coef_c), min = min(coef_c), max = max(coef_c)) %>%
  pivot_longer(cols = c(mean, min, max)) %>%
  ungroup()  %>%
  pivot_wider(values_from = value, names_from = c(name, p))


xtable(df12)


df %>%
  group_by(n, p) %>%
  mutate(bin_c = coef_c >= 4.3) %>%
  summarise(bin_c = mean(bin_c)) %>%
  pivot_wider(values_from = bin_c, names_from = p)

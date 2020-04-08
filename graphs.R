library(mvtnorm)
library(tmvtnorm)
library(MultiRNG)
library(doSNOW)
library(purrr)
library(ggplot2)
library(patchwork)

sigma <- matrix(c(4, 0, 0, 4), nrow=2)
mu <- c(2, 2)
mu_star <- c(2, 2)
n <- 1000
x <- rmvnorm(n, mean=mu, sigma=sigma)
x_star_n <- rmvnorm(100, mean=mu_star, sigma=sigma)
x_star_t <- rtmvnorm(100, mean=mu_star, sigma=sigma)
x_star_u <- draw.d.variate.uniform(100,2,sigma)

T_n <- ((x_star_n[, 1] - colMeans(x)[1])/sqrt(var(x[, 1]))) - ((x_star_n[, 2] - colMeans(x)[2])/sqrt(var(x[, 2])))
T_t <- ((x_star_t[, 1] - colMeans(x)[1])/sqrt(var(x[, 1]))) - ((x_star_t[, 2] - colMeans(x)[2])/sqrt(var(x[, 2])))
T_u <- ((x_star_u[, 1] - colMeans(x)[1])/sqrt(var(x[, 1]))) - ((x_star_u[, 2] - colMeans(x)[2])/sqrt(var(x[, 2])))


df <- data.frame(T = T_n)
pn <- ggplot(df, aes(x = T)) + 
  geom_histogram(aes(y =..density..),
                 breaks = seq(-4, 4, by = 1), 
                 colour = "black") +
  stat_function(fun = dnorm, args = list(mean = mean(df$T), sd = sd(df$T))) + 
  theme_minimal() + 
  ggtitle("T distribution")


qq_n <- data.frame(t=T_n) %>% 
  ggplot(aes(sample=t)) +
  stat_qq() +
  theme_minimal() +
  ggtitle("QQ Plot")


pn | qq_n
ggsave("1_normal.jpeg", width=6, height=4.5, unit="in")



df <- data.frame(T = T_t)
pt <- ggplot(df, aes(x = T)) + 
  geom_histogram(aes(y =..density..),
                 breaks = seq(-4, 4, by = 1), 
                 colour = "black") +
  stat_function(fun = dnorm, args = list(mean = mean(df$T), sd = sd(df$T))) + 
  theme_minimal() + 
  ggtitle("T distribution")

qq_t <- data.frame(t=T_t) %>% 
  ggplot(aes(sample=t)) +
  stat_qq() +
  theme_minimal() +
  ggtitle("QQ Plot")

pt | qq_t
ggsave("2_tdist.jpeg", width=6, height=4.5, unit="in")


df <- data.frame(T = T_u)
pu <- ggplot(df, aes(x = T)) + 
  geom_histogram(aes(y =..density..),
                 breaks = seq(-1, 1, by = 0.2), 
                 colour = "black") +
  stat_function(fun = dnorm, args = list(mean = mean(df$T), sd = sd(df$T))) + 
  theme_minimal() + 
  ggtitle("T distribution")

qq_u <- data.frame(t=T_u) %>% 
  ggplot(aes(sample=t)) +
  stat_qq() +
  theme_minimal() +
  ggtitle("QQ Plot")

pu | qq_u
ggsave("3_uniform.jpeg", width=6, height=4.5, unit="in")


shapiro.test(T_n)
shapiro.test(T_t)
shapiro.test(T_u)

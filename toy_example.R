# This is the code for the example we use in the book chapter
# Copyright Richard Torkar, Robert Feldt, Carlo A. Furia
#
# We will use both Rethinking and brms so that people can see both approaches
library(rethinking)
library(brms)
library(tidyverse)
library(ggthemes)
library(extraDistr)
library(ggpubr)

data(Howell1)
d <- Howell1
d2 <- d[ d$age >= 18 , ]

p_mu <- ggplot(data.frame(x = c(100, 250)), aes(x = x)) +
  stat_function(fun = dnorm, args=list(mean=181, sd = 20)) + theme_tufte() +
  xlab(expression(alpha)) +
  ylab("") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

p_beta <- ggplot(data.frame(x = c(-50, 50)), aes(x = x)) +
  stat_function(fun = dnorm, args=list(mean=0, sd = 10)) + theme_tufte() +
  xlab(expression(beta)) +
  ylab("") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

p_cauchy <- ggplot(data.frame(x = c(0, 100)), aes(x = x)) +
  stat_function(fun = dhcauchy, args=list(sigma=10)) + theme_tufte() +
  xlab(expression(sigma)) +
  ylab("") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

sample_mu <- rnorm(1e4, 181, 20)
sample_sigma <- rhcauchy(1e4, 10)
sample_beta <- rnorm(1e4, 0, 10)
prior_h <- rnorm(1e4, sample_mu + sample_beta, sample_sigma)
prior_h <- data.frame(length=prior_h)

p_combo <- ggplot(data.frame(prior_h), aes(x=length)) +  geom_density() + xlim(0,360) + 
  xlab(expression(mu, " and ", sigma)) +
  ylab("") +
  theme_tufte() 

ggarrange(ggarrange(p_mu, p_beta, p_cauchy, ncol=3),
          p_combo, nrow=2)

model <- ulam(
  alist(
    height ~ normal(mu, sigma),
    mu <- alpha + beta_w * weight,
    alpha ~ normal(180, 20),
    beta_w ~ normal(0,10),
    sigma ~ cauchy(0,10)
  ), data = d2, chains = 4, cores = 4, iter = 1e4
)

precis(model)
posterior <- extract.samples(model)
weight.seq <- seq(from=25, to=70, length.out=352)
sim <- sim(model, data=list(weight=weight.seq))
height.PI <- apply( sim , 2 , PI , prob=0.95 )

a <- mean(posterior$alpha)
b <- mean(posterior$beta_w)

ggplot(foo, aes(x=weight, y=height)) + geom_point(color="grey90") + 
  geom_abline(intercept=a, slope=b) + 
  theme_tufte()

# Or why not in brms (height is model against 1 (which is alpha) and the predictor weight)
model <- brm(bf(height ~ 1 + weight), 
             prior = c(prior(normal(180,20), class = "Intercept"),
                               prior(normal(0,10), class="b"),
                               prior(cauchy(0,10), class="sigma")),
             data=d2, chains=4, cores=4, iter=1e4)

weight_seq <- tibble(weight = seq(from = 25, to = 70, by = 1))

mu <-
  fitted(model,
         summary = F,
         newdata = weight_seq) %>%
  as_tibble() %>%
  mutate(Iter = 1:20000) %>%
  select(Iter, everything())


pred_height <-
  predict(model,
          newdata = weight_seq) %>%
  as_tibble() %>%
  bind_cols(weight_seq)

mu_summary <-
  fitted(model, 
         newdata = weight_seq) %>%
  as_tibble() %>%
  bind_cols(weight_seq)

d2 %>%
  ggplot(aes(x = weight)) +
  geom_ribbon(data = pred_height, 
              aes(y = Estimate, ymin = Q2.5, ymax = Q97.5),
              fill = "grey85") +
  geom_ribbon(data = mu_summary, 
              aes(y = Estimate, ymin = Q2.5, ymax = Q97.5),
              fill = "grey60") +
  geom_line(data = mu_summary, aes(y = Estimate)) +
  geom_point(aes(y = height), shape = 1, size = 1.5, alpha = 2/5) +
  coord_cartesian(xlim = range(d2$weight),
                  ylim = range(d2$height)) +
  theme(panel.grid = element_blank()) + ylab("height") +
  theme_tufte()

pred <- predict(model)

## Species-area relationships

source("03_setup.R")

# 1. Alpha diversity ------------------------------------------------------

# Let's create a Fascia variable (100m height altitudinal belts), that we
# have for gamma diversity, for alpha diversity: 
speciesrich$Fascia <- as.factor(trunc(speciesrich$quota/100)*100)

# Let's see how the species-area relationship is structured by the altitudinal
# belts:
ggplot(speciesrich, aes(x = as.numeric(size), y = alpha, color = Fascia)) +
  geom_point() +
  geom_smooth(se = F, method = "lm", formula = "y ~ I(log1p(x))") +
  labs(x = "Plot size", y = "Alpha diversity") +
  theme_classic() +
  scale_color_viridis_d()

# It seems there is a strong effect of the Fascia on the slope of the SA
# relationship, generally diminishing with altitude. Let's see how does this
# decay of the slope of the SAR relates with altitude.

# 1. Make a different log-linear model for each Fascia value, and extract
# the coefficients and confidence intervals.
m_df_SA_alpha <- 
  tapply(speciesrich, 
         speciesrich$Fascia, 
         function(df){
           m <- lm(alpha ~ I(log1p(size2)), data = df)
           c <- coef(m)[2]
           ci <- confint(m)[2,]
           d  <- data.frame("size" = c, "lwr_ci" = ci[1], "upr_ci" = ci[2])
           return(d)
         }
  ) %>% 
  Reduce(rbind, .) %>% 
  data.frame(., fascia = unique(speciesrich$Fascia))

# The R2 of the different models (easily obtained modifying the previous
# tapply) are between 75 and 99%, without trend along the fascia.

# Is any slope non-significant?
all(sign(m_df_SA_alpha$lwr_ci) == sign(m_df_SA_alpha$upr_ci)) # TRUE
# This means that the lower and upper confidence intervals never cross the 0 line,
# so all the slopes are statistically significant.

# Let's explore three possible types of decay of the slope with altitude:
forms <- 
  alist("linear_decay" = size ~ as.numeric(fascia),
        "log_decay" = size ~ I(log(as.numeric(fascia))),
        "geom_decay" = size ~ I(1/(as.numeric(fascia)))
  )

# Extract Akaike weights to assess the validity of each type.
lapply(forms, lm, data = m_df_SA_alpha) %>% sapply(AIC) %>% phytools::aic.w()
# The logarithmic decay has 92% of the Akaike weights.

# We can extract the LaTeX output of the relationship with:
equatiomatic::extract_eq(lm(forms$log_decay, data = m_df_SA_alpha), use_coefs = T)

# Let's take it as the valid one, and plot it with:
m_df_SA_alpha %>% 
  ggplot(aes(x = fascia, y = size, group = 1)) +
  geom_point() +
  geom_line() +
  geom_smooth(method = "lm", formula = "y ~ log(x)") +
  theme_classic()

# 2. Gamma diversity ------------------------------------------------------

# Again, explore how the SAR is structured by altitudinal belts:
ggplot(gamma_fascia, aes(x = as.numeric(size), y = gamma, color = as.factor(Fascia))) +
  geom_point() +
  geom_line() +
  geom_smooth(se = F, method = "lm", formula = "y ~ I(log1p(x))") +
  labs(x = "Plot size", y = "Gamma diversity", color = "Altitudinal\nbelt") +
  scale_color_viridis_d() +
  # guides(color = )
  theme_classic()

# The structure by altitude is even stronger. Let's replicate the previous analyses
# for gamma diversity:
m_df_SA_gamma <- 
  tapply(gamma_fascia, 
         gamma_fascia$Fascia, 
         function(df){
           m <- lm(gamma ~ I(log1p(size2)), data = df)
           c <- coef(m)[2]
           ci <- confint(m)[2,]
           d  <- data.frame("size" = c, "lwr_ci" = ci[1], "upr_ci" = ci[2])
           return(d)
         }
    ) %>% 
    Reduce(rbind, .) %>% 
    data.frame(., fascia = unique(gamma_fascia$Fascia))

# R2 of the models is between 94 and 99%, without trend along the gradient.

# And again validate the different types of decay
lapply(forms, lm, data = m_df_SA_gamma) %>% sapply(AIC) %>% phytools::aic.w()
# The geometric decay has 90% of the Akaike weight. 

# Extract the LaTeX formula with:
equatiomatic::extract_eq(lm(forms$log_decay, data = m_df_SA_gamma), use_coefs = T)

# Let's represent it with:
m_df_SA_gamma %>% 
  ggplot(aes(x = fascia, y = size, group = 1)) +
  geom_pointrange(aes(ymin = lwr_ci, ymax = upr_ci)) +
  geom_line() +
  geom_smooth(method = "lm", formula = "y ~ I(1/x)") +
  theme_classic()


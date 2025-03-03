
source("03_setup.R")
# Visualization customs for a geom_violin type of plot
#source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")


# Visualization of species richness along the elevational gradient structured
# by plot area:
speciesrich  %>% 
  mutate(size = as.numeric(size)) %>%
  arrange(size) %>% 
  # filter(size != 0.015, size != 0.03125) %>% 
  mutate(size = as.factor(size)) %>% 
  ggplot(aes (x = quota, y = alpha, col = size, group = size)) +
  # geom_point() +
  geom_smooth(formula =  y ~ poly(x,1), method = "lm", size = 1.5, alpha = 0.20) +
  # scale_y_discrete(limits = factor(0:80))+
  theme_bw() +
  xlab("Elevation") +
  ylab("Alpha diversity") +
  scale_color_brewer(palette = "Set1") 

# We can clearly see that the grain size (or area) affects the correlation 
# between alpha-diversity and elevation: fine grain size is not affected by 
# elevation, but the negative correlaiton between elevation and alpha diversity
# becomes stronger on increasing grain size. We can explore what is the "scale
# of effect"  

# We had divided our species richness by grain size, which will come in 
# handy to run separate analyses. The object in question is sr_split.



# 1. Pooled data model ----------------------------------------------------

# Let's first fit a general model of the relationship we have seen on the 
# initial plot: alpha diversity explained by an interaction of plot size
# and altitude:
m <- glm.nb(alpha ~ quota + I(log(size2, base=2)), data = speciesrich) 
# library(glmmTMB)
# mx <- glmmTMB::glmmTMB(alpha ~ offset(log(size2, base=2)) + quota, data = speciesrich, family = nbinom1) 
# my <- glmmTMB::glmmTMB(alpha ~ offset(log(size2, base=2)), data = speciesrich, family = nbinom1) 
# 1 - deviance(mx)/deviance(my) 

# Let's first extract some important information from the model: fitted values,
# residuals, etc:
norm.res <- broom::augment(m)

# And let's check the R2:
1 - m$deviance/m$null.deviance
# 0.85, very good!
AIC(m)

# First let's try to be sure that there is no problem of spatial autocorrelation
# and of normality of residuals. 
# We have to obtain a geometric distance matrix 
D <- st_distance(coords)
# And get rid of it's annoying "units" and metadata
class(D) <- "matrix"

# The following functions come from the DHARMa package
m %>% 
  # Obtain simulated standardized residuals
  simulateResiduals(n = 1e3) %>%
  # Each observation across a size range comes from the same "location", 
  # so there are ties in the geographic coordinates. We have to account for these
  # with the following line:
  recalculateResiduals(group = rep(1:83, each = 9)) %>%
  # And now we can test spatial autocorrelation of residuals through Moran's I:
  testSpatialAutocorrelation(distMat = D)

# Moran's I is extremely small and p-value around 0.5 so nothing to worry about

# Now let's see the residuals structure: add simulated residuals to the
# database of residuals we already had:
norm.res$std_res <- 
  residuals(
    DHARMa::simulateResiduals(m, n = 1e3, integerResponse = F, refit = F, plot = F),
    quantileFunction = qnorm
  )

# Now create a qqplot of standardized vs theoretical residuals:
ggplot(data = norm.res, mapping = aes(sample = std_res)) + 
  qqplotr::stat_qq_point()+
  theme_classic()+
  qqplotr::stat_qq_line()+
  qqplotr::stat_qq_band(alpha=0.3)

# And check whether the residuals are flat and have constant variance 
# against the fitted values
ggplot(data=norm.res, mapping = aes(x = .fitted, 
                                    y = std_res)) + 
  geom_point(col="black", size=2)+
  geom_hline(yintercept =0,linetype = 2, size=1.2)+
  theme_classic()+ 
  labs(x="Fitted",y="Residuals")+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

# The assumptions of the linear model are generally well supported. 
# Nothing to worry about!

# Finally let's explore the variance explained by each of variables in our
# model: quota and size by themselves, and their interaction. 
# First we have to create a function that returns the deviance explained
# from the glm.nb model:
NB_glm_r2 <- 
  function(form, df = speciesrich){
  m <- glm.nb(form, data = df)
  1 - m$deviance/m$null.deviance
}

# And feed it to the domir function in the domir package, that requires a 
# formula and the function that extracts the R2 from it:
domir::domir(alpha ~ quota + log(size2), NB_glm_r2)
domir::domir(alpha ~ quota * log(size2), NB_glm_r2)

# The R2 basically doesnt change between the two models. However, when we 
# add the interaction, we see that its term captures a large part of the variance
# initially explained by size alone. Let's compare models varying both the
# interaction and the polynomial degree of the altitude:

m  <- glm.nb(alpha ~ quota + I(log(size2, 2)), data = speciesrich)
m2 <- glm.nb(alpha ~ poly(quota, 2) + I(log(size2, base = 2L)), data = speciesrich)
m3 <- glm.nb(alpha ~ poly(quota, 3) + I(log(size2, base = 2L)), data = speciesrich)

m_int  <- glm.nb(alpha ~ quota * I(log(size2, 2)), data = speciesrich)
m_int2 <- glm.nb(alpha ~ poly(quota, 2) * I(log(size2, base = 2L)), data = speciesrich)
m_int3 <- glm.nb(alpha ~ poly(quota, 3) * I(log(size2, base = 2L)), data = speciesrich)

# Assess their relative evidence with AIC
AIC(m, m2, m3, m_int, m_int2, m_int3) %>% .[, 2] %>% phytools::aic.w()

# Apparently only those of 3rd degree have substantial support. Let's test
# if the interaction is statistically significant:
anova(m3, m_int3)

# No, so for prediction purposes we don't need the interaction term 
# between quota and plot size. However, we saw by the domir function
# That there is a large share of variation explained by the interaction
# of the two variables. Hence, let's make a plot of the variance explained
# by each variable independently and of shared variance.


# First, however, let's just make a visualization of the model fit 
# against the data:
fit_vals_m3 <- data.frame("quota" = speciesrich$quota, 
                          "fitted" = fitted(m_int3),
                          "alpha" = speciesrich$alpha,
                          "size" = as.factor(speciesrich$size2))
fit_vals_m3 %>% 
  ggplot(aes(x = quota, y = fitted, color = size)) +
  geom_point(aes(y = alpha), alpha = .5) +
  geom_line(linewidth = 1.2) +
  theme_classic() +
  scale_color_viridis_d()

# Store the different polynomial terms of quota within speciesrich:
# speciesrich$polyquota1 <- poly(speciesrich$quota, 1)
# speciesrich$polyquota2 <- poly(speciesrich$quota, 2)[,2]
# speciesrich$polyquota3 <- poly(speciesrich$quota, 3)[,3]

# Let's obtain the decomposition of explained variance by each term
# using dominance analysis. We must use the manual exponentiation
# instead of the poly() function as in the previous lines because the
# algorithm doesn't converge:
r2 <- domir::domir(alpha ~ 
                     I(quota**1) * log(size2) +
                     I(quota**2) * log(size2) +
                     I(quota**3) * log(size2), 
                   NB_glm_r2)

# Let's turn the results it in a suitable dataframe
r2_df <-
  round(r2$General_Dominance, 2) %>% 
  cbind(., c("Linear", "Log-size", "Quadratic", "Cubic", "interaction1", "interaction2", "interaction3")) %>%
  as.data.frame()
colnames(r2_df) <- c("var_explained", "predictor")
r2_df$var_explained <- as.numeric(r2_df$var_explained)

# Aggregate the variance explained by all the interactions of quota and log-size
r2_df <- rbind(r2_df[1:4, ], 
               data.frame(
                 var_explained = sum(r2_df$var_explained[5:7]), 
                 predictor =  "Interaction")
)
# Give the adequate order to the factor (important for visualization)
r2_df$predictor <- factor(r2_df$predictor, c("Interaction", "Log-size", "Linear", "Quadratic", "Cubic"))

# Visualize:
ggplot(r2_df,
       aes(x = 1, y = var_explained,
           fill = predictor)) +
  geom_col(position = "dodge") +
  labs(x = "Predictor", y = "Variance explained") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20, face = "bold"),
    axis.text.y = element_text(size = 18)
  ) +
  labs(fill = "Predictor variable") +
  scale_fill_viridis_d()

# It is obvious that the accumulated interactions explain most of the variance.
# PENDING: CHECK THE .fct ARGUMENT IN domin()

# Now we know how the species diversity reacts as a whole, explainable as:
summary(m_int3)
saveRDS(summary(m_int3), "./results/RDS/summary_m_int3.RDS")
saveRDS(summary(m3), "./results/RDS/summary_m3.RDS")

# The interaction terms are statistically insignificant, but they point anyway
# to a negative relationship between plot size and the effect of the linear term,
# and a positive relationship between plot size and the effect of higher-order
# terms. Let's explore this in more detail and see whether the relationship is
# actually linear and how strong


# 2. Plot size-divided models ---------------------------------------------

# First let's define a function to extract the coefficients
extract_coefs <- 
  function(mod, var){
    eff <- coef(mod)[[var]]
    l_ci <- confint(mod)[var, 1]
    u_ci <- confint(mod)[var, 2]
    df <- data.frame(eff, l_ci, u_ci)
    rownames(df) <- names(coef(mod))[[var]]
    df
  }

# And apply this for our previous model:
sapply(1:6, extract_coefs, mod = m_int3) %>% 
  t() %>% 
  set_rownames(., c("Intercept", "linear_quota", "quad_quota", "cube_quota", "log_size", "interaction_quota.size")) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "predictor")

# It works, returning us our effect sizes and confidence intervals.
# Now, if we make individual models for each plot size value, we can extract 
# this information easily

sizes <- unique(speciesrich$size2)

# Store the summaries of the models:
summary_lists <- 
  lapply(sizes,
         function(x){
           df <- speciesrich[speciesrich$size2 == x, ]
           m <- glm.nb(alpha ~ poly(quota, 3), data = df)
           summary(m)
           }
  )
saveRDS(summary_lists, "./results/RDS/summary_lists.RDS")

# Extract the effect sizes with our function:
df_effect_sizes <-
  lapply(sizes,
       function(x){
         df <- speciesrich[speciesrich$size2 == x, ]
         m <- glm.nb(alpha ~ poly(quota, 3), data = df)
         res <- t(sapply(1:4, extract_coefs, mod = m))
         res <- apply(res, 2, unlist)
         res <- as.data.frame(res)
         res$var  <- factor(c("Intercept", "Linear", "Quadratic", "Cubic"), 
                            levels = c("Intercept", "Linear", "Quadratic", "Cubic")
                            )
         res$size <- x
         return(res)
       }) %>% 
  Reduce(rbind, .)

# Visualize how do the effect sizes for each term change by plot-size:
df_effect_sizes %>% 
  as_tibble() %>% 
  ggplot(aes(x = log(size,2L), y = eff,
             alpha = ifelse(sign(l_ci) == sign(u_ci), T, F))) +
  geom_pointrange(aes(ymin = l_ci, ymax = u_ci), linewidth = 1) +
  geom_line(group = 1) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_classic() +
  labs(x = "Plot size (log2 scale)", y = "Effect size", alpha = "Significant",
       title = "Effects of altitude on\nalpha diversity across plot sizes") +
  scale_color_viridis_d() +
  theme(
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 18),
    plot.title = element_text(size = 20, face = "bold", hjust = .5),
    strip.text = element_text(size = 18, face = "bold"),
    strip.background = element_rect(fill = "oldlace")
  ) +
  facet_wrap(~var, scales = "free_y")

# We see what the signs in summary(m_int3) were telling us: there is an increasingly
# more negative effect of quota on the linear term across plot sizes, and 
# more positive effects of all the other terms. However, these are, as the
# summary told us, pretty weak.

# Hence, it is safe to assume that the m3 model, with no interaction between
# plot size and quota is a good, parsimonious, representation of our data:
# The effects of quota are quite consistent across plot sizes, with very 
# weak trends across their different linear, quadratic and cubic components.

# The striking pattern seen in our first exploratory visualization, that 
# there is a downward trend across elevation that is much stronger for
# large plot sizes, is very likely due to the fact that the negative binomial
# is a parametric distribution: the higher the mean, the higher the variance,
# so the 

d2_AIC_alpha <-   
  map(
    sr_split, 
    function(x){
    # x <- sr_split[[9]]
    # deh <- list()
    
      # Let's define a battery of models with increasing complexity:
      # Either linear, quadratic or cubic
      glm_1alpha <- glm (formula = alpha ~ poly(quota, 1), data = x, family = poisson)
      glm_2alpha <- glm (formula = alpha ~ poly(quota, 2), data = x, family = poisson)
      glm_3alpha <- glm (formula = alpha ~ poly(quota, 3), data = x, family = poisson)
      
      # Select the model by AIC differences
      if(summary(glm_1alpha)$aic - 2 <= summary(glm_2alpha)$aic){
        glm_alpha = glm_1alpha
      } else if(summary(glm_2alpha)$aic - 2 <= summary(glm_3alpha)$aic){
          glm_alpha = glm_2alpha
      }else{glm_alpha = glm_3alpha}
      
      # Obtain summary statistics from the fitted model
      sum_glm <- summary(glm_alpha)
      coef_model <- sum_glm$coefficients 
      names_coef <- 
        c("Int", "p_interc", "Elev1", "p1", "Elev2", "p2", "Elev3", "p3") %>%
        .[1:(nrow(coef_model)*2)]
      
      coef_model <- 
        coef_model %>%
        t() %>%
        as.data.frame() %>%  
        rowwise() %>% 
        map(., function(y) y[c(1,4)]) %>% 
        do.call(c,.) %>% 
        setNames(., names_coef)
      
      # deh <- list()
      
      # Obtain an evaluation from the model based on the deviance squared
      # and a summary of model coefficients
      deh <- 
        c(D2 = modEvA::Dsquared(glm_alpha), AIC = sum_glm$aic, coef_model) %>%
        t() %>%
        as.data.frame() %>%
        cbind(size = x$size %>% unique())
      
      return(deh) 
    
    }
  )
  
table <-  
  d2_AIC_alpha %>% 
  do.call(bind_rows, .) %>% 
  mutate(size = grain_sizes2) %>% 
  arrange(size)

#coef_alpha <- read.csv("coef_alpha.csv")

#write.csv(table, "coef_alpha.csv")
    
     AIC_alpha <-  table %>% mutate(across(where(is.numeric), .fns =  ~ round(.x, 3))) %>%
      mutate(p1 = case_when(
           p1 > 0.05 ~ "ns",
           p1 > 0.01 & p1 <= 0.05 ~ "*",
           p1 > 0.001  & p1 <= 0.01 ~ "**",
           p1 <= 0.001 ~"***" )) %>% 
       
         mutate(p2 = case_when(
           p2 > 0.05 ~ "ns",
           p2 > 0.01 & p2 <= 0.05 ~ "*",
           p2 > 0.001  & p2 <= 0.01 ~ "**",
           p2 <= 0.001 ~"***" )) %>%
      
       
       mutate(p3 = case_when(
         p3 > 0.05 ~ "ns",
         p3 > 0.01 & p3 <= 0.05 ~ "*",
         p3 > 0.001  & p3 <= 0.01 ~ "**",
         p3 <= 0.001 ~"***" )) %>%
       
       
       mutate(p_interc = case_when(
         p_interc > 0.05 ~ "ns",
         p_interc > 0.01 & p_interc <= 0.05 ~ "*",
         p_interc > 0.001  & p_interc <= 0.01 ~ "**",
         p_interc <= 0.001 ~"***" )) %>%
       replace(is.na(.), "") %>% 
       
       gt(rowname_col = "size") %>% 
         
        tab_header (title = md(paste("**Alpha diversity Vs Elevation**", sep = ""))) %>% 
         
         tab_style(style = cell_text(style = "italic"),
                   location = cells_stub()) %>% 
         cols_label(
           D2 = "D2", 
           p1 = "",
           p2 = "", 
           p3 = "",
           p_interc = "",
           AIC = "AIC") %>% 
       show()
       #gtsave(filename = "GLMs Alpha.png")
     
  
     
     predict_alpha <-   map(sr_split, function(x){
     #x = sr_split[[1]]
       
       glm_1alpha <- glm (formula = alpha ~ poly(quota, 1), data = x, family = poisson)
       glm_2alpha <- glm (formula = alpha ~ poly(quota, 2), data = x,family = poisson)
       glm_3alpha <- glm (formula = alpha ~ poly(quota, 3), data = x, family = poisson)
       
       if(summary(glm_1alpha)$aic - 2 <= summary(glm_2alpha)$aic){
         glm_alpha = glm_1alpha
       } else if(summary(glm_2alpha)$aic - 2 <= summary(glm_3alpha)$aic){
         glm_alpha = glm_2alpha
       }else{glm_alpha = glm_3alpha}
       
       #newdata <- list()
       
       newdata<- seq(min(x$quota), max(x$quota), length = 100) %>% 
         cbind.data.frame(quota = ., size = rep(x$size %>% unique(), 100))
      
       newdata$alpha <- predict(glm_alpha, newdata = newdata, type = "response")
       
       return(newdata)
       
     }) 
     
     
    predict <-  predict_alpha %>% do.call(bind_rows, .) %>% 
       mutate(size = as.numeric(size)) %>%
       arrange(size) %>% 
       filter(size != 0.015, size != 0.03125) %>% 
       mutate(size = as.factor(size))
    
    #write.csv(predict, "predict_alpha.csv")
     
    #predict <- read.csv("predict_alpha.csv")
    
     mypalette <- c("#E0DA34", "#99B81E", "#37D46B", "#3AE0C7", "#1F65B5", "#8034D1", "#E00DB6", 
                    "#E6094F", "#DB881B")
sr_plot <-   speciesrich %>% 
       mutate(size = as.numeric(size)) %>%
       arrange(size) %>% 
       filter(size != 0.015, size != 0.03125) %>%
       mutate(size = as.factor(size)) %>% 
       ggplot(aes (x = quota, y = alpha, col = size, group = size)) +
       geom_point(alpha = 0.2) +
  guides(color = guide_legend(title = "Sample size")) +
       geom_line( data = predict, mapping = aes(x = quota, col = size, group = size, y = alpha), size = .9, alpha = 0.60) +
       theme_bw()+
       labs(x = "Elevation m a.s.l.", y = "Alpha diversity") +
       scale_colour_brewer(palette = "Paired") 

    # ggsave(filename = "alphadive.tiff", plot =  sr_plot, 
     #              width = 180, height = 160, units = "mm", dpi = 300, compression="lzw")
     
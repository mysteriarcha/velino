
library(sf) 
library(tidyverse)
library(vegan)
library(gt)
library(modEvA)


######## D2 e AIC


d2_AIC_gamma <-   map(gam_split, function(x){
 # x = gam_split[[3]]
  
  
  glm_1gamma <- glm (formula = gamma ~ poly(Fascia, 1), data = x, family = poisson ("identity"))
  glm_2gamma <- glm (formula = gamma ~ poly(Fascia, 2), data = x, family = poisson ("identity"))
  glm_3gamma <- glm (formula = gamma ~ poly(Fascia, 3), data = x, family = poisson ("identity"))
  
  if(summary(glm_1gamma)$aic - 2 <= summary(glm_2gamma)$aic){
    glm_gamma = glm_1gamma
  } else if(summary(glm_2gamma)$aic - 2 <= summary(glm_3gamma)$aic){
    glm_gamma = glm_2gamma
  }else{glm_gamma = glm_3gamma}
  
  sum_glm <- summary(glm_gamma)
  
  coef_model <- sum_glm$coefficients 
  
  nomi_coef <- c("Int", "p_interc", "Elev1", "p1", "Elev2", "p2", "Elev3", "p3") %>%
    .[1:(nrow(coef_model)*2)]
  
  coef_model <- coef_model %>%
    t() %>% as.data.frame() %>%  
    rowwise() %>% 
    map(., function(y) y[c(1,4)]) %>% 
    do.call(c,.) %>% 
    setNames(., nomi_coef)
 
  
  deh <- c(D2 = modEvA::Dsquared(glm_gamma), AIC = sum_glm$aic, coef_model) %>%
    t() %>%
    as.data.frame() %>%
    cbind(size = x$size %>% unique())
  
  return(deh) 
}) %>% do.call(bind_rows, .)  %>%
  arrange (factor(size, levels = c ('0.0625', '0.125', '0.25', 
                                    '0.5', '1', '2', '4', '8', '16')))  

AIC_gamma <-  d2_AIC_gamma %>% mutate(across(where(is.numeric), .fns =  ~ round(.x, 3))) %>%
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
  
    tab_header (title = md(paste("**Gamma diversity Vs Elevational belts**", sep = ""))) %>% 
  
  tab_style(style = cell_text(style = "italic"),
            location = cells_stub()) %>% 
  cols_label(
    D2 = "D2", 
    p1 = "",
    p2 = "", 
    p3 = "",
    p_interc = "", 
    AIC = "AIC") %>% 
  cols_hide(Elev3) %>% 
  show()#%>% 

#gtsave("gamma Vs Belts.png")

predict_gamma <-   map(gam_split, function(x){
   #x = gam_split[[1]]
  
  glm_1gamma <- glm (formula = gamma ~ poly(Fascia, 1), data = x, family = poisson ("identity"))
  glm_2gamma <- glm (formula = gamma ~ poly(Fascia, 2), data = x, family = poisson ("identity"))
  glm_3gamma <- glm (formula = gamma ~ poly(Fascia, 3), data = x, family = poisson ("identity"))
  
  if(summary(glm_1gamma)$aic - 2 <= summary(glm_2gamma)$aic){
    glm_gamma = glm_1gamma
  } else if(summary(glm_2gamma)$aic - 2 <= summary(glm_3gamma)$aic){
    glm_gamma = glm_2gamma
  }else{glm_gamma = glm_3gamma}
  
  newdata <- seq(min(x$Fascia), max(x$Fascia), length = 100) %>% 
    cbind.data.frame(Fascia = ., size = rep(x$size %>% unique(), 100))
  
  newdata$gamma <- predict(glm_gamma, newdata = newdata, type = "response")
  
  return(newdata)
  
}) %>% do.call(bind_rows, .) %>% 
  mutate(size = as.numeric(size)) %>%
  arrange(size) %>% 
  filter(size != 0.015, size != 0.03125) %>% 
  mutate(size = as.factor(size))

 gamma_plot <- gamma_fascia %>% 
  mutate(size = as.numeric(size)) %>%
  arrange(size) %>% 
  filter(size != 0.015, size != 0.03125) %>%
  mutate(size = as.factor(size)) %>% 
  ggplot(aes (x = Fascia, y = gamma, col = size, group = size)) +
  geom_point(alpha = 0.2) +
  guides(color = guide_legend(title = "Sample size")) +
  geom_line( data = predict_gamma, mapping = aes(x = Fascia, col = size, group = size, y = gamma),
             size = .9, alpha = 0.60) +
  theme_bw()+
  labs(x = "Elevation m a.s.l.", y = "Gamma diversity") +
  scale_colour_brewer(palette = "Paired") 
ggsave(filename = "gammadive.tiff", plot =  gamma_plot, 
       width = 180, height = 160, units = "mm", dpi = 300, compression="lzw")

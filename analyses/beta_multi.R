

library(tidyverse)
library(betapart)
library(sf)
library(purrr)
library(gt)
library(betareg)
library(AICcmodavg)



####Data

ril_vel <- read.csv("velino_firstcycle.csv", row.names = 1) %>% 
  t() %>% as.data.frame() %>% 
  add_column(ril = word(rownames(.), 1, sep = '_'), 
             size = word(rownames(.), 2, sep = '_') ) 

coords <- st_read("coords_rilievi.shp")

vel_env <- read.csv ("velino_environment.csv", row.names = 1, sep = ";") %>% 
  cbind(coords$geometry)




##calcolo beta.multi per fascia

bm_ril_fascia <- vel_env %>% 
  rownames_to_column("ril") %>% 
  select(ril, Fascia) %>% 
  left_join(ril_vel,.) %>% 
  dplyr::select(-ril) %>% 
  group_by(Fascia, size) %>% 
  group_split() %>% 
  map_dfr(~ beta.multi(.x %>% dplyr::select(-Fascia, -size)) %>% 
            c(., fascia = .x$Fascia %>% unique(), size = .x$size %>% unique())) 


pivot_bm <- bm_ril_fascia %>% 
  pivot_longer(cols = starts_with("beta"), names_to = "beta", values_to = "value")



####beta.multi Vs quota 

group_sizes <- bm_ril_fascia %>% 
  filter(size != 0.015, size != 0.03125) %>% 
  group_by(size) %>% 
  group_split()  
 


d2_AIC_SOR <-   map(group_sizes, function(x){
#x <- group_sizes[[1]]
 
 #print(x$size %>% unique())
   
 x <- x %>% 
    drop_na() 
  

   
  glm_1SOR <- MASS::glm.nb(formula = beta.SOR ~ poly(fascia, 1), data = x) 
  glm_2SOR <- MASS::glm.nb(formula = beta.SOR ~ poly(fascia, 2), data = x)
  glm_3SOR <- MASS::glm.nb(formula = beta.SOR ~ poly(fascia, 3), data = x)

  
  if(AICc(glm_1SOR) - 2 <= AICc(glm_2SOR)){
    glm_SOR = glm_1SOR 
  } else if(AICc(glm_2SOR) - 2 <= AICc(glm_3SOR)){
    glm_SOR = glm_2SOR 
  }else{glm_SOR = glm_3SOR}
  
  sum_SOR <- summary(glm_SOR)
  
  coef_model <- sum_SOR$coefficients
  
  nomi_coef <- c("Intercept", "p_interc", "Elev1", "p1", "Elev2", "p2", "Elev3", "p3") %>%
    .[1:(nrow(coef_model)*2)]
  
  coef_model <- coef_model %>%
    t() %>% as.data.frame() %>%  
    rowwise() %>% 
    map(., function(y) y[c(1,4)]) %>% 
    do.call(c,.) %>% 
    setNames(., nomi_coef)
  
  deh <- c(D2 = modEvA::Dsquared(glm_SOR), AIC = sum_SOR$aic, coef_model) %>%
    t() %>%
    as.data.frame() %>%
    cbind(size = x$size %>% unique())
  
  return(deh) 
}) %>% do.call(bind_rows, .) %>% 

  filter(size != 0.015, size != 0.03125) %>%
  arrange (factor(size, levels = c ('0.0625', '0.125', '0.25', 
                                    '0.5', '1', '2', '4', '8', '16')))
  

AIC_SOR <-  d2_AIC_SOR %>% mutate(across(where(is.numeric), .fns =  ~ round(.x, 3))) %>%
  mutate(p1 = case_when(
    p1 > 0.05 ~ "ns",
    p1 > 0.01 & p1 <= 0.05 ~ "*",
    p1 > 0.001  & p1 <= 0.01 ~ "**",
    p1 <= 0.001 ~"***" )) %>% 
  
 #mutate(p2 = case_when(
  #p2 > 0.05 ~ "ns",
   #p2 > 0.01 & p2 <= 0.05 ~ "*",
   #p2 > 0.001  & p2 <= 0.01 ~ "**",
   #p2 <= 0.001 ~"***" )) %>%
  
  
 # mutate(p3 = case_when(
  #  p3 > 0.05 ~ "ns",
   # p3 > 0.01 & p3 <= 0.05 ~ "*",
    #p3 > 0.001  & p3 <= 0.01 ~ "**",
    #p3 <= 0.001 ~"***" )) %>%
  
  
  mutate(p_interc = case_when(
    p_interc > 0.05 ~ "ns",
    p_interc > 0.01 & p_interc <= 0.05 ~ "*",
    p_interc > 0.001  & p_interc <= 0.01 ~ "**",
    p_interc <= 0.001 ~"***" )) %>%
  replace(is.na(.), "") %>% 
  
  gt(rowname_col = "size") %>% 
  
  tab_header (title = md(paste("**Total multisite beta Vs Elevation**", sep = ""))) %>% 
  
  tab_style(style = cell_text(style = "italic"),
            location = cells_stub()) %>% 
  cols_label(
    D2 = "D2", 
    p1 = "",
    #p2 = "", 
    #p3 = "",
    p_interc = "",
    AIC = "AIC") 
show(AIC_SOR)





predict_betamulti_SOR <-   map(group_sizes, function(x){
  
 # x= group_sizes [[1]]
 
  glm_1SOR <- MASS::glm.nb(formula = beta.SOR ~ poly(fascia, 1), data = x) 
  glm_2SOR <- MASS::glm.nb(formula = beta.SOR ~ poly(fascia, 2), data = x)
  glm_3SOR <- MASS::glm.nb(formula = beta.SOR ~ poly(fascia, 3), data = x)
  
  if(summary(glm_1SOR)$aic - 2 <= summary(glm_2SOR)$aic){
    glm_multi = glm_1SOR
  } else if(summary(glm_2SOR)$aic - 2 <= summary(glm_3SOR)$aic){
    glm_multi = glm_2SOR
  }else{glm_multi = glm_3SOR}
  
  #newdata <- list()
  
  newdata <- seq(min(x$fascia), max(x$fascia), length = 100) %>% 
    cbind.data.frame(fascia = ., size = rep(x$size %>% unique(), 100))
  
  newdata$beta.SOR <- predict(glm_multi, newdata = newdata, type = "response")
  
  return(newdata)
  
}) %>% do.call(bind_rows, .) %>% 
  mutate(size = as.numeric(size)) %>%
  arrange(size) %>% 
  filter(size != 0.015, size != 0.03125) %>% 
  mutate(size = as.factor(size))


plot_predict_multi <-   bm_ril_fascia %>% 
  mutate(size = as.numeric(size)) %>%
  arrange(size) %>% 
  filter(size != 0.015, size != 0.03125) %>%
  mutate(size = as.factor(size)) %>% 
  ggplot(aes (x = fascia, y = beta.SOR, col = size, group = size)) +
  geom_line( data = predict_betamulti_SOR,mapping = aes(x = fascia, col = size, group = size, y = beta.SOR), size = .9, alpha = 0.60) +
  geom_point(alpha = 0.2) +
  guides(color = guide_legend(title = "Sample size")) +
   theme_bw()+
  labs(x = "Elevation m a.s.l.", y = "Total beta multiple-site diversity (Sørensen index)") +
  scale_colour_brewer(palette = "Paired") 


ggsave(filename = "betamulti_predict.tiff", plot =  plot_predict_multi, 
       width = 180, height = 160, units = "mm", dpi = 300, compression="lzw")





############# D2 AIC BETA SIM
d2_AIC_SIM <-   map(group_sizes, function(x){
  #x <- group_sizes[[1]]
  
  #print(x$size %>% unique())
  
  x <- x %>% 
    drop_na() 
  
  
  
  glm_1SIM <- MASS::glm.nb(formula = beta.SIM ~ poly(fascia, 1), data = x) 
  glm_2SIM <- MASS::glm.nb(formula = beta.SIM ~ poly(fascia, 2), data = x)
  glm_3SIM <- MASS::glm.nb(formula = beta.SIM ~ poly(fascia, 3), data = x)
  
  
  if(AICc(glm_1SIM) - 2 <= AICc(glm_2SIM)){
    glm_SIM = glm_1SIM 
  } else if(AICc(glm_2SIM) - 2 <= AICc(glm_3SIM)){
    glm_SIM = glm_2SIM 
  }else{glm_SIM = glm_3SIM}
  
  sum_SIM <- summary(glm_SIM)
  
  coef_model <- sum_SIM$coefficients
  
  nomi_coef <- c("Intercept", "p_interc", "Poly1", "p1", "Poly2", "p2", "Poly3", "p3") %>%
    .[1:(nrow(coef_model)*2)]
  
  coef_model <- coef_model %>%
    t() %>% as.data.frame() %>%  
    rowwise() %>% 
    map(., function(y) y[c(1,4)]) %>% 
    do.call(c,.) %>% 
    setNames(., nomi_coef)
  
  deh <- c(D2 = modEvA::Dsquared(glm_SIM), AIC = sum_SIM$aic, coef_model) %>%
    t() %>%
    as.data.frame() %>%
    cbind(size = x$size %>% unique())
  
  return(deh) 
}) %>% do.call(bind_rows, .) %>% 
  
  filter(size != 0.015, size != 0.03125) %>%
  arrange (factor(size, levels = c ('0.0625', '0.125', '0.25', 
                                    '0.5', '1', '2', '4', '8', '16')))


AIC_SIM <-  d2_AIC_SIM %>% mutate(across(where(is.numeric), .fns =  ~ round(.x, 3))) %>%
  mutate(p1 = case_when(
    p1 > 0.05 ~ "ns",
    p1 > 0.01 & p1 <= 0.05 ~ "*",
    p1 > 0.001  & p1 <= 0.01 ~ "**",
    p1 <= 0.001 ~"***" )) %>% 
  
  # mutate(p2 = case_when(
  #   p2 > 0.05 ~ "ns",
  #   p2 > 0.01 & p2 <= 0.05 ~ "*",
  #   p2 > 0.001  & p2 <= 0.01 ~ "**",
  #   p2 <= 0.001 ~"***" )) %>%
  
  
  #mutate(p3 = case_when(
  #  p3 > 0.05 ~ "ns",
  #  p3 > 0.01 & p3 <= 0.05 ~ "*",
#  p3 > 0.001  & p3 <= 0.01 ~ "**",
#  p3 <= 0.001 ~"***" )) %>%


mutate(p_interc = case_when(
  p_interc > 0.05 ~ "ns",
  p_interc > 0.01 & p_interc <= 0.05 ~ "*",
  p_interc > 0.001  & p_interc <= 0.01 ~ "**",
  p_interc <= 0.001 ~"***" )) %>%
  replace(is.na(.), "") %>% 
  
  gt(rowname_col = "size") %>% 
  
  tab_header (title = md(paste("**Turnover multisite beta Vs Elevation**", sep = ""))) %>% 
  
  tab_style(style = cell_text(style = "italic"),
            location = cells_stub()) %>% 
  cols_label(
    D2 = "D2", 
    p1 = "",
    # p2 = "", 
    #p3 = "",
    p_interc = "",
    AIC = "AIC") 
show(AIC_SIM)



predict_betamulti_SIM <-   map(group_sizes, function(x){
  
  # x= group_sizes [[1]]
  
  glm_1SIM <- MASS::glm.nb(formula = beta.SIM ~ poly(fascia, 1), data = x) 
  glm_2SIM <- MASS::glm.nb(formula = beta.SIM ~ poly(fascia, 2), data = x)
  glm_3SIM <- MASS::glm.nb(formula = beta.SIM ~ poly(fascia, 3), data = x)
  
  if(summary(glm_1SIM)$aic - 2 <= summary(glm_2SIM)$aic){
    glm_multi = glm_1SIM
  } else if(summary(glm_2SIM)$aic - 2 <= summary(glm_3SIM)$aic){
    glm_multi = glm_2SIM
  }else{glm_multi = glm_3SIM}
  
  #newdata <- list()
  
  newdata <- seq(min(x$fascia), max(x$fascia), length = 100) %>% 
    cbind.data.frame(fascia = ., size = rep(x$size %>% unique(), 100))
  
  newdata$beta.SIM <- predict(glm_multi, newdata = newdata, type = "response")
  
  return(newdata)
  
}) %>% do.call(bind_rows, .) %>% 
  mutate(size = as.numeric(size)) %>%
  arrange(size) %>% 
  filter(size != 0.015, size != 0.03125) %>% 
  mutate(size = as.factor(size))


plot_predict_multi_SIM <-   bm_ril_fascia %>% 
  mutate(size = as.numeric(size)) %>%
  arrange(size) %>% 
  filter(size != 0.015, size != 0.03125) %>%
  mutate(size = as.factor(size)) %>% 
  ggplot(aes (x = fascia, y = beta.SIM, col = size, group = size)) +
  geom_point(alpha = 0.2) +
  guides(color = guide_legend(title = "Sample size")) +
  geom_line( data = predict_betamulti_SIM, mapping = aes(x = fascia, col = size, group = size, y = beta.SIM), size = .9, alpha = 0.60) +
  theme_bw()+
  labs(x = "Elevation m a.s.l.", y = "Beta multiple-site diversity turnover (Simpson index)") +
  scale_colour_brewer(palette = "Paired") 


ggsave(filename = "betamulti_SIM_predict.tiff", plot =  plot_predict_multi_SIM, 
       width = 180, height = 160, units = "mm", dpi = 300, compression="lzw")

#ggplot

bm_fasc <- bm_ril_fascia %>% 
  mutate(size = as.numeric(size)) %>% 
  filter(size != 0.015, size != 0.03125) %>% 
  mutate(size = as.factor(size)) %>% 
  ggplot(aes(y = beta.SOR, x = fascia, col = size, group = size))+
  geom_point() + 
  theme_bw() +
  geom_smooth(formula = y ~ poly(x,1), method = "lm", size = 1, alpha = 0.12) +
  scale_color_brewer(palette = "Pastel1") +
  labs ( x = "Elevation", y = "Beta Total")

ggsave(filename = "betamul Vs Elevation.tiff", plot =  bm_fasc , 
       width = 300, height = 150, units = "mm", dpi = 300, compression="lzw")
 # facet_wrap(~beta, scales = "free_y") 

###################nested 

d2_AIC_SNE <-   map(group_sizes, function(x){
  #x <- group_sizes[[1]]
  
  #print(x$size %>% unique())
  
  x <- x %>% 
    drop_na() 
  
  
  
  glm_1SNE <- MASS::glm.nb(formula = beta.SNE ~ poly(fascia, 1), data = x) 
  glm_2SNE <- MASS::glm.nb(formula = beta.SNE ~ poly(fascia, 2), data = x)
  glm_3SNE <- MASS::glm.nb(formula = beta.SNE ~ poly(fascia, 3), data = x)
  
  
  if(AICc(glm_1SNE) - 2 <= AICc(glm_2SNE)){
    glm_SNE = glm_1SNE 
  } else if(AICc(glm_2SNE) - 2 <= AICc(glm_3SNE)){
    glm_SNE = glm_2SNE 
  }else{glm_SNE = glm_3SNE}
  
  sum_SNE <- summary(glm_SNE)
  
  coef_model <- sum_SNE$coefficients
  
  nomi_coef <- c("Intercept", "p_interc", "Poly1", "p1", "Poly2", "p2", "Poly3", "p3") %>%
    .[1:(nrow(coef_model)*2)]
  
  coef_model <- coef_model %>%
    t() %>% as.data.frame() %>%  
    rowwise() %>% 
    map(., function(y) y[c(1,4)]) %>% 
    do.call(c,.) %>% 
    setNames(., nomi_coef)
  
  deh <- c(D2 = modEvA::Dsquared(glm_SNE), AIC = sum_SNE$aic, coef_model) %>%
    t() %>%
    as.data.frame() %>%
    cbind(size = x$size %>% unique())
  
  return(deh) 
}) %>% do.call(bind_rows, .) %>% 
  
  filter(size != 0.015, size != 0.03125) %>%
  arrange (factor(size, levels = c ('0.0625', '0.125', '0.25', 
                                    '0.5', '1', '2', '4', '8', '16')))


AIC_SNE <-  d2_AIC_SNE %>% mutate(across(where(is.numeric), .fns =  ~ round(.x, 3))) %>%
  mutate(p1 = case_when(
    p1 > 0.05 ~ "ns",
    p1 > 0.01 & p1 <= 0.05 ~ "*",
    p1 > 0.001  & p1 <= 0.01 ~ "**",
    p1 <= 0.001 ~"***" )) %>% 
  
  # mutate(p2 = case_when(
  #   p2 > 0.05 ~ "ns",
  #   p2 > 0.01 & p2 <= 0.05 ~ "*",
  #   p2 > 0.001  & p2 <= 0.01 ~ "**",
  #   p2 <= 0.001 ~"***" )) %>%
  
  
  #mutate(p3 = case_when(
  #  p3 > 0.05 ~ "ns",
  #  p3 > 0.01 & p3 <= 0.05 ~ "*",
#  p3 > 0.001  & p3 <= 0.01 ~ "**",
#  p3 <= 0.001 ~"***" )) %>%


mutate(p_interc = case_when(
  p_interc > 0.05 ~ "ns",
  p_interc > 0.01 & p_interc <= 0.05 ~ "*",
  p_interc > 0.001  & p_interc <= 0.01 ~ "**",
  p_interc <= 0.001 ~"***" )) %>%
  replace(is.na(.), "") %>% 
  
  gt(rowname_col = "size") %>% 
  
  tab_header (title = md(paste("**Nestedness multisite beta Vs Elevation**", sep = ""))) %>% 
  
  tab_style(style = cell_text(style = "italic"),
            location = cells_stub()) %>% 
  cols_label(
    D2 = "D2", 
    p1 = "",
    # p2 = "", 
    #p3 = "",
    p_interc = "",
    AIC = "AIC") 
show(AIC_SNE)



predict_betamulti_SNE <-   map(group_sizes, function(x){
  
  # x= group_sizes [[1]]
  
  glm_1SNE <- MASS::glm.nb(formula = beta.SNE ~ poly(fascia, 1), data = x) 
  glm_2SNE <- MASS::glm.nb(formula = beta.SNE ~ poly(fascia, 2), data = x)
  glm_3SNE <- MASS::glm.nb(formula = beta.SNE ~ poly(fascia, 3), data = x)
  
  if(summary(glm_1SNE)$aic - 2 <= summary(glm_2SNE)$aic){
    glm_multi = glm_1SNE
  } else if(summary(glm_2SNE)$aic - 2 <= summary(glm_3SNE)$aic){
    glm_multi = glm_2SNE
  }else{glm_multi = glm_3SNE}
  
  #newdata <- list()
  
  newdata <- seq(min(x$fascia), max(x$fascia), length = 100) %>% 
    cbind.data.frame(fascia = ., size = rep(x$size %>% unique(), 100))
  
  newdata$beta.SNE <- predict(glm_multi, newdata = newdata, type = "response")
  
  return(newdata)
  
}) %>% do.call(bind_rows, .) %>% 
  mutate(size = as.numeric(size)) %>%
  arrange(size) %>% 
  filter(size != 0.015, size != 0.03125) %>% 
  mutate(size = as.factor(size))


plot_predict_multi_NE <-   bm_ril_fascia %>% 
  mutate(size = as.numeric(size)) %>%
  arrange(size) %>% 
  filter(size != 0.015, size != 0.03125) %>%
  mutate(size = as.factor(size)) %>% 
  ggplot(aes (x = fascia, y = beta.SNE, col = size, group = size)) +
  geom_point(alpha = 0.2) +
  guides(color = guide_legend(title = "Sample size")) +
  geom_line( data = predict_betamulti_SNE, mapping = aes(x = fascia, col = size, group = size, y = beta.SNE), size = .9, alpha = 0.60) +
  theme_bw()+
  labs(x = "Elevation m a.s.l.", y = "Beta multiple-site diversity nestedness") +
  scale_colour_brewer(palette = "Paired") 


ggsave(filename = "betamulti_SNE_predict.tiff", plot =  plot_predict_multi_NE, 
       width = 180, height = 160, units = "mm", dpi = 300, compression="lzw")

#ggplot

bm_fasc <- bm_ril_fascia %>% 
  mutate(size = as.numeric(size)) %>% 
  filter(size != 0.015, size != 0.03125) %>% 
  mutate(size = as.factor(size)) %>% 
  ggplot(aes(y = beta.SOR, x = fascia, col = size, group = size))+
  geom_point() + 
  theme_bw() +
  geom_smooth(formula = y ~ poly(x,1), method = "lm", size = 1, alpha = 0.12) +
  scale_color_brewer(palette = "Pastel1") +
  labs ( x = "Elevation", y = "Beta Total")

ggsave(filename = "betamul Vs Elevation.tiff", plot =  bm_fasc , 
       width = 300, height = 150, units = "mm", dpi = 300, compression="lzw")
# facet_wrap(~beta, scales = "free_y") 





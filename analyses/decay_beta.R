


coef_sp <- map (beta_decay_size, function(x){
  
  #x = beta_decay_size [[8]]
  
  # y = beta_dist_split[[1]]
  #deh <- list()
  
  glmfit1 <- glm (formula = value ~ poly(dist_sp, 1), data = x, family = 'binomial')
  glmfit2 <- glm (formula = value ~ poly(dist_sp, 2), data = x, family = 'binomial')
  
  if(summary(glmfit1)$aic - 2 <= summary(glmfit2)$aic){
    glm_beta = glmfit1
  } else {glm_beta = glmfit2}
  
  sum_glm <- summary(glm_beta)
  
  coef_model <- sum_glm$coefficients 
  
  nomi_coef <- c("Int", "p_interc", "Elev1", "p1", "Elev2", "p2") %>%
    .[1:(nrow(coef_model)*2)]
  
  coef_model <- coef_model %>%
    t() %>% as.data.frame() %>%  
    rowwise() %>% 
    map(., function(y) y[c(1,4)]) %>% 
    do.call(c,.) %>% 
    setNames(., nomi_coef)
  # deh <- list()
  
  deh <- c(D2 = modEvA::Dsquared(glm_beta), AIC = sum_glm$aic, coef_model) %>%
    t() %>%
    as.data.frame() %>%
    cbind(size = x$size %>% unique(), beta = x$beta %>% unique())
  
  
  return(deh)   
  
}) %>% do.call(bind_rows,.)




coef_sp$title_beta <-  ifelse (coef_sp$beta == 'beta.sim', 'Beta Turnover',
                                ifelse(coef_sp$beta == 'beta.sne', 'Beta Nestedness', 'Beta Total')) 





coef_sp %>% 
  select(beta, title_beta) %>% 
  distinct(beta, .keep_all = T) %>% 
  group_by(beta) %>% 
  group_split() %>% 
  
  map(., function(x){
    
    co_d2 <- coef_sp %>% 
      filter(beta == x$beta) %>% 
      arrange (factor(size, levels = c ('0.0625', '0.125', '0.25', 
                                        '0.5', '1', '2', '4', '8', '16'))) %>% 
      mutate(across(where(is.numeric), .fns =  ~ round(.x, 3))) %>% 
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
      
      mutate(p_interc = case_when(
        p_interc > 0.05 ~ "ns",
        p_interc > 0.01 & p_interc <= 0.05 ~ "*",
        p_interc > 0.001  & p_interc <= 0.01 ~ "**",
        p_interc <= 0.001 ~"***" )) %>%
      replace(is.na(.), "") %>% 
      select(- beta) %>%
      
      gt(rowname_col = "size") %>% 
      
      tab_header (title = md(paste("**", x$title_beta, " Vs ", "Spatial distance", "**", sep = ""))) %>% 
      
      tab_style(style = cell_text(style = "italic"),
                location = cells_stub()) %>% 
      cols_label(
        D2 = "D2", 
        p1 = "",
        p2 = "", 
        p_interc = "",
        AIC = "AIC") %>% 
      
      cols_hide (c(title_beta))   
    # gtsave(filename = paste(x$title_beta, " Vs ", x$title_dist, ".png"))
    
    
  })



predict_beta <- map (beta_decay_size, function(x){
  beta_dist_split <- x %>%
    group_by (beta) %>% 
    group_split() %>% 

  glmfit1 <- glm (formula = value ~ poly(dist_sp, 1), data = x, family = 'binomial')
  glmfit2 <- glm (formula = value ~ poly(dist_sp, 2), data = x, family = 'binomial')
  
  
  
  if(summary(glmfit1)$aic - 2 <= summary(glmfit2)$aic){
    glm_beta = glmfit1
  } else {glm_beta = glmfit2}
  
  newdata <- seq(min(x$dist_sp), max(x$dist_sp), length = 100) %>% 
    cbind.data.frame(dist_sp = ., size = rep(x$size %>% unique(),100))
  
  newdata$value <- predict(glm_beta, newdata = newdata, type = "response")
  
  return(newdata)
  
}) %>% do.call(bind_rows, .) %>% 
  
  
  mutate(size = as.numeric(size)) %>%
  arrange(size) %>% 
  filter(size != 0.015, size != 0.03125) %>% 
  mutate(size = as.factor(size))

beta_plot <- gamma_fascia %>% 
  mutate(size = as.numeric(size)) %>%
  arrange(size) %>% 
  filter(size != 0.015, size != 0.03125) %>%
  mutate(size = as.factor(size)) %>% 
  ggplot(aes (x = dist_value, y = beta, col = size, group = size)) +
  geom_point(alpha = 0.2) +
  geom_line( data = predict_gamma, mapping = aes(x = Fascia, col = size, group = size, y = gamma),
             size = .9, alpha = 0.60) +
  theme_bw()+
  labs(x = "Elevation", y = "Gamma") +
  scale_colour_brewer(palette = "Paired") 
ggsave(filename = "gammadive.tiff", plot =  gamma_plot, 
       width = 180, height = 160, units = "mm", dpi = 300, compression="lzw")



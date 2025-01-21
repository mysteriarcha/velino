library(raster)
library(tidyverse)
library(vegan)
library(gt)
library(modEvA)
library(tmap)
library(sf)


ril_vel <- read.csv("velino_firstcycle.csv", row.names = 1) %>% 
   t() %>% as.data.frame() %>% 
  add_column(ril = word(rownames(.), 1, sep = '_'), 
             size = word(rownames(.), 2, sep = '_') ) 

coords <- st_read("coords_rilievi.shp") %>% 
  st_transform("EPSG:32632")

vel_env <- read.csv ("velino_environment.csv", row.names = 1, sep = ";") %>% 
  cbind(coords$geometry)


ril_fascia <- vel_env %>% 
  rownames_to_column("ril") %>% 
  select(ril, quota) %>% 
  left_join(ril_vel) 
  

  speciesrich <- ril_fascia %>% 
    select(-ril) %>%
    filter(size != 0.015, size != 0.03125) %>% 
    mutate(alpha = select(., `Acer campestre`:`Xeranthemum inapertum`) %>% rowSums(na.rm = TRUE)) %>% 
    select(size, quota, alpha)
    
 
 #source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

  speciesrich  %>% 
    mutate(size = as.numeric(size)) %>%
    arrange(size) %>% 
    filter(size != 0.015, size != 0.03125) %>% 
    mutate(size = as.factor(size)) %>% 
    ggplot(aes (x = quota, y = alpha, col = size, group = size)) +
  #  geom_point() +
    geom_smooth(formula =  y ~ poly(x,1), method = "lm", size = 1.5, alpha = 0.20) +
 #   scale_y_discrete(limits = factor(0:80))+
    theme_bw()+
    xlab("Elevation") +
    ylab("Alpha diversity") +
    scale_color_brewer(palette = "Set1") 
  
  

 ##split SR x grain size
  
   sr_split <-  speciesrich %>% 
    group_by(size) %>% 
    group_split() 
   
    
  ########D2 e AIC
     

  d2_AIC_alpha <-   map(sr_split, function(x){
  
    #x <- sr_split[[9]]
    #deh <- list()
    
    glm_1alpha <- glm (formula = alpha ~ poly(quota, 1), data = x, family = poisson)
    glm_2alpha <- glm (formula = alpha ~ poly(quota, 2), data = x, family = poisson)
    glm_3alpha <- glm (formula = alpha ~ poly(quota, 3), data = x, family = poisson)
    
    
    if(summary(glm_1alpha)$aic - 2 <= summary(glm_2alpha)$aic){
      glm_alpha = glm_1alpha
    } else if(summary(glm_2alpha)$aic - 2 <= summary(glm_3alpha)$aic){
        glm_alpha = glm_2alpha
    }else{glm_alpha = glm_3alpha}
    
    sum_glm <- summary(glm_alpha)
    
    coef_model <- sum_glm$coefficients 
    
    names_coef <- c("Int", "p_interc", "Elev1", "p1", "Elev2", "p2", "Elev3", "p3") %>%
      .[1:(nrow(coef_model)*2)]
    
    coef_model <- coef_model %>%
      t() %>% as.data.frame() %>%  
      rowwise() %>% 
      map(., function(y) y[c(1,4)]) %>% 
      do.call(c,.) %>% 
      setNames(., names_coef)
    
 # deh <- list()
    
    deh <- c(D2 = modEvA::Dsquared(glm_alpha), AIC = sum_glm$aic, coef_model) %>%
      t() %>%
      as.data.frame() %>%
      cbind(size = x$size %>% unique())
    
    return(deh) 
  
  })
  
table <-  d2_AIC_alpha %>% do.call(bind_rows, .) %>% 
    arrange (factor(size, levels = c ('0.0625', '0.125', '0.25', 
                                    '0.5', '1', '2', '4', '8', '16')))  

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
sr_plot <-   sr %>% 
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
     
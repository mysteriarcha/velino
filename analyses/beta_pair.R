
library(tidyverse)
library(betapart)
library(sf)
library(purrr)
library(gt)


ril_vel <- read.csv("velino_firstcycle.csv", row.names = 1) %>% 
  t() %>% as.data.frame() %>% 
  add_column(ril = word(rownames(.), 1, sep = '_'), 
             size = word(rownames(.), 2, sep = '_') ) 

coords <- st_read("coords_rilievi.shp")

vel_env <- read.csv ("velino_environment.csv", row.names = 1, sep = ";") %>% 
  cbind(coords$geometry)




####beta.pair per fascia
  
  bp_ril_fascia <- vel_env %>% 
              rownames_to_column("ril") %>% 
              select(ril, Fascia) %>% 
              left_join(ril_vel) %>% 
              dplyr::select(-ril) %>% 
              group_by(Fascia, size) %>% 
              group_split() %>% 
              map(~ beta.pair(.x %>% dplyr::select(-Fascia, -size)) %>% do.call(cbind, .) %>% 
                     cbind(., fascia = .x$Fascia %>% unique(), size = .x$size %>% unique())) %>% 
                     do.call(rbind,.) %>% 
    as.data.frame() 

####ggplot
#bp_ril_fascia %>% 
#  mutate_across.(where(is.character), as.numeric) %>%
#  filter(size != 0.015, size != 0.03125) %>% 
#  ggplot(aes(y = beta.sor, x = fascia, col = size, group = size))+
#  scale_color_viridis_c() +
#  geom_smooth()

geom_distance <- sf::st_as_sf(vel_env) %>% sf::st_distance()
veg_dist <- ril_vel %>% select(-c(ril, size)) %>% tapply(., INDEX = ril_vel$ril, FUN = sum) %>% vegdist()

tmp <- expand.grid(ril_vel$ril, ril_vel$size)

which(paste0(tmp[,1], tmp[,2]), paste0(ril_vel$ril, ril_vel$size))

make_veg_dist_by_size <- function(my_size){
  
  ind <- which(ril_vel$size == my_size)
  
  veg_dist <- ril_vel %>%
    .[ind, ] %>% 
    dplyr::select(-c(ril, size)) %>% 
    tapply(., INDEX = ril_vel$ril[ind], FUN = sum) %>% 
    vegdist(method = "euclidean")
  return(veg_dist)
}

my_scales <- unique(ril_vel$size)
scale_effects <- vector("list", length(my_scales))
for(i in seq_along(scale_effects)){
  scale_effects[[i]] <- mantel(geom_distance, make_veg_dist_by_size(my_scales[[i]]))
}

geom_vs_veg <-
  lapply(scale_effects, "[[", "statistic") %>%
  Reduce(c, .) %>% 
  data.frame(scale = as.numeric(my_scales), 
             correlation = .,
             significance = sapply(scale_effects, "[[", "signif") < 0.05)

lm(correlation ~ log(scale,2), geom_vs_veg) %>% summary()

geom_vs_veg %>% 
  ggplot(aes(x = round(log(scale,2),0), y = correlation, alpha = significance))+
  geom_point() +
  geom_path(group = 1) +
  scale_x_continuous(breaks = round(log(as.numeric(my_scales), 2), 0))+
  labs(x = "Sampling Scale (log2-scale)", y = "Cor. Vegetational vs Geographic distances") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20, face = "bold")
  )

##beta.pair X grain size + decay
        
bp_tot <-  vel_env %>% 
          rownames_to_column("ril") %>% 
          select(ril,quota,geometry) %>% 
          left_join(ril_vel) %>% 
            dplyr::select(-ril) %>% 
            group_by(size) %>% 
            group_split()  

    beta_decay <-   map(bp_tot, function(x){
        
        beta_p <- beta.pair(x %>% dplyr::select(-size,-quota,-geometry)) %>% 
          do.call(cbind,.) %>% 
          cbind(., size = x [496] %>% unique())
        
          dist_h <- dist(x [1]) %>% 
                     as.numeric() 
          
          dist_sp <- st_as_sf(x [2]) %>% 
                     st_distance() %>% 
                      as.dist()  %>% 
                      as.numeric()
           
                    
          beta_dist <- cbind(beta_p, dist_h, dist_sp) %>% 
            
          
       return(beta_dist)
         }) %>% do.call(rbind,.) %>% unique()
    
    hist(beta_decay$beta.sor)
    
    ##pivot_longer for plotting throught ggplot 
    
    beta_pivo <- beta_decay %>% 
                pivot_longer(cols = starts_with("beta"), names_to = "beta", values_to = "value") %>% 
                 mutate(size = as.numeric (size))
    
    
  hist(beta_pivo$value)  
    
  bp_pdist <- beta_decay %>% 
   mutate(size = as.numeric(size)) %>%
    arrange(size) %>% 
    filter(size != 0.015, size != 0.03125) %>% 
    mutate(size = as.factor(size)) %>% 
   ggplot (aes(x = dist_sp, y = beta.sor, group = size, col = size)) +
    theme_bw ()+
     geom_smooth( formula = y ~ poly(x,2), size = 1, alpha = 0.12) +
    scale_color_brewer(palette = "Pastel1")+
    labs(x = "Spatial distance", y = "Beta Total") 
  #facet_wrap(~beta, scales = "free_y")
  
  #ggsave(filename = "beta.tot_hdist.tiff", plot =  bp_pdist  , 
   #      width = 210, height = 150, units = "mm", dpi = 300, compression="lzw")
    

   
   
   #####D2 dist altitudinale e spaziale Vs beta.sor
   
  beta_decay_size <- beta_decay %>% 
                          group_by(size) %>% 
                           group_split() 
   

  
 pivot_beta_dist <- map (beta_decay_size, function(x){
    
    pivot <- x %>% 
               pivot_longer(starts_with(("beta")))
    
    pivot_2 <- pivot %>% 
               pivot_longer(starts_with("dist"), names_repair = c("unique"))

    colnames(pivot_2) <- c("size", "beta", "beta_value", "dist", "dist_value") 
   
  return (pivot_2)
    
    })
  
 
 
   coef_1_2 <- map (pivot_beta_dist, function(x){
  
    
     
       beta_dist_split <- x %>%
         group_by(beta, dist) %>% 
         group_split() 
      
       
       map(., function(y){ 
         
   
         #deh <- list()
         glmfit1 <- glm (formula = beta_value ~ poly(dist_value, 1), data = y, family = 'binomial')
         glmfit2 <- glm (formula = beta_value ~ poly(dist_value, 2), data = y, family = 'binomial')
         
         if(summary(glmfit1)$aic - 2 <= summary(glmfit1)$aic (glmfit2)){
           glm_beta = glmfit1
         } else {glm_beta = glmfit2}
                  
         sum_glm <- summary(glm_beta)
         
         coef_model <- sum_glm$coefficients 
         
         names_coef <- c("Int", "p_interc", "Elev1", "p1", "Elev2", "p2") %>%
           .[1:(nrow(coef_model)*2)]
         
         coef_model <- coef_model %>%
           t() %>% as.data.frame() %>%  
           rowwise() %>% 
           map(., function(y) y[c(1,4)]) %>% 
           do.call(c,.) %>% 
           setNames(., names_coef)
     
         
         deh <- c(D2 = modEvA::Dsquared(glm_beta), AICc = sum_glm$aic, coef_model) %>%
           t() %>%
           as.data.frame() %>%
           cbind(size = x$size %>% unique(), beta = y$beta %>% unique(), dist = y$dist %>% unique())
         
                    
           return(deh)   
                          
       }) %>% do.call(bind_rows,.)
  }) %>% do.call(bind_rows,.)

  

    coef_1_2$title_beta <-  ifelse (coef_1_2$beta == 'beta.sim', 'Beta Turnover',
                             ifelse(coef_1_2$beta == 'beta.sne', 'Beta Nestedness', 'Beta Total')) 
    
    coef_1_2$title_dist <- ifelse (coef_1_2$dist == 'dist_h', ' Altitudinal distance', 'Spatial distance') 
    

                                 
   coef_1_2 %>% 
     select(beta, dist, title_beta, title_dist) %>% 
     distinct(beta, dist, .keep_all = T) %>% 
    group_by(beta, dist) %>% 
    group_split() %>% 
     
     map(., function(x){
       
   co_d2 <- coef_1_2 %>% 
     filter(size != 0.015, size != 0.03125) %>%
     filter(beta == x$beta & dist == x$dist) %>% 
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
     select(- beta, - dist) %>%
       
     gt(rowname_col = "size") %>% 
       
     tab_header (title = md(paste("**", x$title_beta, " Vs ", x$title_dist, "**", sep = ""))) %>% 
     
     tab_style(style = cell_text(style = "italic"),
               location = cells_stub()) %>% 
     cols_label(
       D2 = "D2", 
       p1 = "", 
       p2 = "",
       p_interc = "",
       AIC = "AIC") %>% 
     
     cols_hide (c(title_beta, title_dist, Elev2,p2))     
    # gtsave(filename = paste(x$title_beta, " Vs ", x$title_dist, ".png"))
  
   
    })


   predict_beta <- map (pivot_beta_dist, function(x){
    
      beta_dist_split <- x %>%
       group_by (beta, dist) %>% 
       group_split() %>% 
       
      map (., function (y) {
      
       glmbeta <- glm (formula = beta_value ~ poly(dist_value, 1), data = y, family = 'binomial')
    
          newdata <- seq(min(y$dist_value), max(y$dist_value), length = 100) %>% 
         cbind.data.frame(dist_value = ., size = rep(y$size %>% unique(),100))
          
          newdata$value <- predict(glmbeta, newdata = newdata, type = "response")
          newdata$dist <- y$dist %>% unique()
          newdata$beta <- y$beta %>% unique()
          return(newdata)
      }) %>% do.call (bind_rows,.)
      
      return(beta_dist_split)
       
   }) %>% do.call(bind_rows, .) %>% 
     
     
     mutate(size = as.numeric(size)) %>%
     arrange(size) %>% 
     filter(size != 0.015, size != 0.03125) %>% 
     mutate(size = as.factor(size))
   
  new_data <-  beta_decay %>% 
    filter (size != 0.015, size != 0.03125) %>%
  arrange (factor(size, levels = c ('0.0625', '0.125', '0.25', 
                                    '0.5', '1', '2', '4', '8', '16'))) 
   
 h_beta <-    predict_beta %>% 
     filter(dist == "dist_h" & beta == "beta.sor")
     
     
    beta_plot <- ggplot(h_beta, aes (x = dist_value, y = value, col = size, group = size)) +
     geom_line(data = h_beta, mapping = aes(x = dist_value, y = value, col = size, group = size), size = .9, alpha = 0.60) +
       geom_point(data = new_data, mapping = aes( x = dist_h,y = beta.sor,group=size, col=size), size = .2, alpha = 0.1) +
     theme_bw()+
     labs(x = "Altitudinal distance", y = "Beta Total") +
     scale_colour_brewer(palette = "Paired") 
   
   #ggsave(filename = "beta_h.tiff", plot =  beta_plot, 
    #      width = 210, height = 150, units = "mm", dpi = 300)
   

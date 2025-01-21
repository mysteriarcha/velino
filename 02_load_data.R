## Script to load the data required for the analysis of Velino vegetation


#############################################################################
# 1. General use dataframes -----------------------------------------------

## Vegetation data (Plot x Species presence-absence matrix)
ril_vel <- read.csv("velino_firstcycle.csv", row.names = 1)
ril_vel %<>%
  t() %>%
  as.data.frame() %>% 
  # Add columns to have the plot identity and the area of each sample
  add_column(
    ril = word(rownames(.), 1, sep = '_'), 
    size = word(rownames(.), 2, sep = '_')
    )

# Load coordinates for the spatial analyses
coords <- 
  st_read("coords_rilievi.shp") %>% 
  # The EPSG:32632 is as the EPSG:4362 but in meters instead of decimal degrees
  st_transform("EPSG:32632")

# Load environmental variables and put them into an sf format
vel_env <- 
  read.csv ("velino_environment.csv", row.names = 1, sep = ";") %>% 
  cbind(coords$geometry) %>% 
  sf::st_as_sf()

# Whenever we want to get rid of the geographic metadata we can just
# apply an st_drop_geometry() on the vel_env object


#############################################################################
### 2. Alpha diversity dataframes -------------------------------------------

ril_fascia <-
  vel_env %>%
  as.data.frame() %>%
  rownames_to_column("ril") %>% 
  dplyr::select(ril, quota) %>% 
  left_join(ril_vel)


speciesrich <- 
  ril_fascia %>% 
  dplyr::select(-ril) %>%
  filter(size != 0.015, size != 0.03125) %>% 
  mutate(alpha = {  
    dplyr::select(., `Acer campestre`:`Xeranthemum inapertum`) %>% 
      rowSums(na.rm = TRUE)
    }
  ) %>% 
  dplyr::select(size, quota, alpha)

#############################################################################
### 3. Beta diversity dataframes --------------------------------------------

bp_ril_fascia <- 
  vel_env %>%
  as.data.frame() %>%
  rownames_to_column("ril") %>% 
  dplyr::select(ril, Fascia) %>% 
  left_join(ril_vel) %>% 
  dplyr::select(-ril) %>% 
  group_by(Fascia, size) %>% 
  group_split() %>% 
  map(~ {
    beta.pair(
      .x %>% 
        dplyr::select(-Fascia, -size)
      ) %>% 
      do.call(cbind, .) %>% 
      cbind(., fascia = .x$Fascia %>% unique(), size = .x$size %>% unique()
      )
    }
  ) %>% 
  do.call(rbind, .) %>% 
  as.data.frame() 


#############################################################################
### 4. Gamma diversity dataframes -------------------------------------------

ril_gamma <- 
  vel_env %>% 
  rownames_to_column("ril") %>% 
  dplyr::select(ril, Fascia) %>% 
  left_join(ril_vel) %>% 
  st_drop_geometry()

gamma_fascia <-
  ril_gamma %>% 
  dplyr::select(-ril) %>% 
  group_by(size, Fascia) %>% 
  summarise(across(.fns =  ~ max(.x,na.rm = T)), .groups = "keep") %>%
  ungroup() %>% 
  filter(size != 0.015, size != 0.03125) %>% 
  mutate(gamma = dplyr::select(., `Acer campestre`:`Xeranthemum inapertum`) %>% rowSums(na.rm = TRUE)) %>% 
  dplyr::select(size,  Fascia, gamma) 

gam_split <-  
  gamma_fascia %>% 
  group_by(size) %>% 
  group_split() 


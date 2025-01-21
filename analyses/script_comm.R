

library(tidyverse)
library(dplyr)
library(betapart)
library(devtools)

# Carico i dati -----------------------------------------------------------


ril_vel <- read.csv("velino_primociclo.csv", row.names = 1) #carico la matrice di comunita per i rilievi del 1 ciclo
vel_env <- read.csv ("velino_env.csv", row.names = 1, sep = ";") #carico i dati stazionali e climatici per i rilievi 


# Perparazione dei dati per le analisi ------------------------------------

ril_vel_t <- ril_vel %>%  
t() %>% as.data.frame() %>%  #con t faccio la trasposta della matrice e sposto i rilievi sulle righe e le specie sulle colonne
  add_column(ril = word(rownames(.), 1, sep = '_'), # la colonna ril mi sara utile in seguito per il join
             size = word(rownames(.), 2, sep = '_') ) #con add_column aggiugo una colonna con il codice identificatico dei rilievi e una con la rispettiva taglia


#Join tra i due df

ril_quota <- vel_env %>% 
  rownames_to_column("ril") %>%  #da righe a colonne
  select(ril, quota) %>% #seleziono solo il nome dei rilievi e la quota 
  left_join(ril_vel_t) #faccio il join con il df di comunita


# Alpha diversita lungo il gradiente altitudinale -------------------------



#Calcolo l'alpha diversita lungo tutto il gradiente altitudinale e per tutte le taglie dei plot

SR_tot <- ril_quota %>% 
  select(-ril) %>% #rimuovo il nome dei rilievi dalla matrice
  mutate(alpha = select(., `Acer campestre`:`Xeranthemum inapertum`) %>% rowSums(na.rm = TRUE)) %>% #creo una nuova colonna con l'alpha div, 
  #calcolata come somma dei valori delle righe selezionando solo le colonne riportanti la presenza/assenza delle specie
  select(size,quota, alpha) 


#Calcolo l'alpha diversita solo per i rilievi al di sopra dei 1482m e con taglia uguale a 16m

SR_middle <- ril_quota %>% 
  select(-ril) %>%
  filter(quota > 1482 & size == 16) %>% #seleziono solo le quote al di sopra dei 1482m e i plot con taglia 16m
  mutate(alpha = select(., `Acer campestre`:`Xeranthemum inapertum`) %>% rowSums(na.rm = TRUE)) %>% 
  select(quota, alpha) 


# Beta multi per fascia altitudinale  -------------------------------------



#Calcolo la beta div per fascia altitudinale 

bm_ril_fascia <- vel_env %>% 
  rownames_to_column("ril") %>% #colonna ril utile per il join e per aggiungere la fascia alla matrice di comunita
  select(ril,Fascia) %>% #seleziono le colonne 
  left_join(ril_vel_t,.) %>% #faccio il join
  dplyr::select(-ril)  #rimuovo la colonna ril che non mi e piu utile dopo il join


#Calcolo la beta multi per tutte le fasce e taglie (NB: la funzione restituisce 3 valori di beta)

bm_tot <- bm_ril_fascia %>% 
  group_by(Fascia, size) %>% #raggruppo i rilievi per fascia altitudinale e per taglia
  group_split() %>% 
  map_dfr(~ beta.multi(.x %>% dplyr::select(-Fascia, -size)) %>% #calcolo la beta con la funzione beta.multi
            c(., fascia = .x$Fascia %>% unique(), size = .x$size %>% unique())) 


###Beta multi calcolata solo per la fascia 2100 m

bm_2100 <- bm_ril_fascia %>% 
  filter(Fascia == 2100) %>%
  select(-Fascia)  %>% 
  group_by(size) %>% 
  group_split() %>% 
  map_dfr(~ beta.multi(.x %>% dplyr::select(-size)) %>% #calcolo la beta con la funzione beta.multi
            c(., size = .x$size %>% unique()))


# ggplot ------------------------------------------------------------------

#alpha div lungo il gradiente altitudinale 

alpha <- SR %>% 
  mutate(size = as.factor(size)) %>% #da caratteriale a fattoriale, trasformazione utile per ordinare i rilievi dalla taglia piu piccola
  arrange(size) %>% #ordino i rilievi dalla taglia piu piccola alla piu grande
  filter(size != 0.015, size != 0.03125) #filtro le taglie piu piccole 


#alpha plot

  ggplot(alpha, aes (x = quota, y = alpha, col = size, group = size)) + #ggplot 
  #geom_point() + #aggiungo i valori di ricchezza dei rilievi
  geom_smooth(formula =  y ~ poly(x,1), method = "lm", size = 1.5, alpha = 0.20) +
  theme_bw()+ #cambio lo sfondo del grafico
  scale_color_brewer(palette = "Pastel1") + ## scala di colore 
  xlab ('Elevation') +
  ylab ('Alpha')
  
#beta multi lungo le fasce altitudinali  
  
#beta multi plot, plotto la componente beta SOR, cioe la beta multi totale

 bm <- bm_tot %>% 
   mutate(size = as.factor(size)) %>% ##da caratteriale a fattoriale, trasformazione utile per ordinare i rilievi dalla taglia piu piccola
   arrange(size) %>% ##ordino i rilievi dalla taglia piu piccola alla piu grande
   filter(size != 0.015, size != 0.03125) %>% 
   drop_na #rimuovo le righe con NaN

   ggplot(bm, aes(y = beta.SOR, x = fascia,  col = size, group = size))+
      geom_point() + 
      theme_bw() +
      geom_smooth(formula = y ~ poly(x,1), method = "lm", size = 1, alpha = 0.09) +
      scale_color_brewer(palette = "Pastel1") +
      labs ( x = "Elevation", y = "Beta Total")
  

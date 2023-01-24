
## mapa exploratorio dos sitios

initial_map_function <- function (df1, df2) { 
  
  # mapa mundi
  world <- ne_countries(scale = "medium", returnclass = "sf")
  
  require(ggplot2)
  library(ggrepel)
  
  # cortar o mapa para ver a america do Sul e parte da central
  wm <- ggplot() + 
    geom_sf (data=world, size = 0.1, 
             fill= "gray90",colour="gray90") +
    coord_sf (xlim = c(-50, -30),  ylim = c(-28, -1), expand = FALSE) +
    theme_bw() + xlab(NULL) + ylab(NULL) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),axis.ticks.x=element_blank(),
          axis.text.y = element_blank(),axis.ticks.y=element_blank(),
          title = element_text(size=8)) 
  
  jitter <- position_jitter(width = 0.2, height = 0.5)
  
  bentos_coord <- wm + geom_point(data=df1,aes (x= Lon, y=Lat),
                                  stroke=1,shape=1, size=1, 
                                  position = jitter,col="red") 
  
  peixes_coord <- bentos_coord + geom_point(data=df2,aes (x= Lon, y=Lat),
                                            shape=19, size=0.1, 
                                            position = jitter)
  
  peixes_coord ## mostre o mapa
}


###  kernel function from Cooke et al. 2019 (check it out here : https://github.com/03rcooke/hyper_pca/commit/2b7df79a30242d3d479e75382a8865df3f5a6f7d)

cl <- function(df, prob) {
  dx <- diff(df$x[1:2])
  dy <- diff(df$y[1:2])
  sz <- sort(df$z)
  c1 <- cumsum(sz) * dx * dy
  approx(c1, sz, xout = 1 - prob)$y
}

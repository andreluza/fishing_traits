

# ===========================================
#                                  
#             FISHING TRAITS

#               ><(((º> 

# ===========================================

# load kernel density function cl
source("R/packages.R")
source("R/functions_kernel.R")

# create directory to host resutls
dir.create ("output")

# trait data
# load your data
data <- read.csv(here ("data","Table_R_oficial.csv"), header=TRUE, sep = ",",
                 fileEncoding="latin1")

# Community
comunidade <- data[,c("Species","Sambaquis","Colonial","Present","All_periods")]
rownames(comunidade) <- comunidade[,"Species"]
comunidade <- comunidade[,which (colnames(comunidade) != "Species")]

#Traits
all_traits_spp <- data[,c(1,7:9, 12:14)]
rownames(all_traits_spp) <- all_traits_spp[,"Species"]
all_traits_spp <- all_traits_spp[,which (colnames(all_traits_spp) != "Species")]
all_traits_spp$Max_size_cm <- as.numeric(all_traits_spp$Max_size_cm)

# bind size based on fishbase
all_traits_spp[which(is.na(all_traits_spp$Max_size_cm)),"Max_size_cm"] <-  16.3

### standardizing maximum size
all_traits_spp <- cbind (all_traits_spp,
                         max_size_std = log(all_traits_spp [,"Max_size_cm"]))   # using log  


# average body size
mean(as.numeric(data$Max_size_cm),na.rm=T)
sd(as.numeric(data$Max_size_cm),na.rm=T)
# per period
mean(as.numeric(data$Max_size_cm[which(data$Sambaquis ==1)]),na.rm=T) #sambaquis
sd(as.numeric(data$Max_size_cm[which(data$Sambaquis ==1)]),na.rm=T)
mean(as.numeric(data$Max_size_cm[which(data$Colonial ==1)]),na.rm=T) #colonial
sd(as.numeric(data$Max_size_cm[which(data$Colonial ==1)]),na.rm=T)
mean(as.numeric(data$Max_size_cm[which(data$Present ==1)]),na.rm=T) #modern
sd(as.numeric(data$Max_size_cm[which(data$Present ==1)]),na.rm=T)

### ordering traits
all_traits_spp$Vertical_position <- ordered (all_traits_spp$Vertical_position)
all_traits_spp$Habitat_use.1 <- ordered (all_traits_spp$Habitat_use.1)
all_traits_spp$Mobility <- ordered (all_traits_spp$Mobility)
all_traits_spp$Group_size <- ordered (all_traits_spp$Group_size)

# gawdis  
# help here: https://cran.r-project.org/web/packages/gawdis/vignettes/gawdis.html
gawdist.traits <- gawdis(all_traits_spp [,which (colnames(all_traits_spp ) %in% 
                                                   c("Max_size_cm","center", "scale") != T)],
                         w.type ="analytic")  #GAWDIS


#----------------------------------

# PCoA
## principal component analysis, and extract the two first scores
#  transformation to make the gower matrix as euclidean. nf= number of axis
gawdis.euclid <- quasieuclid(gawdist.traits)

# the complete trait space
PCoA <- dudi.pco (gawdis.euclid, scannf=F, nf=10)

# percentage of inertia explained by the two first axes
(Inertia2<-(PCoA$eig[1]+PCoA$eig[2]+PCoA$eig[3]) /(sum(PCoA$eig)))

## only the frst axis
Inertia.first <- (PCoA$eig[1]) /(sum(PCoA$eig))
## only the frst axis
Inertia.scnd <- (PCoA$eig[2]) /(sum(PCoA$eig))
## only the frst axis
Inertia.trd <- (PCoA$eig[3]) /(sum(PCoA$eig))

# axes explication 
exp_axis <- c(Inertia.first,Inertia.scnd , Inertia.trd)

# extracting scores (coordinates)
PCoA_scores <- PCoA$li

## bind taxon and species name
PCoA_scores <- data.frame (PCoA$li[,1:3], 
                           Species = rownames(all_traits_spp))

##Extracting scores per period (trait space slices)
#Sambaquis
scores_sambaquis<- PCoA_scores[which(rownames(PCoA_scores) %in% rownames(comunidade)[which(comunidade$Sambaquis>0)]),]
scores_sambaquis$Period <- "Sambaquis"

#Colonial
scores_colonial<- PCoA_scores[which(rownames(PCoA_scores) %in% rownames(comunidade)[which(comunidade$Colonial>0)]),]
scores_colonial$Period <- "Colonial"

#Current
scores_current<- PCoA_scores[which(rownames(PCoA_scores) %in% rownames(comunidade)[which(comunidade$Present>0)]),]
scores_current$Period <- "Current"

# All periods
scores_all_periods<- PCoA_scores[which(rownames(PCoA_scores) %in% rownames(comunidade)[which(comunidade$All_periods>0)]),]
scores_all_periods$Period <- "All periods"


# PCoA Scores All Periods
PCoA_scores_all <- rbind(scores_sambaquis, 
                         scores_colonial, 
                         scores_current, 
                         scores_all_periods)



## transforming SCORES to the long format
PCoA_scores_long <- melt(PCoA_scores_all, id.vars=c("Period", "Species"))

# and then to the wide format
PCoA_scores_wide <- dcast(PCoA_scores_long, Period+Species ~ variable, value="value")

## extract the loadings (correction of traits with each axis)
quantitative_traits<- all_traits_spp
quantitative_traits$Vertical_position <- as.numeric(quantitative_traits$Vertical_position ) #Vertical position
quantitative_traits$Habitat_use.1 <- as.numeric(quantitative_traits$Habitat_use.1 ) # Habitat use
quantitative_traits$Mobility <- as.numeric(quantitative_traits$Mobility )#Mobility
quantitative_traits$Group_size <- as.numeric(quantitative_traits$Group_size )#Group size

#Trophic category
quantitative_traits$Trophic_category [which(quantitative_traits$Trophic_category == "PLANK")] <- 1
quantitative_traits$Trophic_category [which(quantitative_traits$Trophic_category == "HERB")] <- 2
quantitative_traits$Trophic_category [which(quantitative_traits$Trophic_category == "INV")] <- 3
quantitative_traits$Trophic_category [which(quantitative_traits$Trophic_category == "PLANK-PISC")] <- 4
quantitative_traits$Trophic_category [which(quantitative_traits$Trophic_category == "OMNI")] <- 5
quantitative_traits$Trophic_category [which(quantitative_traits$Trophic_category == "PISC")] <- 6
quantitative_traits$Trophic_category [which(quantitative_traits$Trophic_category == "MCAR")] <- 7
quantitative_traits$Trophic_category <- as.numeric(quantitative_traits$Trophic_category ) #as.numeric
# apply (quantitative_traits,2,class)   #para ver a classe de todos os traits


# test the correlation between traits
cor(quantitative_traits,use = "complete.obs")
PCoA_loadings <- lapply (seq (1,ncol(PCoA$li)), function (i) 
  
  
        cor (quantitative_traits[,which (colnames(quantitative_traits) %in% 
                                           c("Max_size_cm") != T)], 
             PCoA$li[,i], use = "complete.obs") 
  
        )    #correlation matrix each trait


PCoA_loadings <- do.call(cbind,PCoA_loadings)    # melt    
colnames(PCoA_loadings)<- colnames(PCoA$li)

# fortify (long format)    
# show the correlation of each trait with each axis
PCoA_loadings_wide <- fortify(data.frame (PCoA_loadings))

# bind taxon to data
PCoA_scores_wide$taxon <- data [match (PCoA_scores_wide$Species,data$Species), "Type"]

# scalar to adjust arrow length
#sc_arrow <- 0.5
#PCoA_loadings_wide <- PCoA_loadings_wide*sc_arrow ## ajusting

# kernel density estimation for each period 
# apply to all periods
fator_de_correcao_arrow<-0.5 # fator pra corrigir o comprimento da seta
PCoA_scores_wide_save<-PCoA_scores_wide #  * fator_de_correcao_arrow

# run kernel and produce plots
kde_period <- lapply (c("Sambaquis", "Colonial", "Current","All periods"), function (i) {
  
  # subsetting 
  PCoA_scores_wide <- PCoA_scores_wide_save[which(PCoA_scores_wide_save$Period == i),]
  
  # optimal bandwidth estimation
  hpi_mi_d1 <- Hpi(x = PCoA_scores_wide[,c("A1","A2")])
  
  # kernel density estimation
  est_mi_d1 <- kde(x = PCoA_scores_wide[,c("A1","A2")], 
                   H = hpi_mi_d1, 
                   compute.cont = TRUE)  
  
  # bandwidths for each point
  den_mi_d1 <- list(est_mi_d1$eval.points[[1]], est_mi_d1$eval.points[[2]], 
                    est_mi_d1$estimate)
  names(den_mi_d1) <- c("x", "y", "z")
  dimnames(den_mi_d1$z) <- list(den_mi_d1$x, den_mi_d1$y)
  dcc_mi_d1 <- melt(den_mi_d1$z)
  
  # 0.5 probability kernel
  # run kernel
  cl_50_mi_d1 <- cl(df = den_mi_d1, prob = 0.50)
  # 0.95 probability kernel
  cl_95_mi_d1 <- cl(df = den_mi_d1, prob = 0.95)
  # 0.99 probability kernel
  cl_99_mi_d1 <- cl(df = den_mi_d1, prob = 0.99)
  
  ## PCoA
  # colour palette
  col_pal <- colorRampPalette(c("red", "yellow", "white"))(100)
  
  # plot PCoA
  PCoA_plot_mi_d1 <- ggplot(dcc_mi_d1, aes(x = Var1, y = Var2)) + 
    
    # coloured probabilty background
    geom_raster(aes(fill = value)) +
    scale_fill_gradientn(colours = rev(col_pal), 
                         limits = c(0,11)) +
    
    # points for species
    geom_point(data = PCoA_scores_wide, 
               aes(x = A1,y=A2), size = 0.3, alpha = 0.5, colour = "grey20") +
    
    # probability kernels
    geom_contour(aes(z = value), breaks = cl_50_mi_d1, colour = "grey30", size = 1) +
    geom_contour(aes(z = value), breaks = cl_95_mi_d1, colour = "grey60", size = 1) +
    geom_contour(aes(z = value), breaks = cl_99_mi_d1, colour = "grey70", size = 1) +
    #scale_fill_gradient(low = 'white', high = 'red') +
    coord_equal() +
    theme_classic()+
    
    # add arrows
    geom_segment(data = PCoA_loadings_wide[,c("A1", "A2")]*fator_de_correcao_arrow, aes(x = 0, y = 0, 
                                                                                        xend = A1, 
                                                                                        yend = A2), 
                 arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
    
    # add dashed arrows ends
    geom_segment(data = PCoA_loadings_wide[,c("A1", "A2")]*fator_de_correcao_arrow, aes(x = 0, y = 0, 
                                                                                        xend = -A1, 
                                                                                        yend = -A2), 
                 arrow = arrow(length = unit(0.2, "cm")), lty = 5, colour = "darkgrey") +
    # add arrow labels
    geom_text(data = PCoA_loadings_wide*fator_de_correcao_arrow, aes(x = A1, y = A2, 
                                                                     label = rownames(PCoA_loadings_wide)),
              size = 4, 
              nudge_x = c(0, 0, 0, 0, 0), 
              nudge_y = c(0, 0, 0,0,0)) +
    
    # axis labels - see comp_var
    labs(x = paste ("A1 (",round(exp_axis[1],2)*100,"%)",sep=""), 
         y = paste ("A2 (",round(exp_axis[2],2)*100,"%)",sep="")) +
    # edit plot
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(colour = "black",size =5),
          axis.title = element_text(colour = "black",size =7),
          legend.position = "none"#,
          #text = element_text(size = 20)
    )
  PCoA_plot_mi_d1
  
  ## adding the name of some species to check if this analysis makes sense
  ## affinity between species and axes (based on score value)
  range_of_values <- apply (PCoA_scores_wide [,c("A1", "A2")],2,range)
  # select
  species_to_project <- rbind (
    PCoA_scores_wide [which(PCoA_scores_wide [,"A1"] == range_of_values[1]),],
    PCoA_scores_wide [which(PCoA_scores_wide [,"A1"] == range_of_values[2]),],
    PCoA_scores_wide [which(PCoA_scores_wide [,"A2"] == range_of_values[3]),],
    PCoA_scores_wide [which(PCoA_scores_wide [,"A2"] == range_of_values[4]),],
    
    # species close to ordination origin
    # PCoA_scores_wide [which(abs(PCoA_scores_wide [,"A1"] - 0) == min(abs(PCoA_scores_wide [,"A1"] - 0))),],
    # PCoA_scores_wide [which(abs(PCoA_scores_wide [,"A2"] - 0) == min(abs(PCoA_scores_wide [,"A2"] - 0))),],
    
    # locations with the highest densities of species
    PCoA_scores_wide [which(round (PCoA_scores_wide [,"A1"],2) == round(dcc_mi_d1[which(dcc_mi_d1$value== max(dcc_mi_d1$value)),"Var1"],2)),]
    
    
  )
  
  PCoA_plot_mi_d1 <- PCoA_plot_mi_d1 + 
    geom_text(data = species_to_project, aes(x = A1, y = A2, 
                                             label = Species),
              size = 3, 
              nudge_x = c(0, 0, 0, 0, 0,0), 
              nudge_y = c(0, 0, 0,0,0,0)) + 
    xlim(c(-0.9, 0.9)) + ylim (c(-0.8, 0.8)) + 
    geom_point(data = species_to_project, aes(x = A1, y = A2),size=2)
    
  # all species projected in the plot
  plot_with_labels <-  PCoA_plot_mi_d1 + 
    
    geom_text(data = PCoA_scores_wide, aes(x = A1, y = A2, 
                                           label = Species),
              size = 2, 
              nudge_x = c(0, 0, 0, 0, 0,0), 
              nudge_y = c(0, 0, 0,0,0,0)) + 
    xlim(c(-0.9, 0.9)) + ylim (c(-0.8, 0.8))
  
  
  # plot of group position
  group_position <- ggplot (PCoA_scores_wide,
                            aes (x= A1,y=A2,
                                 colour=taxon,
                                 label = Species)) + 
    geom_point() + 
    geom_text(size=2)
  
  
  # list of results
  res <- list (species_to_project = species_to_project,
               plot = PCoA_plot_mi_d1,
               plot_with_labels = plot_with_labels,
               group_position = group_position,
               density50 = cl_50_mi_d1,
               density95 = cl_95_mi_d1,
               density99 = cl_99_mi_d1,
               hpi_mi_d1 = hpi_mi_d1)
  
  ; # return
  # display plot
  res
})


###############################
### Plots

# plot one particular period
kde_period[[1]]$plot +  xlab ("Axis I") + ylab ("Axis II")+ theme(legend.position = "top")
kde_period[[2]]$plot +  xlab ("Axis I") + ylab ("Axis II")+ theme(legend.position = "top")
kde_period[[3]]$plot +  xlab ("Axis I") + ylab ("Axis II")+ theme(legend.position = "top")
kde_period[[4]]$plot +  xlab ("Axis I") + ylab ("Axis II")+ theme(legend.position = "top")

# find the plot with all species
kde_period[[1]]$plot_with_labels

# and taxon position
kde_period[[1]]$group_position

# arrange plots
grid.arrange(kde_period[[1]]$plot +  xlab ("Axis I") + ylab ("Axis II")+ theme(legend.position = "top"),
             kde_period[[2]]$plot +  xlab ("Axis I") + ylab ("Axis II")+ theme(legend.position = "top"),
             kde_period[[3]]$plot +  xlab ("Axis I") + ylab ("Axis II")+ theme(legend.position = "top"),
             kde_period[[4]]$plot +  xlab ("Axis I") + ylab ("Axis II")+ theme(legend.position = "top"), 
             ncol=4,nrow=1) # 4 plots in the row


#-----------------------------------

# the complete trait space 
# run kernel and produce plots
# subsetting 
PCoA_scores_whole <- PCoA_scores_wide_save

# optimal bandwidth estimation
hpi_mi_d1 <- Hpi(x = PCoA_scores_whole[,c("A1","A2")])

# kernel density estimation
est_mi_d1 <- kde(x = PCoA_scores_whole[,c("A1","A2")], 
                 H = hpi_mi_d1, 
                 compute.cont = TRUE)  

# bandwidths for each point
den_mi_d1 <- list(est_mi_d1$eval.points[[1]], est_mi_d1$eval.points[[2]], 
                  est_mi_d1$estimate)
names(den_mi_d1) <- c("x", "y", "z")
dimnames(den_mi_d1$z) <- list(den_mi_d1$x, den_mi_d1$y)
dcc_mi_d1_whole <- melt(den_mi_d1$z)


# ---------------------------------------

# Map the difference between periods
# first period : sambaquis
# run kernel and produce plots
# subsetting 
PCoA_scores_wide <- PCoA_scores_wide_save[which(PCoA_scores_wide_save$Period == "Sambaquis"),]

# optimal bandwidth estimation
hpi_mi_d1 <- Hpi(x = PCoA_scores_wide[,c("A1","A2")])

# kernel density estimation
est_mi_d1 <- kde(x = PCoA_scores_wide[,c("A1","A2")], 
                 H = hpi_mi_d1, 
                 compute.cont = TRUE)  

# bandwidths for each point
den_mi_d1 <- list(est_mi_d1$eval.points[[1]], est_mi_d1$eval.points[[2]], 
                  est_mi_d1$estimate)
names(den_mi_d1) <- c("x", "y", "z")
dimnames(den_mi_d1$z) <- list(den_mi_d1$x, den_mi_d1$y)
dcc_mi_d1 <- melt(den_mi_d1$z)

## PCoA
# colour palette
col_pal <- colorRampPalette(c("red", "yellow", "white"))(100)

# plot PCoA
PCoA_plot_mi_d1_sambaquis <- ggplot(dcc_mi_d1, aes(x = Var1, y = Var2)) + 
  
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal),
                       limits=c(0,12)) +
  
  # points for species
  geom_point(data = PCoA_scores_wide, 
             aes(x = A1,y=A2), size = 0.3, alpha = 0.5, colour = "grey20") +
  
  coord_equal() +
  theme_classic()+
  
  # add arrows
  geom_segment(data = PCoA_loadings_wide[,c("A1", "A2")]*fator_de_correcao_arrow, aes(x = 0, y = 0, 
                                                                                      xend = A1, 
                                                                                      yend = A2), 
               arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  
  # add dashed arrows ends
  geom_segment(data = PCoA_loadings_wide[,c("A1", "A2")]*fator_de_correcao_arrow, aes(x = 0, y = 0, 
                                                                                      xend = -A1, 
                                                                                      yend = -A2), 
               arrow = arrow(length = unit(0.2, "cm")), lty = 5, colour = "darkgrey") +
  # add arrow labels
  geom_text(data = PCoA_loadings_wide*fator_de_correcao_arrow, aes(x = A1, y = A2, 
                                                                   label = rownames(PCoA_loadings_wide)),
            size = 3.3, 
            nudge_x = c(0, 0, 0, 0, 0), 
            nudge_y = c(0, 0, 0,0,0)) +
  xlim(c(-0.9, 0.9)) + ylim (c(-0.8, 0.8)) + 
  # axis labels - see comp_var
  labs(x = paste ("A1 (",round(exp_axis[1],2)*100,"%)",sep=""), 
       y = paste ("A2 (",round(exp_axis[2],2)*100,"%)",sep="")) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(colour = "black",size =5),
        axis.title = element_text(colour = "black",size =7),
        legend.position = "top"#,
        #text = element_text(size = 20)
  )


# second period: colonial
# subsetting 
PCoA_scores_wide_c <- PCoA_scores_wide_save[which(PCoA_scores_wide_save$Period == "Colonial"),]

# optimal bandwidth estimation
hpi_mi_d1_c <- Hpi(x = PCoA_scores_wide_c[,c("A1","A2")])

# kernel density estimation
est_mi_d1_c <- kde(x = PCoA_scores_wide_c[,c("A1","A2")], 
                   H = hpi_mi_d1_c, 
                   compute.cont = TRUE)  


# bandwidths for each point
den_mi_d1_c <- list(est_mi_d1_c$eval.points[[1]], est_mi_d1_c$eval.points[[2]], 
                    est_mi_d1_c$estimate)
names(den_mi_d1_c) <- c("x", "y", "z")
dimnames(den_mi_d1_c$z) <- list(den_mi_d1_c$x, den_mi_d1_c$y)
dcc_mi_d1_c <- melt(den_mi_d1_c$z)

# plot PCoA
PCoA_plot_mi_d1_colonial <- ggplot(dcc_mi_d1_c, aes(x = Var1, y = Var2)) + 
  
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal),
                       limits=c(0,12)) +
  
  # points for species
  geom_point(data = PCoA_scores_wide, 
             aes(x = A1,y=A2), size = 0.3, alpha = 0.5, colour = "grey20") +
  
  coord_equal() +
  theme_classic()+
  xlim(c(-0.9, 0.9)) + ylim (c(-0.8, 0.8)) + 
  # add arrows
  geom_segment(data = PCoA_loadings_wide[,c("A1", "A2")]*fator_de_correcao_arrow, aes(x = 0, y = 0, 
                                                                                      xend = A1, 
                                                                                      yend = A2), 
               arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  
  # add dashed arrows ends
  geom_segment(data = PCoA_loadings_wide[,c("A1", "A2")]*fator_de_correcao_arrow, aes(x = 0, y = 0, 
                                                                                      xend = -A1, 
                                                                                      yend = -A2), 
               arrow = arrow(length = unit(0.2, "cm")), lty = 5, colour = "darkgrey") +
  
  # axis labels - see comp_var
  labs(x = paste ("A1 (",round(exp_axis[1],2)*100,"%)",sep=""), 
       y = paste ("A2 (",round(exp_axis[2],2)*100,"%)",sep="")) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(colour = "black",size =5),
        axis.title = element_text(colour = "black",size =7),
        legend.position = "top"
  )



# third period: current
# subsetting 
PCoA_scores_wide_b <- PCoA_scores_wide_save[which(PCoA_scores_wide_save$Period == "Current"),]


# optimal bandwidth estimation
hpi_mi_d1_b <- Hpi(x = PCoA_scores_wide_b[,c("A1","A2")])

# kernel density estimation
est_mi_d1_b <- kde(x = PCoA_scores_wide_b[,c("A1","A2")], 
                   H = hpi_mi_d1_b, 
                   compute.cont = TRUE)  


# bandwidths for each point
den_mi_d1_b <- list(est_mi_d1_b$eval.points[[1]], est_mi_d1_b$eval.points[[2]], 
                    est_mi_d1_b$estimate)
names(den_mi_d1_b) <- c("x", "y", "z")
dimnames(den_mi_d1_b$z) <- list(den_mi_d1_b$x, den_mi_d1_b$y)
dcc_mi_d1_b <- melt(den_mi_d1_b$z)

# plot PCoA
PCoA_plot_mi_d1_current <- ggplot(dcc_mi_d1_b, aes(x = Var1, y = Var2)) + 
  
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rev(col_pal),
                       limits=c(0,12)) +
  
  # points for species
  geom_point(data = PCoA_scores_wide, 
             aes(x = A1,y=A2), size = 0.3, alpha = 0.5, colour = "grey20") +
  
  coord_equal() +
  theme_classic()+
  xlim(c(-0.9, 0.9)) + ylim (c(-0.8, 0.8)) + 
  # add arrows
  geom_segment(data = PCoA_loadings_wide[,c("A1", "A2")]*fator_de_correcao_arrow, aes(x = 0, y = 0, 
                                                                                      xend = A1, 
                                                                                      yend = A2), 
               arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  
  # add dashed arrows ends
  geom_segment(data = PCoA_loadings_wide[,c("A1", "A2")]*fator_de_correcao_arrow, aes(x = 0, y = 0, 
                                                                                      xend = -A1, 
                                                                                      yend = -A2), 
               arrow = arrow(length = unit(0.2, "cm")), lty = 5, colour = "darkgrey") +
 
  # axis labels - see comp_var
  labs(x = paste ("A1 (",round(exp_axis[1],2)*100,"%)",sep=""), 
       y = paste ("A2 (",round(exp_axis[2],2)*100,"%)",sep="")) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(colour = "black",size =5),
        axis.title = element_text(colour = "black",size =7),
        legend.position = "top"
  )


# difference between sambquis and colonial
diff_kernels <- den_mi_d1_c$z - den_mi_d1$z # difference
#diff_kernels <- abs(diff_kernels) # absolute diff
diff_kernels_coord <- den_mi_d1
diff_kernels_coord$z <- diff_kernels

# melt to plot
dcc_mi_d1_test <- melt(diff_kernels_coord$z)

# plot PCoA
PCoA_plot_mi_d1_test_sambaqui_colonial <- ggplot(dcc_mi_d1_test, aes(x = Var1, y = Var2)) + 
  
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_colour_gradient2(
    low = ("#F46B2F"),
    mid = "white",
    high = ("#4268cb"),
    midpoint = 0,
    aesthetics = "fill",
    limits=c(-5,7)
  ) +

    # points for species
  geom_point(data = PCoA_scores_wide, 
             aes(x = A1,y=A2), size = 0.3, alpha = 0.5, colour = "grey20") +
  
  coord_equal() +
  theme_classic()+
  xlim(c(-0.9, 0.9)) + ylim (c(-0.7, 0.7)) + 
  
  # add arrows
  geom_segment(data = PCoA_loadings_wide[,c("A1", "A2")]*fator_de_correcao_arrow, aes(x = 0, y = 0, 
                                                                                      xend = A1, 
                                                                                      yend = A2), 
               arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  
  # add dashed arrows ends
  geom_segment(data = PCoA_loadings_wide[,c("A1", "A2")]*fator_de_correcao_arrow, aes(x = 0, y = 0, 
                                                                                      xend = -A1, 
                                                                                      yend = -A2), 
               arrow = arrow(length = unit(0.2, "cm")), lty = 5, colour = "darkgrey") +
  
  # axis labels - see comp_var
  labs(x = paste ("A1 (",round(exp_axis[1],2)*100,"%)",sep=""), 
       y = paste ("A2 (",round(exp_axis[2],2)*100,"%)",sep="")) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(colour = "black",size =5),
        axis.title = element_text(colour = "black",size =7),
        legend.position = "top"
  )


# difference between colonial and current

diff_kernels <- den_mi_d1_b$z - den_mi_d1_c$z # difference
#diff_kernels <- abs(diff_kernels) # absolute diff
diff_kernels_coord <- den_mi_d1
diff_kernels_coord$z <- diff_kernels

# melt to plot
dcc_mi_d1_test2 <- melt(diff_kernels_coord$z)

# plot PCoA
PCoA_plot_mi_d1_colonial_current <- ggplot(dcc_mi_d1_test2, aes(x = Var1, y = Var2)) + 
  
  # coloured probabilty background
  geom_raster(aes(fill = value)) +
  scale_colour_gradient2(
    low = ("#F46B2F"),
    mid = "white",
    high = ("#4268cb"),
    midpoint = 0,
    aesthetics = "fill",
    limits=c(-5,7)
  ) +
  
  # points for species
  geom_point(data = PCoA_scores_wide, 
             aes(x = A1,y=A2), size = 0.3, alpha = 0.5, colour = "grey20") +
  coord_equal() +
  theme_classic()+
  xlim(c(-0.9, 0.9)) + ylim (c(-0.7, 0.7)) + 
  # add arrows
  geom_segment(data = PCoA_loadings_wide[,c("A1", "A2")]*fator_de_correcao_arrow, aes(x = 0, y = 0, 
                                                                                      xend = A1, 
                                                                                      yend = A2), 
               arrow = arrow(length = unit(0.2, "cm")), colour = "black") +
  
  # add dashed arrows ends
  geom_segment(data = PCoA_loadings_wide[,c("A1", "A2")]*fator_de_correcao_arrow, aes(x = 0, y = 0, 
                                                                                      xend = -A1, 
                                                                                      yend = -A2), 
               arrow = arrow(length = unit(0.2, "cm")), lty = 5, colour = "darkgrey") +
  # axis labels - see comp_var
  labs(x = paste ("A1 (",round(exp_axis[1],2)*100,"%)",sep=""), 
       y = paste ("A2 (",round(exp_axis[2],2)*100,"%)",sep="")) +
  # edit plot
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(colour = "black",size =5),
        axis.title = element_text(colour = "black",size =7),
        legend.position = "top"
  )




# plot of the difference between periods
# difference 
diff_periods <- data.frame (rbind (
                  cbind (val = sum(dcc_mi_d1_test2$value<0)/length(dcc_mi_d1_whole$value), 
                         comp = "Loss",
                         var = "Colonial vs Current"),
                  cbind (val = sum(dcc_mi_d1_test2$value>0)/length(dcc_mi_d1_whole$value), 
                         comp = "Gain",
                         var = "Colonial vs Current"),
                  cbind (val = sum(dcc_mi_d1_test2$value==0)/length(dcc_mi_d1_whole$value), 
                         comp = "No change",
                         var = "Colonial vs Current"),
                  
                  cbind (val = sum(dcc_mi_d1_test$value<0)/length(dcc_mi_d1_whole$value), 
                         comp = "Loss",
                         var = "Sambaquis vs Colonial"),
                  cbind (val = sum(dcc_mi_d1_test$value>0)/length(dcc_mi_d1_whole$value), 
                         comp= "Gain",
                         var = "Sambaquis vs Colonial"),
                  cbind (val = sum(dcc_mi_d1_test$value==0)/length(dcc_mi_d1_whole$value), 
                         comp= "No change",
                         var = "Sambaquis vs Colonial"),
                  
                   cbind (val = sum(dcc_mi_d1$value>0)/sum(dcc_mi_d1_whole$value>0), 
                          comp = "NA",
                          var= "Sambaquis"),
                   cbind (val = sum(dcc_mi_d1_c$value>0)/sum(dcc_mi_d1_whole$value>0), 
                          comp = "NA",
                          var = "Colonial"), 
                   cbind (val = sum(dcc_mi_d1_b$value>0)/sum(dcc_mi_d1_whole$value>0),
                          comp = "NA",
                          var = "Current"))
)
diff_periods$var <- factor(diff_periods$var,
                           levels = c("Sambaquis",
                                      "Colonial",
                                      "Current", 
                                      "Sambaquis vs Colonial",
                                      "Colonial vs Current"))

# plot_diff 
# observed area
plot_area <- ggplot (diff_periods [-grep ("vs", diff_periods$var),], 
                     aes (x= var, y=as.numeric(val))) +
  theme_bw()+
  theme (axis.text = element_text(angle=0, vjust = 0.6,size=8),
         axis.title = element_text(angle=0, vjust = 0.6,size=8)) +
  geom_point(size=3) +
  ylab ("Trait space occupancy\n(proportion of the total area)") +
  xlab ("Periods") +
  coord_flip()

# differece
plot_diff <- ggplot (diff_periods [grep ("vs", diff_periods$var),], 
                     aes (x= var, y=as.numeric(val),
                          fill=comp,colour = comp),size=4) +
  theme_bw()+
  theme (axis.text = element_text(angle=0, vjust = 0.6,size=8),
         axis.title = element_text(angle=0, vjust = 0.6,size=8)) +
  geom_point() +
  geom_point(stroke=2,colour = "black",shape=1)+
  scale_color_manual(values=c("#4268cb","#F46B2F","white"))+
  ylab ("Difference in trait space occupancy") +
  xlab ("Periods") +
  coord_flip()


# save results
pdf (here ("output","plot1.pdf"),width = 12,height=8)

# arrange plot
grid.arrange (
  plot_area,
  kde_period[[1]]$plot + 
    theme(legend.position = "none",
                                  plot.title = element_text(size=8)),
  kde_period[[2]]$plot + theme(legend.position = "none",
                                 plot.title = element_text(size=8)
                               ),
  kde_period[[3]]$plot + theme(legend.position = "right",
                                legend.key.size = unit(0.3, 'cm'), #change legend key size
                                legend.key.height = unit(0.3, 'cm'), #change legend key height
                                legend.key.width = unit(0.3, 'cm'), #change legend key width
                                plot.title = element_text(size=8)
                               ) + 
    guides(fill=guide_legend(title="Density")),
  PCoA_plot_mi_d1_test_sambaqui_colonial+theme(legend.position = "none",
                                               plot.title = element_text(size=8))+ggtitle("Sambaqui vs colonial"),
  PCoA_plot_mi_d1_colonial_current+theme(legend.position = "right",
                                         plot.title = element_text(size=8))+
    ggtitle("Colonial vs current"),
  
  plot_diff,
    ncol=6,
    nrow=6,
    layout_matrix = rbind (c(NA,NA,1,1,NA,NA),
                           c(2,2,3,3,4,4),
                           c(2,2,3,3,4,4),
                           c(NA,5,5,6,6,NA),
                           c(NA,5,5,6,6,NA),
                           c(NA,NA,7,7,7,NA))
)
dev.off()


# end
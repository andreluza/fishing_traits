# load an important function
# closest function

closest<-function(xv,sv){
  xv[which(abs(xv-sv)==min(abs(xv-sv)))]}


# trait data
# load your data
data <- read.csv("Tabela_R_oficial.csv", header=TRUE, sep = ",")


#Comunidade
comunidade <- data[,c("Species","Sambaquis","Colonial","Present","All_periods")]
rownames(comunidade) <- comunidade[,"Species"]
comunidade <- comunidade[,which (colnames(comunidade) != "Species")]


#Traits

all_traits_spp <- data[,c(1,7:9, 12:14)]


#View(all_traits_spp)
rownames(all_traits_spp) <- all_traits_spp[,"Species"]
all_traits_spp <- all_traits_spp[,which (colnames(all_traits_spp) != "Species")]



### padronizando o tamanho maximo de cada spp.

all_traits_spp <- cbind (all_traits_spp,
                         max_size_std = log(all_traits_spp [,"Max_size_cm"]), center = T, scale = T)     #using log 


### classificando vertical position, habitat use, mobility e group size

all_traits_spp$Vertical_position <- ordered (all_traits_spp$Vertical_position)

all_traits_spp$Habitat_use.1 <- ordered (all_traits_spp$Habitat_use.1)

all_traits_spp$Mobility <- ordered (all_traits_spp$Mobility)

all_traits_spp$Group_size <- ordered (all_traits_spp$Group_size)


# gawdis  
install.packages("gawdis")
install.packages("ape")
install.packages("foreach",dependencies=T)
# help here: https://cran.r-project.org/web/packages/gawdis/vignettes/gawdis.html


#gowdist.traits <- gowdis(all_traits_spp [,which (colnames(all_traits_spp ) %in% c("Max_size_cm","center", "scale") != T)],  ord="metric")  # tentativa GOWER
require(gawdis)
gawdist.traits <- gawdis(all_traits_spp [,which (colnames(all_traits_spp ) %in% 
                                                   c("Max_size_cm","center", "scale") != T)], w.type ="analytic")  #GAWDIS

# PCoA
## principal component analysis, and extract the two first scores
require(ape)
#install.packages("reshape")
#install.packages("reshape2")
require(reshape)
require(reshape2)

#  transformation to make the gower matrix as euclidean. nf= number of axis
gawdis.euclid <- quasieuclid(gawdist.traits)


# the complete trait space
PCoA <- dudi.pco (gawdis.euclid, scannf=F, nf=10)


# percentage of inertia explained by the two first axes
(Inertia2<-(PCoA$eig[1]+PCoA$eig[2]+PCoA$eig[3]+PCoA$eig[4]+PCoA$eig[5]+PCoA$eig[6]) /(sum(PCoA$eig)))

# finding axes explanation
axes_inertia <- lapply (seq (2,length(PCoA$eig)), function (i) {
  
  
  Inertia<-sum(PCoA$eig[1:i]) /(sum(PCoA$eig))
  ;Inertia # return
  
  
}
)

# variation explained with the addition of axes
threshold <- 0.9
dat_variation <- data.frame (axis = seq(2,length(PCoA$eig)),
            variation = unlist(axes_inertia))
chosen_axes <- seq (1,
                    dat_variation[which(dat_variation$variation == closest(dat_variation$variation,
                                                       threshold)),"axis"]
)

# plot 
require(ggplot2)
ggplot(dat_variation, aes (x=axis, y = variation)) +
  geom_point()+
  theme_bw()


## explanation of chosen axes
sum(Inertia.first <- (PCoA$eig[chosen_axes]) /(sum(PCoA$eig)))


# building the hypervolume based on the previous axes
require(hypervolume)
hypervolume_test <- hypervolume (PCoA$tab[,1:2], method= "gaussian",
             
             kde.bandwidth=estimate_bandwidth(PCoA$tab[,1:2],
                                              method="silverman"))

# dividing periods
scores_sambaquis<- PCoA$tab[which(rownames(PCoA$tab) %in% 
                                       rownames(comunidade)[which(comunidade$Sambaquis>0)]),]
scores_sambaquis$Period <- "Sambaquis"


# test sambaquis
hypervolume_test_sambaquis <- hypervolume (scores_sambaquis[,1:2], method= "gaussian",
                                 
                                 kde.bandwidth=estimate_bandwidth(scores_sambaquis[,1:2],
                                                                  method="silverman"))


# functional analyzes
# generate the permpath
# Take 200 permutations of the 300 data points. Using more cores is faster.
perm_path <- hypervolume_permute("sambquis_perm_100", 
                                 hypervolume_test, 
                                 hypervolume_test_sambaquis, 
                                 n = 2, 
                                 cores = 8)

# overlap_test 
overlap_test <- hypervolume_overlap_test (hypervolume_test,
                          hypervolume_test_sambaquis)


#PCoA <- cmdscale(gawdist.traits,k=3, eig=TRUE)

## percentage of explanation of PCoA axes
#PCoA $values$Eigenvalues[1:3] / sum(PCoA $values$Eigenvalues)



# extracting scores (coordinates)
PCoA_scores <- PCoA$li

## bind taxon and species name
PCoA_scores <- data.frame (PCoA$li[,1:3], 
                           Species = rownames(all_traits_spp))
rownames(PCoA_scores)

##Extracting scores per period (trait space slices)
#Sambaquis
scores_sambaquis<- PCoA_scores[which(rownames(PCoA_scores) %in% rownames(comunidade)[which(comunidade$Sambaquis>0)]),]
scores_sambaquis$Period <- "Sambaquis"
rownames(scores_sambaquis)
table (comunidade$Sambaquis>0)
View(scores_sambaquis)

#Colonial
scores_colonial<- PCoA_scores[which(rownames(PCoA_scores) %in% rownames(comunidade)[which(comunidade$Colonial>0)]),]
scores_colonial$Period <- "Colonial"
View(scores_colonial)

#Current
scores_current<- PCoA_scores[which(rownames(PCoA_scores) %in% rownames(comunidade)[which(comunidade$Present>0)]),]
scores_current$Period <- "Current"
View(scores_current)

# All periods
scores_all_periods<- PCoA_scores[which(rownames(PCoA_scores) %in% rownames(comunidade)[which(comunidade$All_periods>0)]),]
scores_all_periods$Period <- "All periods"
View(scores_all_periods)



# PCoA Scores All Periods
PCoA_scores_all <- rbind(scores_sambaquis, 
                         scores_colonial, 
                         scores_current, 
                         scores_all_periods)
View(PCoA_scores_all)




## transforming SCORES to the long format
PCoA_scores_long <- melt(PCoA_scores_all, id.vars=c("Period", "Species"))
View(PCoA_scores_long)

# and then to the wide format
PCoA_scores_wide <- dcast(PCoA_scores_long, Period+Species ~ variable, value="value")
View(PCoA_scores_wide)

## extract the loadings (correction of traits with each axis)

quantitative_traits<- all_traits_spp
View(all_traits_spp)

#Vertical position

quantitative_traits$Vertical_position <- as.numeric(quantitative_traits$Vertical_position )

# Habitat use
quantitative_traits$Habitat_use.1 <- as.numeric(quantitative_traits$Habitat_use.1 )

#Mobility
quantitative_traits$Mobility <- as.numeric(quantitative_traits$Mobility )

#Group size
quantitative_traits$Group_size <- as.numeric(quantitative_traits$Group_size )

#Trophic category
quantitative_traits$Trophic_category [which(quantitative_traits$Trophic_category == "PLANK")] <- 1
quantitative_traits$Trophic_category [which(quantitative_traits$Trophic_category == "HERB")] <- 2
quantitative_traits$Trophic_category [which(quantitative_traits$Trophic_category == "INV")] <- 3
quantitative_traits$Trophic_category [which(quantitative_traits$Trophic_category == "PLANK-PISC")] <- 4
quantitative_traits$Trophic_category [which(quantitative_traits$Trophic_category == "OMNI")] <- 5
quantitative_traits$Trophic_category [which(quantitative_traits$Trophic_category == "PISC")] <- 6
quantitative_traits$Trophic_category [which(quantitative_traits$Trophic_category == "MCAR")] <- 7

quantitative_traits$Trophic_category <- as.numeric(quantitative_traits$Trophic_category ) #as.numeric

apply (quantitative_traits,2,class)   #para ver a classe de todos os traits



# test  the correlation between traits

cor(quantitative_traits,use = "complete.obs")



### 

PCoA_loadings <- lapply (seq (1,ncol(PCoA$li)), function (i) 
  
  
        cor (quantitative_traits[,which (colnames(quantitative_traits) %in% 
                                           c("Max_size_cm" , "center", "scale") != T)], 
             PCoA$li[,i], use = "complete.obs") 
  
        )    #correlation matrix each trait


PCoA_loadings <- do.call(cbind,PCoA_loadings)    # melt    
colnames(PCoA_loadings)<- colnames(PCoA$li)



require(ggplot2)
# fortify (long format)    
# show the correlation of each trait with each axis
PCoA_loadings_wide <- fortify(data.frame (PCoA_loadings))


# bind taxon to data
PCoA_scores_wide$taxon <- data [match (PCoA_scores_wide$Species,data$Species), "Type"]


# scalar to adjust arrow length
#sc_arrow <- 0.5
## ajusting
#PCoA_loadings_wide <- PCoA_loadings_wide*sc_arrow

# kernel density estimation for each period 
require(ks);require(hypervolume);require(vegan);require(reshape)
View(all_traits_spp)
# load kernel density function cl
source("functions_kernel.R")

# apply to all periods

fator_de_correcao_arrow<-0.5 # fator pra corrigir o comprimento da seta
PCoA_scores_wide_save<-PCoA_scores_wide



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
    labs(x = paste ("A1 (",round(exp_axis[1],2),"%)",sep=""), 
         y = paste ("A2 (",round(exp_axis[2],2),"%)",sep="")) +
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


# 12 warning messages2

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


library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)


grid.arrange(kde_period[[1]]$plot +  xlab ("Axis I") + ylab ("Axis II")+ theme(legend.position = "top"),
             kde_period[[2]]$plot +  xlab ("Axis I") + ylab ("Axis II")+ theme(legend.position = "top"),
             kde_period[[3]]$plot +  xlab ("Axis I") + ylab ("Axis II")+ theme(legend.position = "top"),
             kde_period[[4]]$plot +  xlab ("Axis I") + ylab ("Axis II")+ theme(legend.position = "top"), 
             ncol=4,nrow=1) # 4 plots na linha


grid.arrange(kde_period[[1]]$plot +  xlab ("Axis I") + ylab ("Axis II")+ theme(legend.position = "top"),
             kde_period[[2]]$plot +  xlab ("Axis I") + ylab ("Axis II")+ theme(legend.position = "top"),
             kde_period[[3]]$plot +  xlab ("Axis I") + ylab ("Axis II")+ theme(legend.position = "top"),
             kde_period[[4]]$plot +  xlab ("Axis I") + ylab ("Axis II")+ theme(legend.position = "top"), 
             ncol=1,nrow=4) # 4 plots na coluna


######################################

# Map the difference


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
  
  
  # second period
  # subsetting 
  PCoA_scores_wide_b <- PCoA_scores_wide_save[which(PCoA_scores_wide_save$Period == k),]
  
  
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
  
  # 0.5 probability kernel
  # run kernel
  cl_50_mi_d1_b <- cl(df = den_mi_d1_b, prob = 0.50)
  # 0.95 probability kernel
  cl_95_mi_d1_b <- cl(df = den_mi_d1_b, prob = 0.95)
  # 0.99 probability kernel
  cl_99_mi_d1_b <- cl(df = den_mi_d1_b, prob = 0.99)
  
  
  
  # difference
# test <- dcc_mi_d1_b$value - dcc_mi_d1$value
# test <- cbind (dcc_mi_d1_b,
#                value=test)[,-3]
# 
# # 0.5 probability kernel
# # run kernel
# cl_50_mi_d1_t <- cl(df = test, prob = 0.50)
# # 0.95 probability kernel
# cl_95_mi_d1_t <- cl(df = test, prob = 0.95)
# # 0.99 probability kernel
# cl_99_mi_d1_t <- cl(df = test, prob = 0.99)
  
  
  
  
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
    geom_contour(aes(z = value), breaks = cl_50_mi_d1_b, colour = "grey30", size = 1) +
    geom_contour(aes(z = value), breaks = cl_95_mi_d1_b, colour = "grey60", size = 1) +
    geom_contour(aes(z = value), breaks = cl_99_mi_d1_b, colour = "grey70", size = 1) +
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
    labs(x = paste ("A1 (",round(exp_axis[1],2),"%)",sep=""), 
         y = paste ("A2 (",round(exp_axis[2],2),"%)",sep="")) +
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






library(ggplot2)

# function to calculate E, V, B fields
# returns dataframe
# Args: vectors <source1> and <source2> containing x-y coordinates (fractions of plot width/height) 
# and field strength (equal to k*q for electrostatic fields or µ0*I/2/pi for magnetic fields).
# Usage example:
# df_field <- calc_fields(c('x'=0.33, 'y'=0.5, 's'=-0.3), c('x'=0.66, 'y'=0.5, 's'=1))

calc_fields <- function (source1, source2) {
  # Vector normalization scalar (makes field lines longer or shorter). Recommended value: 0.012
  norm_scale <- 0.012
  
  # Resolution of coordinate grid. Recommended value: c('x'=900, 'y'=600)
  res <- c('x'=900, 'y'=600)
  
  # Physical distance scale (plot width in m)
  dist_scale <- 9
  
  source1['x'] <- source1['x']*dist_scale
  source2['x'] <- source2['x']*dist_scale
  source1['y'] <- source1['y']*dist_scale/res['x']*res['y']
  source2['y'] <- source2['y']*dist_scale/res['x']*res['y']
  
  fun_Ex <- function(x, y, X, Y, k)       {  k*(x-X)/((x-X)^2+(y-Y)^2)^3 }
  fun_Ey <- function(x, y, X, Y, k)       {  k*(y-Y)/((x-X)^2+(y-Y)^2)^3 }
  fun_mag <- function(Field_x, Field_y)   { sqrt(Field_x^2 + Field_y^2) }
  fun_log_mag <- function(Field_x, Field_y) { log(sqrt(Field_x^2 + Field_y^2)) }
  fun_V <- function(x, y, X, Y, k)        {  k/sqrt(((x-X)^2+(y-Y)^2)) }
  fun_Bx <- function(x, y, X, Y, k)       { k*(Y-y)/((x-X)^2+(y-Y)^2) }
  fun_By <- function(x, y, X, Y, k)       { k*(x-X)/((x-X)^2+(y-Y)^2) }
  
  
  # combine vector elements into dataframe pairs
  df_field <- expand.grid(seq(0, dist_scale, length=res['x']), seq(0, dist_scale/res['x']*res['y'], length=res['y']))
  
  names(df_field) <- c("x", "y")
  
  # Field coordinates and magnitudes
  
  df_field$E1x <- fun_Ex(df_field$x, df_field$y, source1['x'], source1['y'], source1['s'])
  df_field$E1y <- fun_Ey(df_field$x, df_field$y, source1['x'], source1['y'], source1['s'])
  df_field$E1 <- fun_log_mag(df_field$E1x, df_field$E1y)
  
  df_field$E2x <- fun_Ex(df_field$x, df_field$y, source2['x'], source2['y'], source2['s'])
  df_field$E2y <- fun_Ey(df_field$x, df_field$y, source2['x'], source2['y'], source2['s'])
  df_field$E2 <- fun_log_mag(df_field$E2x, df_field$E2y)
  
  df_field$Ex <- df_field$E1x + df_field$E2x
  df_field$Ey <- df_field$E1y + df_field$E2y
  df_field$E <- fun_log_mag(df_field$Ex, df_field$Ey)
  
  df_field$V1 <-fun_V(df_field$x, df_field$y, source1['x'], source1['y'], source1['s'])
  df_field$V2 <-fun_V(df_field$x, df_field$y, source2['x'], source2['y'], source2['s'])
  df_field$V <- df_field$V1 + df_field$V2
  
  df_field$E1normx <- df_field$E1x/fun_mag(df_field$E1, df_field$E1y)
  df_field$E1normy <- df_field$E1y/fun_mag(df_field$E1, df_field$E1y)
  
  df_field$Enormx <- norm_scale*dist_scale*df_field$Ex/fun_mag(df_field$Ex, df_field$Ey)
  df_field$Enormy <- norm_scale*dist_scale*df_field$Ey/fun_mag(df_field$Ex, df_field$Ey)
  
  df_field$B1x <- fun_Bx(df_field$x, df_field$y, source1['x'], source1['y'], source1['s'])
  df_field$B1y <- fun_By(df_field$x, df_field$y, source1['x'], source1['y'], source1['s'])
  df_field$B1 <- fun_mag(df_field$B1x, df_field$B1y)
  
  df_field$B2x <- fun_Bx(df_field$x, df_field$y, source2['x'], source2['y'], source2['s'])
  df_field$B2y <- fun_By(df_field$x, df_field$y, source2['x'], source2['y'], source2['s'])
  df_field$B2 <- fun_mag(df_field$B2x, df_field$B2y)
  
  df_field$Bx <- df_field$B1x + df_field$B2x
  df_field$By <- df_field$B1y + df_field$B2y
  df_field$B <- fun_mag(df_field$Bx, df_field$By)
  
  df_field$B1normx <- norm_scale*dist_scale*df_field$B1x/fun_mag(df_field$B1x, df_field$B1y)
  df_field$B1normy <- norm_scale*dist_scale*df_field$B1y/fun_mag(df_field$B1x, df_field$B1y)
  
  df_field$B2normx <- norm_scale*dist_scale*df_field$B2x/fun_mag(df_field$B2x, df_field$B2y)
  df_field$B2normy <- norm_scale*dist_scale*df_field$B2y/fun_mag(df_field$B2x, df_field$B2y)
  
  df_field$Bnormx <- norm_scale*dist_scale*df_field$Bx/fun_mag(df_field$Bx, df_field$By)
  df_field$Bnormy <- norm_scale*dist_scale*df_field$By/fun_mag(df_field$Bx, df_field$By)
  
  return (df_field)
}



### Plots

df_field <- calc_fields(c('x'=0.33, 'y'=0.5, 's'=1), c('x'=0.66, 'y'=0.5, 's'=1))

# E field plot
ggplot(df_field, aes(x = x, y = y)) + 
  stat_contour(aes(z = E, color=..level..),breaks=-11:8, size = 0.6) +  
  scale_colour_gradient2(name = "log |E|", low = "blue", mid = "green", high = "red",
                         midpoint = 0, space = "Lab", 
                         na.value = "grey50", guide = "colourbar") +
  geom_segment(data = df_field[sample(nrow(df_field), 10000), ], 
               aes(x = x -Enormx/2, y = y - Enormy/2, xend = x + Enormx/2, yend = y + Enormy/2), 
               size = 0.2, alpha = 0.6) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.08, 0.74),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  coord_fixed() + #preserve aspect ratio
  ggtitle("E-field from two charges: q1 = q2")

ggsave('E_1_1.png', height=4.5, width=6)


# E-field and V-field plot
ggplot() + 
  stat_contour(data = df_field, aes(x = x, y = y, z = V, color=..level..), breaks=seq(0,1.5,0.05)^3, size = 0.6) +  
  scale_colour_gradient2(name = "V", low = "blue", mid = "green", high = "red",
                         midpoint = 0, space = "Lab", 
                         na.value = "grey50", guide = "colourbar") +
  geom_segment(data = df_field[sample(nrow(df_field), 10000), ], 
               aes(x = x -Enormx/2, y = y - Enormy/2, xend = x + Enormx/2, yend = y + Enormy/2), 
               size = 0.2, alpha = 0.6) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.08, 0.74),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  coord_fixed() + 
  ggtitle("V-field and E-field of two charges: q1 = q2")

ggsave('V_1_1.png', height=4.5, width=6)


# B-field plot
ggplot() + 
  stat_contour(data = df_field, aes(x = x, y = y, z = B, color=..level..), breaks=seq(0,1.2,0.05)^3, size = 0.6) +  
  scale_colour_gradient2(name = "|B|", low = "blue", mid = "green", high = "red",
                         midpoint = 0, space = "Lab", 
                         na.value = "grey50", guide = "colourbar") +
  geom_segment(data = df_field[sample(nrow(df_field), 8000), ], 
               aes(x = x - Bnormx/2, y = y - Bnormy/2, xend = x + Bnormx/2, yend = y + Bnormy/2), 
               size = 0.2, alpha = 0.6) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.08, 0.74),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  coord_fixed() + 
  ggtitle("B-field of ∞ parallel wires: I1 = I2")
ggsave('B_1_1.png', height=4.5, width=6)






#---- EXTRA PLOTS

df_field <- calc_fields(c('x'=0.33, 'y'=0.5, 's'=-0.3), c('x'=0.66, 'y'=0.5, 's'=1))

# E field plot
ggplot(df_field, aes(x = x, y = y)) + 
  stat_contour(aes(z = E, color=..level..),breaks=-7:8, size = 0.6) +  
  scale_colour_gradient2(name = "log |E|", low = "blue", mid = "green", high = "red",
                         midpoint = 0, space = "Lab", 
                         na.value = "grey50", guide = "colourbar") +
  geom_segment(data = df_field[sample(nrow(df_field), 10000), ], 
               aes(x = x -Enormx/2, y = y - Enormy/2, xend = x + Enormx/2, yend = y + Enormy/2), 
               size = 0.2, alpha = 0.6) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.08, 0.74),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  coord_fixed() + #preserve aspect ratio
  ggtitle("E-field from two charges: q1 = -0.3q2")
ggsave('E_-03_1.png', height=4.5, width=6)


# E-field and V-field plot
ggplot() + 
  stat_contour(data = df_field, aes(x = x, y = y, z = V, color=..level..), breaks=seq(-1, 1, 0.05)^5, size = 0.6) +  
  scale_colour_gradient2(name = "V", low = "blue", mid = "green", high = "red",
                         midpoint = 0, space = "Lab", 
                         na.value = "grey50", guide = "colourbar") +
  geom_segment(data = df_field[sample(nrow(df_field), 10000), ], 
               aes(x = x -Enormx/2, y = y - Enormy/2, xend = x + Enormx/2, yend = y + Enormy/2), 
               size = 0.2, alpha = 0.6) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.08, 0.74),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  coord_fixed() + 
  ggtitle("V-field and E-field of two charges: q1 = -0.3q2")
ggsave('V_-03_1.png', height=4.5, width=6)


# B-field plot
ggplot() + 
  stat_contour(data = df_field, aes(x = x, y = y, z = B, color=..level..), breaks=seq(0,1.2,0.05)^3, size = 0.6) +  
  scale_colour_gradient2(name = "|B|", low = "blue", mid = "green", high = "red",
                         midpoint = 0, space = "Lab", 
                         na.value = "grey50", guide = "colourbar") +
  geom_segment(data = df_field[sample(nrow(df_field), 8000), ], 
               aes(x = x - Bnormx/2, y = y - Bnormy/2, xend = x + Bnormx/2, yend = y + Bnormy/2), 
               size = 0.2, alpha = 0.6) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.08, 0.74),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  coord_fixed() + 
  ggtitle("B-field of ∞ parallel wires: I1 = -0.3I2")
ggsave('B_-03_1.png', height=4.5, width=6)




#---- EXTRA PLOTS

df_field <- calc_fields(c('x'=0.33, 'y'=0.5, 's'=0.3), c('x'=0.66, 'y'=0.5, 's'=1))

# E field plot
ggplot(df_field, aes(x = x, y = y)) + 
  stat_contour(aes(z = E, color=..level..),breaks=-7:12, size = 0.6) +  
  scale_colour_gradient2(name = "log |E|", low = "blue", mid = "green", high = "red",
                         midpoint = 0, space = "Lab", 
                         na.value = "grey50", guide = "colourbar") +
  geom_segment(data = df_field[sample(nrow(df_field), 10000), ], 
               aes(x = x -Enormx/2, y = y - Enormy/2, xend = x + Enormx/2, yend = y + Enormy/2), 
               size = 0.2, alpha = 0.6) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.08, 0.74),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  coord_fixed() + #preserve aspect ratio
  ggtitle("E-field from two charges: q1 = 0.3q2")
ggsave('E_03_1.png', height=4.5, width=6)


# E-field and V-field plot
ggplot() + 
  stat_contour(data = df_field, aes(x = x, y = y, z = V, color=..level..), breaks=seq(0,1.5,0.05)^3, size = 0.6) +  
  scale_colour_gradient2(name = "V", low = "blue", mid = "green", high = "red",
                         midpoint = 0, space = "Lab", 
                         na.value = "grey50", guide = "colourbar") +
  geom_segment(data = df_field[sample(nrow(df_field), 10000), ], 
               aes(x = x -Enormx/2, y = y - Enormy/2, xend = x + Enormx/2, yend = y + Enormy/2), 
               size = 0.2, alpha = 0.6) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.08, 0.74),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  coord_fixed() + 
  ggtitle("V-field and E-field of two charges: q1 = 0.3q2")

ggsave('V_03_1.png', height=4.5, width=6)


# B-field plot
ggplot() + 
  stat_contour(data = df_field, aes(x = x, y = y, z = B, color=..level..), breaks=seq(0,1.2,0.05)^3, size = 0.6) +  
  scale_colour_gradient2(name = "|B|", low = "blue", mid = "green", high = "red",
                         midpoint = 0, space = "Lab", 
                         na.value = "grey50", guide = "colourbar") +
  geom_segment(data = df_field[sample(nrow(df_field), 8000), ], 
               aes(x = x - Bnormx/2, y = y - Bnormy/2, xend = x + Bnormx/2, yend = y + Bnormy/2), 
               size = 0.2, alpha = 0.6) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.08, 0.74),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  coord_fixed() + 
  ggtitle("B-field of ∞ parallel wires: I1 = 0.3I2")
ggsave('B_03_1.png', height=4.5, width=6)

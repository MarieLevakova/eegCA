# Functions to plot estimated parameters of VECM model for EEG data
# An EEG setup with 64 channels is expected

# Set up coordinates of EEG channels
channels <- data.frame(no = 1:64,
                       name = c(paste("A", 1:32, sep = ""), paste("B", 1:32, sep = "")),
                       x = numeric(64),
                       y = numeric(64))

channels[48,3:4] <- c(0,0)
channels[33,3:4] <- c(0,0.8)
channels[37,3:4] <- c(0,0.6)
channels[38,3:4] <- c(0,0.4)
channels[47,3:4] <- c(0,0.2)
channels[32,3:4] <- c(0,-0.2)
channels[31,3:4] <- c(0,-0.4)
channels[30,3:4] <- c(0,-0.6)
channels[29,3:4] <- c(0,-0.8)
channels[28,3:4] <- c(0,-1)
channels[15,3:4] <- c(-0.8,0)
channels[14,3:4] <- c(-0.6,0)
channels[13,3:4] <- c(-0.4,0)
channels[12,3:4] <- c(-0.2,0)
channels[49,3:4] <- c(0.2,0)
channels[50,3:4] <- c(0.4,0)
channels[51,3:4] <- c(0.6,0)
channels[52,3:4] <- c(0.8,0)

channels[43,3] <- 0.8*cos(pi/10)
channels[43,4] <- 0.8*sin(pi/10)

channels[42,3] <- 0.8*cos(2*pi/10)
channels[42,4] <- 0.8*sin(2*pi/10)

channels[35,3] <- 0.8*cos(3*pi/10)
channels[35,4] <- 0.8*sin(3*pi/10)

channels[34,3] <- 0.8*cos(4*pi/10)
channels[34,4] <- 0.8*sin(4*pi/10)

channels[1,3] <- 0.8*cos(6*pi/10)
channels[1,4] <- 0.8*sin(6*pi/10)

channels[2,3] <- 0.8*cos(7*pi/10)
channels[2,4] <- 0.8*sin(7*pi/10)

channels[7,3] <- 0.8*cos(8*pi/10)
channels[7,4] <- 0.8*sin(8*pi/10)

channels[8,3] <- 0.8*cos(9*pi/10)
channels[8,4] <- 0.8*sin(9*pi/10)

channels[16,3] <- 0.8*cos(11*pi/10)
channels[16,4] <- 0.8*sin(11*pi/10)

channels[23,3] <- 0.8*cos(12*pi/10)
channels[23,4] <- 0.8*sin(12*pi/10)

channels[24,3] <- cos(12*pi/10)
channels[24,4] <- sin(12*pi/10)

channels[25,3] <- 0.8*cos(13*pi/10)
channels[25,4] <- 0.8*sin(13*pi/10)

channels[27,3] <- 0.8*cos(14*pi/10)
channels[27,4] <- 0.8*sin(14*pi/10)

channels[64,3] <- 0.8*cos(16*pi/10)
channels[64,4] <- 0.8*sin(16*pi/10)

channels[62,3] <- 0.8*cos(17*pi/10)
channels[62,4] <- 0.8*sin(17*pi/10)

channels[60,3] <- 0.8*cos(18*pi/10)
channels[60,4] <- 0.8*sin(18*pi/10)

channels[61,3] <- cos(18*pi/10)
channels[61,4] <- sin(18*pi/10)

channels[53,3] <- 0.8*cos(19*pi/10)
channels[53,4] <- 0.8*sin(19*pi/10)

channels[45,3] <- 0.4*cos(pi/10)
channels[45,4] <- 0.8*sin(pi/10)

channels[40,3] <- 0.4*cos(2*pi/10)
channels[40,4] <- 0.8*sin(2*pi/10)

channels[36,3] <- 0.4*cos(3*pi/10)
channels[36,4] <- 0.8*sin(3*pi/10)

channels[3,3] <- 0.4*cos(7*pi/10)
channels[3,4] <- 0.8*sin(7*pi/10)

channels[5,3] <- 0.4*cos(8*pi/10)
channels[5,4] <- 0.8*sin(8*pi/10)

channels[10,3] <- 0.4*cos(9*pi/10)
channels[10,4] <- 0.8*sin(9*pi/10)

channels[18,3] <- 0.4*cos(11*pi/10)
channels[18,4] <- 0.8*sin(11*pi/10)

channels[21,3] <- 0.4*cos(12*pi/10)
channels[21,4] <- 0.8*sin(12*pi/10)

channels[26,3] <- 0.4*cos(13*pi/10)
channels[26,4] <- 0.8*sin(13*pi/10)

channels[63,3] <- 0.4*cos(17*pi/10)
channels[63,4] <- 0.8*sin(17*pi/10)

channels[58,3] <- 0.4*cos(18*pi/10)
channels[58,4] <- 0.8*sin(18*pi/10)

channels[55,3] <- 0.4*cos(19*pi/10)
channels[55,4] <- 0.8*sin(19*pi/10)

channels[44,3] <- 0.6*cos(pi/10)
channels[44,4] <- 0.8*sin(pi/10)

channels[41,3] <- 0.6*cos(2*pi/10)
channels[41,4] <- 0.8*sin(2*pi/10)

channels[46,3] <- 0.2*cos(pi/10)
channels[46,4] <- 0.8*sin(pi/10)

channels[39,3] <- 0.2*cos(2*pi/10)
channels[39,4] <- 0.8*sin(2*pi/10)

channels[6,3] <- 0.6*cos(8*pi/10)
channels[6,4] <- 0.8*sin(8*pi/10)

channels[9,3] <- 0.6*cos(9*pi/10)
channels[9,4] <- 0.8*sin(9*pi/10)

channels[4,3] <- 0.2*cos(8*pi/10)
channels[4,4] <- 0.8*sin(8*pi/10)

channels[11,3] <- 0.2*cos(9*pi/10)
channels[11,4] <- 0.8*sin(9*pi/10)

channels[17,3] <- 0.6*cos(11*pi/10)
channels[17,4] <- 0.8*sin(11*pi/10)

channels[22,3] <- 0.6*cos(12*pi/10)
channels[22,4] <- 0.8*sin(12*pi/10)

channels[19,3] <- 0.2*cos(11*pi/10)
channels[19,4] <- 0.8*sin(11*pi/10)

channels[20,3] <- 0.2*cos(12*pi/10)
channels[20,4] <- 0.8*sin(12*pi/10)

channels[59,3] <- 0.6*cos(18*pi/10)
channels[59,4] <- 0.8*sin(18*pi/10)

channels[54,3] <- 0.6*cos(19*pi/10)
channels[54,4] <- 0.8*sin(19*pi/10)

channels[57,3] <- 0.2*cos(18*pi/10)
channels[57,4] <- 0.8*sin(18*pi/10)

channels[56,3] <- 0.2*cos(19*pi/10)
channels[56,4] <- 0.8*sin(19*pi/10)

# A function to plot a matrix as a heatmap in a nice way
plot_matrix <- function(mat, main = "", limits = c(min(mat), max(mat)),
                        x.lab = "Input channel", y.lab = "Output channel"){

  # Rearrange the matrix into a data frame with x and y coordinates and the value to be plotted
  melted_mat <- data.frame(from = rep(1:ncol(mat), each = nrow(mat)),
                           to = 1:nrow(mat),
                           value = matrix(mat, ncol = 1, byrow = F))

  the_plot <- ggplot(data = melted_mat, aes(x = from, y = 65-to, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, space = "Lab",
                         name = "Value", limits = limits) +
    theme_light()  +
    scale_y_continuous(breaks = 65 - c(1, 9, 17, 25, 33, 41, 49, 57),
                       labels = c("A1", "A9", "A17", "A25", "B1", "B9", "B17", "B25"),
                       expand = expansion(mult = c(0, 0))) +
    scale_x_continuous(expand = expansion(mult = c(0, 0))) +
    theme(axis.text.x = element_text(angle = 90)) +
    xlab(x.lab) + ylab(y.lab) + ggtitle(main)
  if(nrow(mat) == ncol(mat)){ # If mat is a square matrix (Pi expected), add channel labels on the x-axis too
    the_plot <- the_plot +
      scale_x_continuous(breaks = c(1, 9, 17, 25, 33, 41, 49, 57),
                         labels = c("A1", "A9", "A17", "A25", "B1", "B9", "B17", "B25"),
                         expand = expansion(mult = c(0, 0)))
  }
  the_plot
}

# A function to plot Pi as a network
plot_network <- function(mat, channel.sel = 1:64, main = "", text = "", excluded = NULL,
                         reference = NULL, limits = c(min(mat), max(mat))){

  amplitude <- max(c(abs(mat), limits)) # maximum abs. value of a connection

  n <- 100 #???

  # Rearrange the matrix into a data frame with columns "from", "to" and "strength"
  mat.frame <- data.frame(from = rep(1:ncol(mat), each = nrow(mat)),
                          to = rep(1:nrow(mat), ncol(mat)),
                          strength = matrix(mat, ncol = 1, byrow = F))
  # Add exact coordinates of the channels
  mat.frame$from.x <- channels$x[mat.frame$from]
  mat.frame$from.y <- channels$y[mat.frame$from]
  mat.frame$to.x <- channels$x[mat.frame$to]
  mat.frame$to.y <- channels$y[mat.frame$to]
  # Re-order the rows from the weakest to the strongest connection (the strongest will be plotted last)
  mat.frame <- mat.frame[order(abs(mat.frame$strength)),]

  the_plot <- ggplot(channels) +
    # Plot the contour of the head, nose and ears
    geom_ellipse(aes(x0 = 0, y0 = 0.8, a = 0.1, b = 0.35, angle = 0)) +
    geom_ellipse(aes(x0 = -1, y0 = -0.1, a = 0.05, b = 0.25, angle = 0)) +
    geom_ellipse(aes(x0 = 1, y0 = -0.1, a = 0.05, b = 0.25, angle = 0)) +
    geom_circle(aes(x0 = 0, y0 = -0.1, r = 1), fill = "white") +
    # Plot the arrows
    geom_segment(data = mat.frame,
                 aes(x = from.x, y = from.y, xend = to.x, yend = to.y,
                     colour = strength,
                     size = abs(strength/(amplitude))),
                 arrow = arrow(length = unit(0.03, "npc"))) +
    # Add dots at electrode locations
    geom_point(aes(x = x, y = y), alpha = 0.25) +
    # Mark excluded channels by crosses
    geom_point(data = channels[excluded,], aes(x = x, y = y), shape = 4, size = 3) +
    theme_void() +
    scale_colour_gradient2(low = "blue", high = "red", mid = "white",
                           midpoint = 0, space = "Lab",
                           name = "Value", limits = range(c(mat.frame$strength, limits))) +
    scale_size_continuous(range = c(0, 2)) +
    guides(size = FALSE)

  the_plot
}

# A function to plot drifts, loadings, weights of cointegration vectors as levels
# at the channels on the scalp
plot_levels <- function(coefs, main = "", limits = range(coefs), excluded = NULL){
  amplitude <- max(c(abs(coefs), limits))
  the_plot <- ggplot(channels) +
    # Plot the head contour
    geom_ellipse(aes(x0 = 0, y0 = 0.8, a = 0.1, b = 0.35, angle = 0)) +
    geom_ellipse(aes(x0 = -1, y0 = -0.1, a = 0.05, b = 0.25, angle = 0)) +
    geom_ellipse(aes(x0 = 1, y0 = -0.1, a = 0.05, b = 0.25, angle = 0)) +
    geom_circle(aes(x0 = 0, y0 = -0.1, r = 1), fill = "white") +
    # Plot channel locations filled with a colour corresponding to "coefs"
    geom_circle(aes(x0 = x, y0 = y, r=0.06, fill = coefs)) +
    # Cross-out the exluded channels
    geom_point(data = channels[excluded,], aes(x = x, y = y), shape = 4, size = 2) +
    theme_void() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0, space = "Lab",
                         name = "Value", limits = range(c(coefs, limits)))
  the_plot
}

# A function to combine 6 plots for the 6 datasets into 1 plot
plot_6subsets <- function(p1, p2, p3, p4, p5, p6){

  # Extract legend from one of the plots
  legend1 <- get_legend(p1 + theme(legend.position = "bottom"))

  # Assemble first row of the plot
  p_row1 <- plot_grid(p1 + theme(legend.position = "none"),
                      p3 + theme(legend.position = "none"),
                      p5 + theme(legend.position = "none"),
                      labels = c("Fixation", "Stimulation", "Masking"),
                      label_x = 0.5, hjust = 0.5, label_size = 10, nrow = 1, greedy = TRUE) +
    theme(plot.margin = unit(c(1, 0, -1, 2.5), "cm"))

  # Assemble second row of the plot
  p_row2 <- plot_grid(p2 + theme(legend.position = "none"),
                      p4 + theme(legend.position = "none"),
                      p6 + theme(legend.position = "none"), nrow = 1, greedy = TRUE) +
    theme(plot.margin = unit(c(-1.5, 0, 0, 2.5), "cm"))

  # Glue the two rows together
  p_nolegend <- plot_grid(p_row1, NULL, p_row2, ncol = 1, align = 'vh', rel_heights = c(1, -0.15, 1),
                          labels = c("Correct", "", "Incorrect"),
                          label_x = 0, label_y = 0.5, hjust = 0, label_size = 10, greedy = TRUE)

  # Add the legend
  p <- plot_grid(p_nolegend, legend1, ncol = 1, align = 'vh', rel_heights = c(1, 0.1))
  p
}

# Plotting the Renyi Divergences
library(latex2exp)
library(lattice)

# Get data from ALICE and adjust
# so that we can plot.
test <- read.csv('~/Desktop/results.txt')[, -1]
test <- test[, 1:101]
test <- as.matrix(test)

# x and y coordinates
theta_x <-  ((1:151) - 75) / 50
theta_1 <-  ((1:101) - 50) / 50

# re-label the columns and rows
# such that the plot is visually acceptable
colnames(test) <- theta_1
colnames(test) <- rep('', length(theta_1))
colnames(test)[c(1,25,50,
                 75,100)] <- theta_1[c(1,25,50,
                                       75,100)]
colnames(test)[1] <- -1

rownames(test) <- rep('', length(theta_x))
rownames(test)[c(1,25,50,75,
                 100,125,150)] <- theta_x[c(1,25,50,
                                            75,100,
                                            125,150)]
rownames(test)[1] <- -1.5


# Defining custom color palette
colors <- c("green", 
            colorRampPalette(c("orange", 
                               "red"))(100))
breaks <- c(seq(min(test, na.rm= TRUE), 
                1, length.out = 2), 
            seq(1.01, max(test, na.rm = TRUE), 
                length.out = 100))


# Plot results.
levelplot(test, 
          xlab = TeX(r'($\theta_{x}$)'), 
          ylab = TeX(r'($\theta_{1}$)'),
          col.regions = colors, at = breaks)

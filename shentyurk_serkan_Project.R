rm(list = ls())

library(cellWise)
library(ggplot2)
library(rrcov)

makecolors <- function(res) {
  colors <- rep("grey", length(res$sd))
  id.ol <- which(res$od > res$cutoff.od)
  id.gl <- which(res$sd > res$cutoff.sd)
  id.bl <- which((res$sd > res$cutoff.sd)&(res$od > res$cutoff.od))
  colors[id.ol] <- "orange"
  colors[id.gl] <- "blue"
  colors[id.bl] <- "red"
  out <- list(colors = colors, id.ol = id.ol, id.gl = id.gl, id.bl = id.bl)
}


Biscuit.X <- read.table("biscuitX.txt", header = F, sep = "") 
# 40 NIR spectra of biscuit dough, with measurement every 2 nm from 1200 to 2400 nm.
# n = 40, d = 600 

Biscuit.Y <- read.table("biscuitY.txt", header = F, sep = "")
# the percentages of  constituents in the dough: y1 = fat, y2 = flour, y3 = sucrose, y4 = water
set.seed(0833419) # change here with your student number
rsub <- sample(4)[1] # 3
y.reduced <- Biscuit.Y[, rsub] # your selected response variable - y3 = sucrose


### Work Plan:
# 1. Apply MacroPCA and CPCA.
# 2. Calculate and plot orthogonal distance against score distances.
# 3. Plot the loadings of the variables and compare the values.
# 4. Identify the underlying reasons of outlier.
# 5. Apply LTS estimator and classical linear regression to 
#    selected principle components of MacroPCA and CPCA.
# 6. Calculate and plot standardised LS (LTS) residuals against 
#   Mahalanobis distance (Robust distance - MCD) for classical regression (LTS estimation) 
# 7. Compare outputs

### Part 0 - Initial data investigation

(x_dim = dim(Biscuit.X))
(na_counts_X <- apply(Biscuit.X, 1, function(row) sum(is.na(row))))
(na_counts_Y = sum(is.na(y.reduced))
)


### Part 1 - Apply MacroPCA and CPCA.
## MacroPCA
DDCpars <- list(fastDDC = TRUE, silent = TRUE)
MacroPCApars <- list(DDCpars = DDCpars, scale = TRUE, silent = TRUE)

set.seed(0833419)
macroPCA_result = MacroPCA(Biscuit.X, MacroPCApars = MacroPCApars)
(macroPCA_cumulativeVar <- macroPCA_result$cumulativeVar)
# Create a data frame for plotting
macro_PCA_df <- data.frame(
  Component = 1:10,
  CumulativeVariance = macroPCA_cumulativeVar
)
ggplot(macro_PCA_df, aes(x = Component, y = CumulativeVariance)) +
  geom_line(color = "blue", size = 1) +  # Line plot
  geom_point(color = "red", size = 3) +  # Points on the line
  labs(title = "Cumulative Variance Explained by Components - MacroPCA",
       x = "Number of Components",
       y = "Cumulative Variance") +
  theme_minimal() +  # A nice minimal theme
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center and style title
    axis.title.x = element_text(size = 12, face = "bold"),  # Style x-axis title
    axis.title.y = element_text(size = 12, face = "bold")   # Style y-axis title
  )
# k = 10 explains 90.17% variability
macroPCA_result = MacroPCA(Biscuit.X, k = 10)



## PCA
set.seed(0833419)
cPCA_result <- PcaClassic(Biscuit.X, k = 10, scale = TRUE)
summary(cPCA_result)
cPCA_df <- data.frame(
  Component = 1:10,
  CumulativeVariance = cumsum(cPCA_result$eigenvalues)/cPCA_result$totvar0
)
ggplot(cPCA_df, aes(x = Component, y = CumulativeVariance)) +
  geom_line(color = "blue", size = 1) +  # Line plot
  geom_point(color = "red", size = 3) +  # Points on the line
  labs(title = "Cumulative Variance Explained by Components - cPCA",
       x = "Number of Components",
       y = "Cumulative Variance") +
  theme_minimal() +  # A nice minimal theme
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center and style title
    axis.title.x = element_text(size = 12, face = "bold"),  # Style x-axis title
    axis.title.y = element_text(size = 12, face = "bold")   # Style y-axis title
  )
# k=10 89%


### Part 2 - Orthogonal distance against Score distances
## Distance plot for MacroPCA
set.seed(0833419)
outlierMap(macroPCA_result,title="MacroPCA outlier map",
           labelOut=TRUE, id = 3)

#same plot but alternative
bad_leverage_macro = macroPCA_result$OD > macroPCA_result$cutoffOD & macroPCA_result$SD > macroPCA_result$cutoffSD
orthogonal_macro = macroPCA_result$OD > macroPCA_result$cutoffOD & macroPCA_result$SD < macroPCA_result$cutoffSD
good_leverage_macro  = macroPCA_result$OD < macroPCA_result$cutoffOD & macroPCA_result$SD > macroPCA_result$cutoffSD
colors = rep(NA, 20)
colors[bad_leverage_macro] = 'red'
colors[orthogonal_macro] = 'orange'
colors[good_leverage_macro] = 'blue'
colors[is.na(colors)] = 'gray'
# colors[17] = 'pink' # outlier in regular case
# colors[23] = 'pink' # outlier in regular case
plot(macroPCA_result$SD, macroPCA_result$OD, 
     col = colors, 
     main = 'MacroPCA', 
     xlab = 'Score Distance', 
     ylab = 'Orthogonal Distance', 
     pch = 16)
abline(v = macroPCA_result$cutoffSD, col = "black", lwd = 2, lty = 2)  
abline(h = macroPCA_result$cutoffOD, col = "black", lwd = 2, lty = 2)  

## Distance plot for PCA

res = cPCA_result
outliers <- makecolors(res)
colors <- outliers$colors
plot(res, pch = 16, col = colors, main = "CPCA")

# 17 orthogonal
# 23 good leverage


### Part 3 - Plot the loadings of the variables
## Plot loadings of MacroPCA
unique_names <- unique(rownames(macroPCA_result$loadings))
name_colors <- setNames(colorRampPalette(c("blue", "red"))(length(unique_names)), unique_names)

plot(1, type="n", xlim=c(1, 10), ylim=range(macroPCA_result$loadings), 
     xlab="Component", ylab="Loading", main="MacroPCA Loadings")

for (i in 1:10) {
  col_values <- macroPCA_result$loadings[,i]
  text(rep(i, length(col_values)), col_values, 
       labels = rownames(macroPCA_result$loadings), 
       col = name_colors[rownames(macroPCA_result$loadings)], cex = 0.7)
}
matplot(macroPCA_result$loadings, type = "l", lty = 1, xlab = 'Samples', ylab = 'Loading', main = "MacroPCA Loadings")

## Plot loadings of PCA
unique_names <- unique(rownames(cPCA_result$loadings))
plot(1, type="n", xlim=c(1, 10), ylim=range(cPCA_result$loadings), 
     xlab="Component", ylab="Loading", main="cPCA Loadings")
for (i in 1:10) {
  col_values <- cPCA_result$loadings[,i]
  text(rep(i, length(col_values)), col_values, 
       labels = rownames(cPCA_result$loadings), 
       col = name_colors[rownames(cPCA_result$loadings)], cex = 0.7)
}
matplot(cPCA_result$loadings, type = "l", lty = 1, xlab = 'Variable', ylab = 'Loading', main = "MacroPCA Loadings")



### Part 4 - Investigate the outliers try to make sense
## Outlier investigation for MacroPCA
good_leverage_macro_indices = c(1:40)[good_leverage_macro]
bad_leverage_macro_indices = c(1:40)[bad_leverage_macro]
orthogonal_macro_indices = c(1:40)[orthogonal_macro]
(all_outlier_indices = sort(c(good_leverage_macro_indices, 
                        bad_leverage_macro_indices, 
                        orthogonal_macro_indices))
)

cellMap(macroPCA_result$stdResid, ncolumnsinblock = 5, nrowsinblock = 1, mTitle = 'MacroPCA - Cell Map')



## Outlier investigation for PCA

# 17 orthogonal
# 23 good leverage

scores <- cPCA_result@scores
loadings <- cPCA_result@loadings
reconstructed_data <- scores %*% t(loadings)
residuals <- Biscuit.X - reconstructed_data
standardized_residuals <- scale(residuals)
cellMap(standardized_residuals, ncolumnsinblock = 5, nrowsinblock = 1, mTitle = 'cPCA - Cell Map')

  
  
### Part 5 - LTS or classical Linear Regression Analysis
macroPCA_scores = macroPCA_result$scores
macroPCA_to_regression = data.frame(macroPCA_scores, Y = y.reduced)
cPCA_scores = cPCA_result$scores
cPCA_to_regression = data.frame(cPCA_scores, Y = y.reduced)


## LTS analysis of MacroPCA outputs
set.seed(0833419)
LTSModel_macro <- ltsReg(Y~., data = macroPCA_to_regression)
summary(LTSModel_macro)

plot(LTSModel_macro) 
plot(LTSModel_macro,which="rdiag") # 7 20 23 24

LTSModel_macro$lts.wt
which(unname(LTSModel_macro$lts.wt)==0)
which(unname(LTSModel_macro[["raw.weights"]])==0)

Ind <- which((LTSModel_macro$RD > 4.5)&(abs(LTSModel_macro$resid)<2.2))
#Ind <- 421:471
boxplot(macroPCA_to_regression[,1:10])
for(i in 1:10){ points(x = rep(i, length(Ind)), macroPCA_to_regression[Ind,i], col = "red") }



## LTS analysis of cPCA outputs
set.seed(0833419)
LTSModel_classic <- ltsReg(Y~., data = cPCA_to_regression)
summary(LTSModel_classic)

plot(LTSModel_classic)
plot(LTSModel_classic,which="rdiag") # 21 22 23

LTSModel_classic$lts.wt
which(unname(LTSModel_classic$lts.wt)==0)
which(unname(LTSModel_classic[["raw.weights"]])==0)

Ind <- which((LTSModel_classic$RD > 4.5)&(abs(LTSModel_classic$resid)<2.2))
#Ind <- 421:471
boxplot(cPCA_to_regression[,1:10])
for(i in 1:10){ points(x = rep(i, length(Ind)), cPCA_to_regression[Ind,i], col = "red") }

  


## Linear Regression analysis of MacroPCA outputs
set.seed(0833419)
LSModel_macro <- lm(Y~., data = macroPCA_to_regression)
summary(LSModel_macro)

plot(LSModel_macro, which = 2) # 23 24 - 7 20
plot(LSModel_macro,which=4)
plot(LSModel_macro, which=5)
plot(LSModel_macro, which = 6)

standardlm <- rstandard(LSModel_macro)
MDLS <- sqrt(mahalanobis(macroPCA_to_regression[,1:10], colMeans(macroPCA_to_regression[,1:10]), cov(macroPCA_to_regression[,1:10])))

MDLS_thresh = sqrt(qchisq(0.975,ncol(macroPCA_to_regression)-2))
standardised_ls_threshold_1 = -2.5
standardised_ls_threshold_2 = 2.5

vertical_outliers = (MDLS < MDLS_thresh & standardlm > standardised_ls_threshold_2) | (MDLS < MDLS_thresh & standardlm < standardised_ls_threshold_2)
good_leverage_points = MDLS > MDLS_thresh & standardlm > standardised_ls_threshold_1 & standardlm < standardised_ls_threshold_2
regular_obs = MDLS <  MDLS_thresh & standardlm > standardised_ls_threshold_1 & standardlm < standardised_ls_threshold_2
colors = rep(NA, dim(macroPCA_to_regression)[1])
colors[vertical_outliers] = 'orange'
colors[good_leverage_points] = 'green'
colors[regular_obs] = 'gray'
colors[is.na(colors)] = 'red'

plot(MDLS, standardlm, col = colors, pch = 16, ylab = "Standardised LS Residual", xlab = "Mahalanobis Distance", main = "LS to MacroPCA Outputs")
abline(h = standardised_ls_threshold_1)
abline(h = standardised_ls_threshold_2)
abline(v= MDLS_thresh)


## Linear Regression analysis of cPCA outputs
set.seed(0833419)
LSModel_classic <- lm(Y~., data = cPCA_to_regression)
summary(LSModel_classic)

plot(LSModel_macro, which = 2) # 7 23 24 - 20
plot(LSModel_classic,which=4)
plot(LSModel_classic, which=5)
plot(LSModel_classic, which = 6)

standardlm <- rstandard(LSModel_classic)
MDLS <- sqrt(mahalanobis(cPCA_to_regression[,1:10], colMeans(cPCA_to_regression[,1:10]), cov(cPCA_to_regression[,1:10])))

MDLS_thresh = sqrt(qchisq(0.975,ncol(cPCA_to_regression)-2))
standardised_ls_threshold_1 = -2.5
standardised_ls_threshold_2 = 2.5

vertical_outliers = (MDLS < MDLS_thresh & standardlm > standardised_ls_threshold_2) | (MDLS < MDLS_thresh & standardlm < standardised_ls_threshold_2)
good_leverage_points = MDLS > MDLS_thresh & standardlm > standardised_ls_threshold_1 & standardlm < standardised_ls_threshold_2
regular_obs = MDLS <  MDLS_thresh & standardlm > standardised_ls_threshold_1 & standardlm < standardised_ls_threshold_2
colors = rep(NA, dim(cPCA_to_regression)[1])
colors[vertical_outliers] = 'orange'
colors[good_leverage_points] = 'green'
colors[regular_obs] = 'gray'
colors[is.na(colors)] = 'red'

plot(MDLS, standardlm, col = colors, pch = 16, ylab = "Standardised LS Residual", xlab = "Mahalanobis Distance", main = "LS to cPCA Outputs")
abline(h = standardised_ls_threshold_1)
abline(h = standardised_ls_threshold_2)
abline(v= MDLS_thresh)



## Re-do classical linear regression but remove the outliers
cPCA_to_regression_reduced = cPCA_to_regression[-c(7,20,23,24),]
macroPCA_to_regression_reduced = macroPCA_to_regression[-c(7,20,23,24),]

set.seed(0833419)

LSModel_classic_reduced <- lm(Y~., data = cPCA_to_regression_reduced)
(reduced_classic_summary = summary(LSModel_classic_reduced))
classic_classic_summary = summary(LSModel_classic)
(lts_classic_summary = summary(LTSModel_classic))


plot(1:10, reduced_classic_summary$coefficients[2:11,1], pch = 16, col = 'red', 
     ylim = range(c(reduced_classic_summary$coefficients[2:11,1], 
                    lts_classic_summary$coefficients[2:11,1])),
     xlab = "Components",
     ylab = "Coefficients",
     main = "The Regression Coefficients of cPCA Outputs")
points(1:10+0.1, lts_classic_summary$coefficients[2:11,1], pch = 16, col = 'blue')
points(1:10+0.2, classic_classic_summary$coefficients[2:11,1], pch = 16, col = 'black')
legend("topleft", legend = c("Outliers Omitted\nRegression", "LTS Regression", "Linear Regression"), col = c("red", "blue", "black"), pch = 16)



set.seed(0833419)
LSModel_macro_reduced <- lm(Y~., data = macroPCA_to_regression_reduced)

(reduced_macro_summary = summary(LSModel_macro_reduced))
classic_macro_summary = summary(LSModel_macro)
(lts_macro_summary = summary(LTSModel_macro))

plot(1:10, reduced_macro_summary$coefficients[2:11,1], pch = 16, col = 'red', 
     ylim = range(c(reduced_macro_summary$coefficients[2:11,1], 
                    lts_macro_summary$coefficients[2:11,1])),
     xlab = "Components",
     ylab = "Coefficients",
     main = "The Regression Coefficients of MacroPCA Outputs")
points(1:10+0.1, lts_macro_summary$coefficients[2:11,1], pch = 16, col = 'blue')
points(1:10+0.2, classic_macro_summary$coefficients[2:11,1], pch = 16, col = 'black')
legend("topleft", legend = c("Outliers Omitted\nRegression", "LTS Regression", "Linear Regression"), col = c("red", "blue", "black"), pch = 16)


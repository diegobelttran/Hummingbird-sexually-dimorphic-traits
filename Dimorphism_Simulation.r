# THE EVOLUTION OF DIMORPHIC TRAITS IN ALTITUDINAL AND ENVIRONMENTAL GRADIENTS: AN INTERPLAY BETWEEN NATURAL AND SEXUAL SELECTION IN HUMMINGBIRDS
# DIEGO F. BELTRÁN, MARCELO ARAYA-SALAS, JUAN L. PARRA, F. GARY STILES AND ALEJANDRO RICO-GUEVARA
# MANUSCRIPT SUBMITTED TO PROCEEDINGS OF THE ROYAL SOCIETY B

##################
Simulating multidimensional dimorphism

 
Purpose
Evalute the performance of three commonly used metrics of quantifying sexual dimorphism in multidimensional traits:
Euclidean distance between sexes from all PCs with no correction whatsoever
Euclidean distance between sexes from all PCs weighted by the amount of variance explained by each PC
First PC from a PCA on male-to-female ratio for all morphological variables (all ratios where “standardized” so high values always indicates high dimorphism regardless of whether males or females are bigger)
 

Simulation parameters
The function sim_pca_dimorph (code chunk below) simulated a variety of scenarios for sexual dimorphism with the following parameters:

117 spp (argument spp)
13 variables (argument nvars) drawn from a normal-gamma distribution (gamma was used to avoid negative values)
3 different numbers of dimorphic variables: 10, 3 and 0 (no dimorphism)
Dimorphic variables can be the same for all species or be choosen randomly for each species and are likely to be different among dimorphic species
~ 1/3 of species monomorphic, ~1/3 with low dimorphism and ~1/3 with high dimorphism
Three effect sizes: low = 2, moderate = 3, high = 4 (he effect size controls the number of times the differences between males and females increases in more dimorphic species increases compare to a lower level of dimorphism)
Even in monomorphic species females don’t have exactly the same values than males (some error is added)
All variables are z-transformed within PCA (except for ratios)
The simulation can control if dimorphic variables are within the most varying variables or not, or both are independently determined (argument varying_and_dimorphic)
Variables can covary (argument covar_vars). Covariance decreases steadly from the first variable. For instance, if covar_vars = 3 then variable 2 shows the highest covariance with 1 and 3 the highest covariance with 2.
Three levels of covariance between dimorphic traits where simulated:
No covarying variables (i.e. very low covariance, all PCs explained about the same amount of variance)
Moderate covariance among some variables (which results in PC1 explaining most of the variance) but no association between which variables covary and which are the most varying variables in the data set
Moderate covariance among some variables and those variables are also among the most varying ones
The function returns the three dimorphism metrics and the dimorphism type of each species. Each model was replicated 1000 times. Each replicate returned the p-value of linear regression models for the dimorphism metrics as response (1 at the time) and the dimorphism categories (high, low and monomorphic) as predictor. So 3 p-values for each replicate. Dimorphism categories were encoded as ordered nominal variable.

The statistical power and false positive rate (type I error) were calculated for each dimorphism metric. Power was defined as the proportion of significant p-values (alpha = 0.05) when there is actual dimorphism (10 or 3 dimorphic variables). A good test should have a power close to 1 when the sample size is decent, as in these simulations. False positive rate was measured as the proportion of significant p-values when there is no dimorphism (0 dimorphic variables)
##################

###### BEGIN CODE
# calculate euclidean distance weighted by PC explained variance
weight_dist <- function(X, pca) {
  # all (q-p)^2
  q_p2 <- sapply(1:ncol(X), function(y) ((X[1, y] - X[2, y]) ^ 2)) 
  
  # weight by PCA explained variance
  weight_q_p2 <- q_p2 * (summary(pca)$importance[2, ] / max(summary(pca)$importance[2, ]))  
  
  # square root of sum
  weight_d <- sqrt(sum(weight_q_p2))
  
  return(weight_d)
  }

# extract p values from lm models
lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}

# simulating function
sim_pca_dimorph <- function(spp = 117, nvars = 13, dimorphic_vars = 3, covar_vars = 5, varying_and_dimorphic = "yes", same_dimorphic = TRUE, effect_size = 2){

  dimorphism <- sample(c("monomorphic", "low dimorphism", "high dimorphism"), size = spp, replace = TRUE)
  dimorphism_df <- replicate(n = nvars, rep("monomorphic", spp))
  
  if (dimorphic_vars > 0)
    for( i in 1:spp){
    if (same_dimorphic){
       if (dimorphism[i] == "low dimorphism")
       dimorphism_df[i, 1:dimorphic_vars] <- "low dimorphism"
     
      if (dimorphism[i] == "high dimorphism")
       dimorphism_df[i, 1:dimorphic_vars] <- "high dimorphism"
      } else     {
      if (dimorphism[i] == "low dimorphism")
      dimorphism_df[i, sample(1:13, dimorphic_vars)] <- "low dimorphism"
      if (dimorphism[i] == "high dimorphism")
      dimorphism_df[i, sample(1:13, dimorphic_vars)] <- "high dimorphism"
      }
    } 
  
  sds <- sample(seq(1, 3, length.out = 3), size = nvars, replace = TRUE)
  if (varying_and_dimorphic != "random")
  sds <- sort(sds, decreasing = if(varying_and_dimorphic == "yes") TRUE else FALSE)
  means <- sample(seq(10, 100, length.out = nvars))
  
  males_vars <- list()
    
    for(x in 1:nvars)
    { 
      if(x <= covar_vars & x != 1)  
        males <- males_vars[[x - 1]] + rnorm(n = spp, mean = 0, sd = sds[x - 1]) else
          males <- rnorm(n = spp, mean= means[x], sd=sds[x]) + rgamma(n = spp, shape= 2, scale= 2)
    
        males_vars[[x]] <- males
        
        }
        
  vars <- list()
    
    for(w in 1:nvars)
    { 
      females <- males_vars[[w]]
      
      # monomorphic
    females[dimorphism_df[ , w] == "monomorphic"] <- females[dimorphism_df[ , w] == "monomorphic"] + rnorm(n = sum(dimorphism_df[ , w] == "monomorphic"), mean = 0, sd = sd(males) * 1/8)
  
    # low dimorphism
    females[dimorphism_df[ , w] == "low dimorphism"] <- females[dimorphism_df[ , w] == "low dimorphism"] + rnorm(n = sum(dimorphism_df[ , w] == "low dimorphism"), mean = 0, sd = sd(males) * 1/8 * effect_size) 
  
    # high dimorphism
    females[dimorphism_df[ , w] == "high dimorphism"] <- females[dimorphism_df[ , w] == "high dimorphism"] + rnorm(n = sum(dimorphism_df[ , w] == "high dimorphism"), mean = 0, sd(males) * 1/8 * (effect_size^2))
     # plot(males_vars[[w]], females)
     
    vars[[w]] <- c(males_vars[[w]], females)
        }  
    
    
    # put morpho variables together in a data frame
  morph_dat <- as.data.frame(do.call(cbind, vars))
  
  pca <- prcomp(morph_dat, scale. = TRUE)
  
  sp_labs <- paste0(rep("sp", spp), 1:spp)
  sp_labs <- c(sp_labs, sp_labs)
  
  # regular euclidean distance
  dists <- sapply(unique(sp_labs), function(x)  stats::dist(pca$x[sp_labs == x, ]))
  dist_df <- data.frame(dimorphism = dimorphism, dists = dists)
  dist_df$dimorphism <- factor(dist_df$dimorphism, levels = c("monomorphic", "low dimorphism", "high dimorphism"))
  dist_df$dimorphism <- ordered(dist_df$dimorphism)
  
  # weighted distances
  weight_dists <- sapply(unique(sp_labs), function(x)  weight_dist(X = pca$x[sp_labs == x, ], pca))
  weight_dist_df <- data.frame(dimorphism = dist_df$dimorphism, dists = weight_dists)
  
  # ratios
   ratios <- morph_dat[1:spp, ] / morph_dat[(spp + 1):(spp * 2), ]
  
  # check ratios
  # apply(ratios, 2, range)
  # apply(abs(ratios -1), 2, range)

  ratio_pca <- prcomp(abs(ratios - 1), scale. = FALSE)
  ratio_df <- data.frame(dimorphism = dist_df$dimorphism, dists = scale(ratio_pca$x[, 1]))
  
  return(list(morph_dat = morph_dat, dist_df = dist_df, weight_dist_df = weight_dist_df, ratio_df = ratio_df, pcas = list(raw_vars = summary(pca), ratios = summary(ratio_pca))))
  
}

# replicate model function
# model_call <- call("sim_pca_dimorph", dimorphic_vars = 7, covar_vars = 7, varying_and_dimorphic = FALSE)

repicate_models <- function(model, n = 1000){
  ps <- replicate(n = n, simplify = FALSE, expr = {
  
      # eval model
    dimorphs <- try(eval(model), silent = TRUE)
    
    if (!is(dimorphs, "try-error")){
    # run models
    lm_regular <- lm(dists ~ dimorphism, data = dimorphs$dist_df)
    lm_weight <- lm(dists ~ dimorphism, data = dimorphs$weight_dist_df)
    lm_ratio <- lm(dists ~ dimorphism, data = dimorphs$ratio_df)
    out <- list(regular = dimorphs$dist_df, weighted = dimorphs$weight_dist_df, ratios = dimorphs$ratio_df)
    
    
    mean_sds <- lapply(1:length(out), function(x) {
      mns <- stats::aggregate(formula = dists ~ dimorphism, data = out[[x]], FUN = base::mean)
      mns$sd <- stats::aggregate(formula = dists ~ dimorphism, data = out[[x]], FUN = stats::sd)[, 2]
      mns$data <- names(out)[x]
      
      return(mns)
    })
    
    mean_sds <- do.call(rbind, mean_sds)
    
    mean_sds$p_value <-  c(lmp(lm_regular), NA, NA, 
                           lmp(lm_weight), NA, NA, 
                           lmp(lm_ratio), NA, NA)
    } else
      mean_sds <- NULL  
      
    return(mean_sds)
  })
  
  ps <- do.call(rbind, ps)
  
  p_regular <- sum(ps$p_value[ps$data == "regular" & !is.na(ps$p_value)] < 0.05) / sum(ps$data == "regular" & !is.na(ps$p_value))
  
  p_weight <- sum(ps$p_value[ps$data == "weighted" & !is.na(ps$p_value)] < 0.05) / sum(ps$data == "weighted" & !is.na(ps$p_value))
  
  p_ratio <- sum(ps$p_value[ps$data == "ratios" & !is.na(ps$p_value)] < 0.05) / sum(ps$data == "ratios" & !is.na(ps$p_value))

  mean_sd <- aggregate(cbind(dists, sd) ~ dimorphism + data, data = ps, FUN = mean)
  
  return(list(ps = data.frame(p_regular = p_regular, p_weight = p_weight, p_ratio = p_ratio), mean_sd = mean_sd))  
}


#######
Exploring simulated data
Effect size = 2
This is an example data set generated by the simulation. In this case we used 10 dimorphic variables, 8 covarying variables and an effect size of 2

red = high dimorphism species
yellow = low dimporhism species
gray = monomorphic species
Male-female correlation for each simulated variable:
#######

# BEGIN CODE
set.seed(54)
dimorphs <- sim_pca_dimorph(dimorphic_vars = 10, covar_vars = 8, varying_and_dimorphic = "random", effect_size = 2)

cols <- as.character(rep(dimorphs$dist_df$dimorphism, 2))
cols[cols == "monomorphic"] <- "gray"
cols[cols == "low dimorphism"] <- "orange"
cols[cols == "high dimorphism"] <- "red"

par(mfrow = c(5, 3))
# par(mfrow = c(14, 2))

for (i in 1:ncol(dimorphs$morph_dat))
# for (i in 1:4)
  plot(c(dimorphs$morph_dat[1:(nrow(dimorphs$morph_dat) / 2), i]), c(dimorphs$morph_dat[((nrow(dimorphs$morph_dat) / 2) + 1):((nrow(dimorphs$morph_dat) / 2) * 2), i]), xlab = paste("Variable", i, "males"), ylab = paste("Variable", i, "females"), col = cols)

par(mfrow = c(1, 1))


### PCA variance explained
# On raw morphological data:

as.data.frame(dimorphs$pcas$raw_vars$importance)

as.data.frame(dimorphs$pcas$ratios$importance)



###### Effect size = 3
# Male-female correlation for each simulated variable:

set.seed(54)
dimorphs <- sim_pca_dimorph(dimorphic_vars = 10, covar_vars = 8, varying_and_dimorphic = "random", effect_size = 3)
# as.data.frame(dimorphs$pcas$raw_vars$importance)

cols <- as.character(rep(dimorphs$dist_df$dimorphism, 2))
cols[cols == "monomorphic"] <- "gray"
cols[cols == "low dimorphism"] <- "orange"
cols[cols == "high dimorphism"] <- "red"

par(mfrow = c(5, 3))
# par(mfrow = c(2, 2))

for (i in 1:ncol(dimorphs$morph_dat))
# for (i in 1:4)
  plot(c(dimorphs$morph_dat[1:(nrow(dimorphs$morph_dat) / 2), i]), c(dimorphs$morph_dat[((nrow(dimorphs$morph_dat) / 2) + 1):((nrow(dimorphs$morph_dat) / 2) * 2), i]), xlab = paste("Variable", i, "males"), ylab = paste("Variable", i, "females"), col = cols)

par(mfrow = c(1, 1))
######








################
Simulations
Simulations were run for all possible combinations of the following parameters:
0, 3, and 10 dimorphic variables
3, 8 and 10 covarying variables
3 effect size: low = 2, moderate = 3 and high = 4
dimorphic variables are: 1) among the most varying variables 2) among the least varying variables or 3) randomly selected with regard to the variance of the morphological variables (no assocation between dimorphism and variance)
dimorphic variables are 1) the same across dimorphic species or 2) different among dimorphic species (i.e. randomly selected)
A total of 162 different simulation, each one replicated 10000 times
################

grid <- expand.grid(dimorphic_vars = c(10, 3, 0), 
                    covar_vars = c(10, 8, 3), 
                    varying_and_dimorphic = c("yes", "no", "random"), 
                    same_dimorphic = c(TRUE, FALSE), 
                    effect_size = c(2, 3, 4))

out <- pblapply(1:nrow(grid), cl = 18, function(y){
  reps <- repicate_models(model = call("sim_pca_dimorph", dimorphic_vars = grid$dimorphic_vars[y], covar_vars = grid$covar_vars[y], varying_and_dimorphic = grid$varying_and_dimorphic[y], same_dimorphic = grid$same_dimorphic[y], effect_size = grid$effect_size[y]), n = 10000)  

reps$ps$dimorphic_vars <- grid$dimorphic_vars[y]
reps$ps$covar_vars <- grid$covar_vars[y]
reps$ps$varying_and_dimorphic = grid$varying_and_dimorphic[y]
reps$ps$same_dimorphic = if(grid$same_dimorphic[y]) "yes" else "no"
reps$ps$effect_size = if(grid$effect_size[y] == 2) "low" else if(grid$effect_size[y] == 3) "moderate" else "high"
  
return(reps)
})

saveRDS(out, "./output/10000_replicates_morphological_dimorphism.RDS")









############
Results
The following graph shows the proportion of simulation iterations in which significant differences between males and females were detected (i.e. statistically significant dimorphism) for the 3 metrics (x axis), 3 levels of dimorphism (columns) and number of dimoprhic variables (rows). 
Note that 0 dimorphic variables (first row) means no dimorphism:
############
out <- readRDS("./output/10000_replicates_morphological_dimorphism.RDS")

ps <- lapply(out, "[[", 1)

ps <- do.call(rbind, ps)

agg_p_val_regular <- aggregate(p_regular ~ effect_size  + dimorphic_vars, data = ps, FUN = mean)
agg_p_val_regular$metric <- "regular PCA"

agg_p_val_weight <- aggregate(p_weight ~ effect_size  + dimorphic_vars, data = ps, FUN = mean)
agg_p_val_weight$metric <- "weighted PCA"

agg_p_val_ratio <- aggregate(p_ratio ~ effect_size  + dimorphic_vars, data = ps, FUN = mean)
agg_p_val_ratio$metric <- "ratios"

names(agg_p_val_ratio)[3] <- names(agg_p_val_weight)[3] <- names(agg_p_val_regular)[3] <- "p_vals"

agg_p_vals <- rbind(agg_p_val_regular, agg_p_val_weight, agg_p_val_ratio)

agg_p_vals$effect_size <- factor(agg_p_vals$effect_size, levels = c("low", "moderate", "high"))

ggplot(agg_p_vals, aes(x = metric, y = p_vals, fill = metric)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  scale_fill_viridis_d(begin = 0.3, end = 0.8) +
  facet_grid(dimorphic_vars ~ effect_size) +
  theme_classic(base_size = 20) +
  labs(y = "Proportion of significant iterations", x = "Dimorphism metric") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 





###########
Takeaways
Both regular PCA and weigthed PCA metrics show high performance for detecting dimorphism when it is found (second and third rows in the multipanel graph) or an error rate close to 0.05 when no dimorphism should be found

Male-to-female ratios did have a good error rate but much lower performance in detecting dimorphism
###########



# Session information

## R version 4.1.0 (2021-05-18)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 20.04.2 LTS
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/atlas/libblas.so.3.10.3
## LAPACK: /usr/lib/x86_64-linux-gnu/atlas/liblapack.so.3.10.3
## 
## locale:
##  [1] LC_CTYPE=pt_BR.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=es_CR.UTF-8        LC_COLLATE=pt_BR.UTF-8    
##  [5] LC_MONETARY=es_CR.UTF-8    LC_MESSAGES=pt_BR.UTF-8   
##  [7] LC_PAPER=es_CR.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=es_CR.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
## [1] ggplot2_3.3.5    DT_0.18          NormalGamma_1.1  histogram_0.0-25
## [5] optimx_2021-6.12 knitr_1.38       pbapply_1.5-0   
## 
## loaded via a namespace (and not attached):
##  [1] tidyselect_1.1.1    xfun_0.30           bslib_0.2.5.1      
##  [4] purrr_0.3.4         colorspace_2.0-2    vctrs_0.3.8        
##  [7] generics_0.1.0      viridisLite_0.4.0   htmltools_0.5.2    
## [10] yaml_2.3.5          utf8_1.2.2          rlang_1.0.1        
## [13] jquerylib_0.1.4     pillar_1.6.4        glue_1.6.2         
## [16] withr_2.4.3         DBI_1.1.1           lifecycle_1.0.1    
## [19] stringr_1.4.0       munsell_0.5.0       gtable_0.3.0       
## [22] htmlwidgets_1.5.3   evaluate_0.15       labeling_0.4.2     
## [25] fastmap_1.1.0       fansi_1.0.0         highr_0.9          
## [28] scales_1.1.1        jsonlite_1.8.0      farver_2.1.0       
## [31] digest_0.6.29       stringi_1.7.6       dplyr_1.0.7        
## [34] numDeriv_2016.8-1.1 grid_4.1.0          cli_3.1.0          
## [37] tools_4.1.0         magrittr_2.0.3      sass_0.4.0         
## [40] tibble_3.1.6        crayon_1.5.0        pkgconfig_2.0.3    
## [43] ellipsis_0.3.2      assertthat_0.2.1    rmarkdown_2.13     
## [46] rstudioapi_0.13     R6_2.5.1            compiler_4.1.0


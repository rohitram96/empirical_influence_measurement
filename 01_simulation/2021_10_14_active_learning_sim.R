##### Imports

library(future)
library(furrr)
library(truncnorm)
library(dplyr)
library(purrr)
library(BradleyTerry2)
library(ggplot2)
library(stackoverflow)
library(rstack)
library(poweRlaw)

set.seed(17)

##### Fitted Data from External

mean_worker_accuracy = 0.5380250775976005
std_worker_accuracy = 0.20513411249363095

empirical_hits_data <- c(read.csv('./empirical_hits_dist_5.csv', header = F)$V1)*10

##### Generate user influence 
# From a powerlaw dist

generate_user_influence <- function(n, alpha = 2.016, mmin = 1) {
  if ( !exists('user_infl') || is.null(user_infl) ) {
    user_infl <<- conpl$new()
  }
  user_infl$setXmin(mmin)
  user_infl$setPars(alpha)
  return(dist_rand(user_infl, n))
}

##### Simulate a QSORT run
# Outputs the correlation of the BTscores with the original influence intensities

get_data <- function(args, theta){
  e <- new.env()
  e$efficacy_mean <- args[1]
  e$efficacy_std <- args[2]
  e$budget <- args[3]
  e$ntargets <- args[4]
  e$noise <- args[5]
  
  # We make a vector of size budget, sampling from the empirical hits dist
  e$workers <- c()
  n_workers = 0
  while (length(e$workers) < e$budget) {
    n_workers = n_workers + 1
    e$workers <- c(e$workers, rep(n_workers,quantile(empirical_hits_data,runif(1), FALSE, TRUE, 1)))
  }
  e$workers <- sample(e$workers[1:e$budget])
  
  # Each worker has a certain accuracy, here dubbed efficacy, which we assume is normally distributed (truncated)
  worker_efficacy = rtruncnorm(n=n_workers,a=0.0, b=1, mean=e$efficacy_mean, sd=e$efficacy_std)
  
  # p function representes the probability that we identify i < j, given the ground truth influence intensities
  # The probability of guessing wrong is low for theta that are different, 
  p <- function(i,j){
    return(1/(1+exp(-1*(theta[j] - theta[i]))))
  }
  
  # There are a number of ways this response can be generated. We use A here.
  # A: With worker_efficacy probability, the worker makes an educated guess, else a random guess.
  #       If the worker makes an educated guess, the worker is correct with probabiliy p(i,j), else wrong.
  # B: The worker is correct with probability worker_efficacy
  # C: The worker is correct with probability p(i,j)
  # D: The worker produces the ground_truth always
  
  getResponse = function(worker_id,i,j){
    # A
    # if (runif(1) > worker_efficacy[worker_id]) {
    #   return(runif(1) > 0.5)
    # }
    # return(runif(1) >= p(i,j))
    # B
    # if (runif(1) > worker_efficacy[worker_id]) {
    #   return(theta[i]<theta[j])
    # }
    # return(theta[j]<theta[i])
    # C
    # return(runif(1) >= p(i,j))
    return(runif(1) < p(i,j)) # large positive diff implies high prob of true
    # D
    # return(theta[i]<theta[j])
    
    # # This logic here flips the comparison some of the time
    # return(!xor((1/(1+exp(-1*(theta[j] - theta[i]))) > runif(1)),theta[i]<theta[j]))
    # return(theta[i]<theta[j]) # We can use this in the case that we don't want to look at closeness of intensities
  }
  
  # We set up some variables which are consistent between calls of compare
  e$data <- c()
  e$worker_count <- 0
  e$worker_id <- 0
  
  # compare simulates a worker comparison using the current worker
  compare = function(i,j){
    if (e$budget > 0 && e$worker_id <= n_workers+1){
      e$worker_count <- e$worker_count + 1
      e$worker_id <- e$workers[e$worker_count]
      e$budget <- e$budget - 1
      response <- getResponse(e$worker_id,i,j)
      # e$data <- rbind(e$data, if (response) c(i, j) else c(j, i) )
      e$data <- rbind(e$data, if (response) c(j, i) else c(i, j) )
      return(response)
    }else{
      stop("Ran out!")
    }
  }
  
  # runs quicksort using the simulation compare rather than a direct comparison. 
  # We do this iteratively rather than recursively so that we don't exceed R default stack limit.
  qsort_i = function(arr){
    m <- new.env()
    m$arr <- arr
    swap <- function(i,j){
      temp <- m$arr[i]
      m$arr[i] <- m$arr[j]
      m$arr[j] <- temp
    }
    partition <- function(start, end){
      pivot <- m$arr[end]
      i <- start - 1
      for (j in start:(end-1)) {
        # if (compare(m$arr[j],pivot)){
        if (!compare(m$arr[j],pivot)){
          # if (m$arr[j] < pivot){
          i = i + 1
          swap(i,j)
        }
      }
      swap(i+1,end)
      return(i+1)
    }
    qsort_l <- function(low,high){
      n <- new.env()
      n$top <- -1
      inc <- function(){
        temp <- n$top
        n$top <- n$top + 1
        return(temp)
      }
      
      n$stack <- stack$new()
      n$stack$push(low)
      n$stack$push(high)
      
      while (!n$stack$is_empty()) {
        end = n$stack$pop()
        start =n$stack$pop()
        pi = partition(start,end)
        if (pi -1 > start){
          n$stack$push(start)
          n$stack$push(pi - 1)
        }
        if (pi + 1 < end){
          n$stack$push(pi + 1)
          n$stack$push(end)
        }
      }
      
    }
    tryCatch({
      qsort_l(1,length(m$arr))
    },error = function(c){
      return(m$arr)
    })
    return(m$arr)
  }
  
  # We run quicksort as many times as our budget allows us to, and shuffle the users to sort each time.
  while (e$budget > 0 && e$worker_count < length(e$workers)) {
    qsort_i(sample(1:(e$ntargets)))
  }
  
  return(e$data)
}

#
compute_bradle_terry_scores <- function(data){
  countsToBinomial_ <- function(xtab) {
    ## make square if necessary
    if (nrow(xtab) != ncol(xtab) || !all(rownames(xtab) == colnames(xtab))) {
      dat <- as.data.frame(xtab)
      lev <- union(rownames(xtab), colnames(xtab))
      dat[,1] <- factor(dat[,1], levels = lev)
      dat[,2] <- factor(dat[,2], levels = lev)
      xtab <- tapply(dat[,3], dat[1:2], sum)
      xtab[is.na(xtab)] <- 0
    }
    ##assumes square
    players <- rownames(xtab)
    comb <- t(combn(nrow(xtab), 2))
    won <- xtab[comb]
    lost <- t(xtab)[comb]
    res <- !(won == 0 & lost == 0)
    player1 <- factor(players[comb[,1]], levels = players)[res]
    player2 <- factor(players[comb[,2]], levels = players)[res]
    data.frame(player1, player2, win1 = won[res], win2 = lost[res])
  }
  
  
  tab <- table(data[,1], data[,2])
  tb <- countsToBinomial_(tab)
  names(tb)[1:2] <- c("target1", "target2")
  bt <- BTm(cbind(win1, win2), target1, target2, data = tb)
  theta_hat <- c(c('..1'=0),bt$coefficients)
  return(theta_hat)
}

compute_correlations <- function(theta,theta_hat){
  correlations = c('pearson'=cor(theta, theta_hat, method="pearson"), 'kendall'=cor(theta, theta_hat, method="kendall"), 'spearman'=cor(theta, theta_hat, method="spearman"), 'displacement'=sum(abs(order(theta_hat) - order(theta))))
  return(correlations)
}

# theta represents the underlying intensity of target users
get_theta <- function(ntargets, noise){
  # return(runif(ntargets, min=0, max=(ntargets/noise))) # If theta is uniformly distributed we use this
  return((generate_user_influence(ntargets) - 1)/noise + 1) # If theta is powerlaw distributed we use this
}

# A helper function to attach args to the correlation results
# compute_metrics <- function(args){
#   theta <- get_theta(ntargets = args['ntargets'], noise = args['noise'])
#   data <- get_data(args[1:4], theta = theta)
#   results <- c()
#   for (i in seq(args['min_budget'], nrow(data), args['inc_budget'])){
#     data_t <- data[1:i,c(1,2)]
#     theta_hat <- compute_bradle_terry_scores(data_t)
#     args['max_budget'] <- i
#     results <- rbind(results,c(args, compute_correlations(theta = theta, theta_hat = theta_hat)))
#   }
#   return(results)
# }

###### Setup Future Configuration
# plan(cluster, workers=c(rep('jup2', 48),rep('jup11', 36)), homogeneous=FALSE)
# plan(cluster, workers=c(rep('jupiter2', 48), rep('jupiter12', 36)), homogeneous=FALSE)
nodes <- c(rep('jupiter2', 4), rep('jupiter4', 4))
# plan(list(tweak(remote, workers = nodes), tweak(multisession, workers = 12)), homogenous=FALSE)
# plan(remote, workers=c('heph0'))
# plan(multicore)
# plan(sequential)

compute_metrics <- function(args){
  theta <- get_theta(ntargets = args['ntargets'], noise = args['noise'])
  data <- get_data(args[1:4], theta = theta)
  wrapper <- function(i){
    data_t <- data[1:i,c(1,2)]
    theta_hat <- compute_bradle_terry_scores(data_t)
    args['max_budget'] <- i
    print(i)
    return(c(args, compute_correlations(theta = theta, theta_hat = theta_hat)))
  }
  results <- future_map(seq(args['min_budget'], nrow(data), args['inc_budget']), ~wrapper(.x))
  return(as.data.frame(do.call(rbind, results)))
}


####### A basic example of a single run
# worker_efficacy_mean <- mean_worker_accuracy
# worker_efficacy_std <- std_worker_accuracy
# ntargets <- 500
# max_budget <- 30000
# noise <- 1.0
# min_budget <- 10000
# inc_budget <- 500
# 
# results <- c()
# for (rep in 1:1){
#   print(rep)
#   tmp = compute_metrics(c(worker_efficacy=worker_efficacy_mean,worker_efficacy_std=worker_efficacy_std,max_budget=max_budget, ntargets=ntargets, noise=noise, min_budget=min_budget, inc_budget=inc_budget))
#   results <- rbind(results, tmp)
# }

# saveRDS(results, './results/basic_a.rds')

##### Vary Worker Efficacy for strategy b

# worker_efficacy_mean <- mean_worker_accuracy
# worker_efficacy_std <- std_worker_accuracy
# ntargets <- 500
# max_budget <- 30000
# noise <- 1.0
# min_budget <- 30000
# inc_budget <- 500
# 
# results <- c()
# for (rep in 1:50){
#   print(rep)
#   wrapper_b <- function(worker_efficacy_mean){
#     print(worker_efficacy_mean)
#     return(compute_metrics(c(worker_efficacy=worker_efficacy_mean,worker_efficacy_std=worker_efficacy_std,max_budget=max_budget, ntargets=ntargets, noise=noise, min_budget=min_budget, inc_budget=inc_budget)))
#   }
#   tmp <- future_map(seq(0.05,1.0,0.05), ~wrapper_b(.x))
#   tmp <- as.data.frame(do.call(rbind, tmp))
#   results <- rbind(results, tmp)
# }
# 
# saveRDS(results, './results/worker_efficacy_experiment_b.rds')
# worker_efficacy_experiment_b <- readRDS("~/simulation/results/worker_efficacy_experiment_b.rds")
# 
# worker_efficacy_experiment_b %>%
#   group_by(worker_efficacy) %>%
#   summarise(spearman = mean(spearman)) %>%
#   ggplot(mapping=aes(x=worker_efficacy, y=spearman)) +
#   geom_point() +
#   geom_smooth()


##### Vary Worker Efficacy for strategy a

# worker_efficacy_mean <- mean_worker_accuracy
# worker_efficacy_std <- std_worker_accuracy
# ntargets <- 500
# max_budget <- 30000
# noise <- 1.0
# min_budget <- 30000
# inc_budget <- 500
# 
# results <- c()
# for (rep in 1:50){
#   print(rep)
#   wrapper_b <- function(worker_efficacy_mean){
#     print(worker_efficacy_mean)
#     return(compute_metrics(c(worker_efficacy=worker_efficacy_mean,worker_efficacy_std=worker_efficacy_std,max_budget=max_budget, ntargets=ntargets, noise=noise, min_budget=min_budget, inc_budget=inc_budget)))
#   }
#   tmp <- future_map(seq(0.05,1.0,0.05), ~wrapper_b(.x))
#   tmp <- as.data.frame(do.call(rbind, tmp))
#   results <- rbind(results, tmp)
# }
# 
# saveRDS(results, './results/worker_efficacy_experiment_a_2.rds')
# worker_efficacy_experiment_a <- readRDS("~/simulation/results/worker_efficacy_experiment_a_2.rds")
# 
# worker_efficacy_experiment_a %>%
#   group_by(worker_efficacy) %>%
#   summarise(spearman = mean(spearman)) %>%
#   ggplot(mapping=aes(x=worker_efficacy, y=spearman)) +
#   geom_point() +
#   geom_smooth()


##### Show correlation against worker efficacy for different strategies

# worker_efficacy_experiment_a <- readRDS("~/simulation/results/worker_efficacy_experiment_a_2.rds")
# 
# p_a <- worker_efficacy_experiment_a %>%
#   group_by(worker_efficacy) %>%
#   summarise(spearman = mean(spearman)) %>%
#   mutate(strategy='a')
# 
# worker_efficacy_experiment_b <- readRDS("~/simulation/results/worker_efficacy_experiment_b.rds")
# 
# p_b <- worker_efficacy_experiment_b %>%
#   group_by(worker_efficacy) %>%
#   summarise(spearman = mean(spearman)) %>%
#   mutate(strategy='b')
# 
# bind_rows(p_a,p_b) %>%
#   ggplot(mapping=aes(x=worker_efficacy, y=spearman, color=strategy)) +
#   geom_line()

##### Vary ntargets and budget
plan(multisession)
worker_efficacy_mean <- mean_worker_accuracy
worker_efficacy_std <- std_worker_accuracy
# ntargets <- 500
max_budget <- 40000
noise <- 1.756769
min_budget <- 30000
inc_budget <- 1000

results <- c()
for (rep in 1:5){
  print(rep)
  wrapper_b <- function(ntargets){
    print(ntargets)
    return(compute_metrics(c(worker_efficacy=worker_efficacy_mean,worker_efficacy_std=worker_efficacy_std,max_budget=max_budget, ntargets=ntargets, noise=noise, min_budget=min_budget, inc_budget=inc_budget)))
  }
  tmp <- future_map(seq(300,700,100), ~wrapper_b(.x))
  tmp <- as.data.frame(do.call(rbind, tmp))
  results <- rbind(results, tmp)
  saveRDS(results, './results_extended/2021_10_14_ntarget_budget_experiment_5.rds')
}

# nodes <- c(rep('jupiter2', 4), rep('jupiter4', 4))
# plan(list(tweak(remote, workers = nodes), tweak(multisession, workers = 12)), homogenous=FALSE)
# plan(multisession)
# 
# specific_point <- compute_metrics(c(worker_efficacy=worker_efficacy_mean,worker_efficacy_std=worker_efficacy_std,max_budget=40000, ntargets=500, noise=noise, min_budget=30000, inc_budget=2000))
# saveRDS(results, './results_extended/2021_10_14_specific.rds')
# specific_point %>% select(max_budget, ntargets, spearman)



ntarget_budget_experiment_c <- readRDS("./results/2020_11_11_ntarget_budget_experiment_c_3.rds")

ntarget_budget_experiment_c %>%
  mutate(dollars = max_budget * 0.04 / 10) %>%
  group_by(dollars,ntargets) %>%
  summarise(spearman = mean(spearman)) %>%
  ggplot() +
  geom_contour_filled(mapping=aes(x=dollars, y=ntargets, z = spearman), binwidth = 0.01)


library(metR)
library(ggforce)
craftbrewer_pal <- function (type = "seq", palette = 1, direction = 1) 
{
  pal <- scales:::pal_name(palette, type)
  force(direction)
  function(n) {
    n_max_palette <- RColorBrewer:::maxcolors[names(RColorBrewer:::maxcolors) == palette]
    
    if (n < 3) {
      pal <- suppressWarnings(RColorBrewer::brewer.pal(n, pal))
    } else if (n > n_max_palette){
      rlang::warn(paste(n, "colours used, but", palette, "has only",
                        n_max_palette, "- New palette created based on all colors of", 
                        palette))
      n_palette <- RColorBrewer::brewer.pal(n_max_palette, palette)
      colfunc <- grDevices::colorRampPalette(n_palette)
      pal <- colfunc(n)
    }
    else {
      pal <- RColorBrewer::brewer.pal(n, pal)
    }
    pal <- pal[seq_len(n)]
    if (direction == -1) {
      pal <- rev(pal)
    }
    pal
  }
}
scale_fill_craftfermenter <- function(..., type = "seq", palette = 1, direction = -1, na.value = "grey50", guide = "coloursteps", aesthetics = "fill") {
  type <- match.arg(type, c("seq", "div", "qual"))
  if (type == "qual") {
    warn("Using a discrete colour palette in a binned scale.\n  Consider using type = \"seq\" or type = \"div\" instead")
  }
  binned_scale(aesthetics, "fermenter", ggplot2:::binned_pal(craftbrewer_pal(type, palette, direction)), na.value = na.value, guide = guide, ...)
}
ntarget_budget_experiment_c %>%
  group_by(max_budget,ntargets) %>%
  summarise(spearman = mean(spearman)) %>%
  ggplot() +
  metR::geom_contour_fill(aes(x=max_budget, y=ntargets, z = spearman)) +
  scale_fill_craftfermenter(
    breaks = seq(0.5, 1.0, 0.05), 
    palette = "RdYlGn", 
    direction = 1,
    limits = c(-2,11),
    guide = guide_colorsteps(
      frame.colour = "black", 
      ticks.colour = "black", # you can also remove the ticks with NA
      barwidth=20)
  ) +
  annotate(geom = "point", x = 25000, y = 500, shape='cross', size=5) +
  # geom_segment(x = 25000, xend=25000, y=470, yend=530) +
  # geom_segment(x = 24500, xend=25500, y=500, yend=500) +
  # geom_vline(xintercept = 25000) +
  # geom_hline(yintercept = 500) +
  theme(legend.position = "bottom")

ntarget_budget_experiment_c %>%
  group_by(max_budget,ntargets) %>%
  summarise(spearman = mean(spearman)) %>%
  filter(max_budget == 25000, ntargets == 500)

ntarget_budget_experiment_c %>%
  mutate(HITS = max_budget / 10) %>%
  group_by(HITS,ntargets) %>%
  summarise(spearman = mean(spearman)) %>%
  ggplot() +
  geom_contour_filled(mapping=aes(x=HITS, y=ntargets, z = spearman), binwidth = 0.01)

ntarget_budget_experiment_c %>%
  mutate(expected_runs = max_budget / (2*ntargets*log(ntargets))) %>%
  group_by(expected_runs,ntargets) %>%
  summarise(spearman = mean(spearman)) %>%
  ggplot() +
  geom_point(mapping=aes(x=expected_runs, y=ntargets, colour = spearman))

ntarget_budget_experiment_c %>%
  filter(max_budget == 24000, ntargets == 500) %>%
  ggplot(aes(x=spearman)) +
  geom_boxplot() +
  geom_density()

a <- ntarget_budget_experiment_c %>% filter(max_budget == 24000, ntargets == 500) %>% summarise(spearman = mean(spearman))
s <- ntarget_budget_experiment_c %>% filter(max_budget == 24000, ntargets == 500) %>% summarise(spearman = sd(spearman))
n <- ntarget_budget_experiment_c %>% filter(max_budget == 24000, ntargets == 500) %>% count()
error <- qnorm(0.975)*s/sqrt(n)
left <- a-error
right <- a+error
list(left = left, right = right)
##### Vary Noise for strategy c
# 
# worker_efficacy_mean <- 1.0
# worker_efficacy_std <- 1.0
# ntargets <- 500
# max_budget <- 30000
# # noise <- 1.0
# min_budget <- 10000
# inc_budget <- 10000
# 
# results <- c()
# for (rep in 1:20){
#   print(rep)
#   wrapper_b <- function(noise){
#     print(noise)
#     return(compute_metrics(c(worker_efficacy=worker_efficacy_mean,worker_efficacy_std=worker_efficacy_std,max_budget=max_budget, ntargets=ntargets, noise=noise, min_budget=min_budget, inc_budget=inc_budget)))
#   }
#   tmp <- future_map(c(seq(0.1,3.0,0.1),1.573017), ~wrapper_b(.x))
#   tmp <- as.data.frame(do.call(rbind, tmp))
#   results <- rbind(results, tmp)
# }
# 
# saveRDS(results, './results/noise_experiment_c_4.rds')
# noise_experiment_c <- readRDS("./results/noise_experiment_c_4.rds")
# 
# noise_experiment_c %>%
#   group_by(noise) %>%
#   summarise(spearman = mean(spearman)) %>%
#   ggplot(mapping=aes(x=noise, y=spearman)) +
#   geom_point() +
#   geom_smooth()
# 
# noise_experiment_c %>%
#   ggplot(mapping=aes(x=noise, y=spearman)) +
#   geom_point() +
#   geom_smooth() +
#   facet_wrap(~max_budget)

##### Vary Budget for strategy c with fixed noise
# 
# worker_efficacy_mean <- 1.0
# worker_efficacy_std <- 1.0
# # ntargets <- 500
# max_budget <- 30000
# # noise <- 1.573017
# noise <- 0.5
# min_budget <- 30000
# inc_budget <- 1000
# 
# results <- c()
# for (rep in 1:1){
#   print(rep)
#   wrapper_b <- function(ntargets){
#     print(ntargets)
#     return(compute_metrics(c(worker_efficacy=worker_efficacy_mean,worker_efficacy_std=worker_efficacy_std,max_budget=max_budget, ntargets=ntargets, noise=noise, min_budget=min_budget, inc_budget=inc_budget)))
#   }
#   tmp <- future_map(seq(500,500,100), ~wrapper_b(.x))
#   tmp <- as.data.frame(do.call(rbind, tmp))
#   results <- rbind(results, tmp)
# }

# saveRDS(results, './results/2020_11_10_ntarget_budget_experiment_c_fixed_noise.rds')

##### Show correlation against worker efficacy for different strategies

# budget_experiment_c <- readRDS("./results/2020_11_09_ntarget_budget_experiment_c_fixed_noise.rds")
# 
# budget_experiment_c %>%
#   group_by(max_budget, ntargets) %>%
#   summarise(spearman = mean(spearman)) %>%
#   ggplot(mapping=aes(x=max_budget, y=spearman, colour=ntargets)) +
#   geom_point() +
#   geom_smooth()


###### Show facet grid

# data %>%
#   ggplot(mapping = aes(group=worker_efficacy, y=spearman)) +
#   geom_boxplot() +
#   facet_grid(budget~ntargets) +
#   labs(
#     x = "ntargets",
#     y = "Budget"
#   )

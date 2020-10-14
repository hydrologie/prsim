prsimLD.wave <- function (data, station_id = "Qobs", number_sim = 1, win_h_length = 15,
                          marginal = c("kappa", "empirical"), n_par = 4, n_wave = 100,
                          marginalpar = TRUE, GoFtest = NULL, verbose = TRUE, suppWarn = FALSE, proba=0.9999,
                          ...)
{
  fun_icwt <- function(x) {
    wt.r <- Re(x)
    J <- length(x[1, ])
    dial <- 2 * 2^(0:J * 0.125)
    rec <- rep(NA, (length(x[, 1])))
    for (l in 1:(length(x[, 1]))) {
      rec[l] <- 0.2144548 * sum(wt.r[l, ]/sqrt(dial)[1:length(wt.r[l,
                                                                   ])])
    }
    return(rec)
  }
  rand.kappaLD <- function(numerosita, xi, alfa, k, h) {
    F <- runif(numerosita, min = 1e-10, max = proba)
    x <- invF.kappa(F, xi, alfa, k, h)
    #cat(paste0(" LD function F=  ", F, "\n"))
    return(x)
  }
  if (!is.null(GoFtest)) {
    GoFtest <- toupper(GoFtest)[1]
    if (!(GoFtest %in% c("AD", "KS")))
      stop("'GoFtest' should be either 'NULL', 'AD' or 'KS'.")
  }
  else GoFtest <- "NULL"
  marginal <- marginal[1]
  if (!(marginal %in% c("kappa", "empirical"))) {
    if (!is.character(marginal))
      stop("'marginal' should be a character string.")
    rCDF <- get(paste0("r", marginal), mode = "function")
    CDF_fit <- get(paste0(marginal, "_fit"), mode = "function")
    if (GoFtest == "AD")
      pCDF <- get(paste0("p", marginal), mode = "function")
  }
  op <- options("warn")$warn
  for (l in 1:length(data)) {
    if (nrow(data[[l]])[1] < 730)
      stop("At least one year of data required.")
    if (is.numeric(station_id)) {
      station_id <- colnames(data[[l]])[station_id]
    }
    if (is.na(station_id) || !("Qobs" %in% colnames(data[[l]])))
      stop("Wrong column (name) for observations selected.")
    if (any(class(data[[l]][, 1]) %in% c("POSIXct", "POSIXt"))) {
      data <- data.frame(YYYY = as.integer(format(data[[l]][,
                                                            1], "%Y")), MM = as.integer(format(data[[l]][,
                                                                                                         1], "%m")), DD = as.integer(format(data[[l]][,
                                                                                                                                                      1], "%d")), Qobs = data[[l]][, station_id],
                         timestamp = data[[l]][, 1])
    }
    else {
      if (!all(c("YYYY", "MM", "DD") %in% colnames(data[[l]])))
        stop("Wrong time column names")
      data[[l]] <- data[[l]][, c("YYYY", "MM", "DD", station_id)]
      tmp <- paste(data[[l]]$YYYY, data[[l]]$MM, data[[l]]$DD,
                   sep = " ")
      names(data[[l]]) <- c("YYYY", "MM", "DD", "Qobs")
      data[[l]]$timestamp <- as.POSIXct(strptime(tmp,
                                                 format = "%Y %m %d", tz = "GMT"))
    }
    data[[l]] <- data[[l]][format(data[[l]]$timestamp, "%m %d") !=
                             "02 29", ]
    if (which(format(data[[l]]$timestamp, format = "%j") ==
              "001")[1] > 1) {
      data[[l]] <- data[[l]][-c(1:(which(format(data[[l]]$timestamp,
                                                format = "%j") == "001")[1] - 1)), ]
    }
    if ((nrow(data[[l]])%%365) > 0)
      stop("No missing values allowed. Some days are missing.")
    if (length(which(is.na(data[[l]]$timestamp))) > 0) {
      data[[l]][which(is.na(data[[l]]$timestamp)), ]$Qobs <- mean(data[[l]]$Qobs,
                                                                  na.rm = T)
    }
    data[[l]]$index <- as.numeric(format(data[[l]]$timestamp,
                                         format = "%j"))
    if (length(which(is.na(data[[l]]$index)) > 0)) {
      data[[l]]$index[which(is.na(data[[l]]$index))] <- rep(c(1:365),
                                                            times = length(unique(data[[l]]$YYYY)))[which(is.na(data[[l]]$index))]
    }
  }
  if (verbose)
    cat(paste0("Detrending with (half-)length ", win_h_length,
               "...\n"))
  noise_mat_r <- list()
  for (r in 1:number_sim) {
    ts_wn <- rnorm(n = length(data[[1]]$Qobs), mean = 0,
                   sd = 1)
    wt_noise <- wavCWT(x = ts_wn, wavelet = "morlet", n.scale = n_wave)
    noise_mat_r[[r]] <- as.matrix(wt_noise)
  }
  par_day_list <- marginal_list <- list()
  for (l in 1:length(data)) {
    if (marginal == "empirical") {
      marginal_list[[l]] <- "empirical"
    }
    if (marginal == "kappa") {
      marginal_list[[l]] <- "kappa"
      p_vals <- numeric(365)
      par_day <- matrix(0, nrow = 365, ncol = 4)
      win_length <- c(1:win_h_length)
      for (d in c(1:365)) {
        before <- data[[l]]$index[d + 365 - win_length]
        after <- data[[l]]$index[d + 365 + win_length -
                                   1]
        ids <- c(before, after)
        data_window <- data[[l]]$Qobs[which(data[[l]]$index %in%
                                              ids)]
        ll <- homtest::Lmoments(data_window)
        if (suppWarn) {
          suppressWarnings(test <- try(par.kappa(ll[1],
                                                 ll[2], ll[4], ll[5]), silent = TRUE))
        }
        else {
          test <- try(par.kappa(ll[1], ll[2], ll[4],
                                ll[5]), silent = TRUE)
        }
        if (length(test) > 1) {
          kap_par <- test
          par_day[d, ] <- unlist(kap_par)
          quant <- sort(data_window)
          thresh <- kap_par$xi + kap_par$alfa * (1 -
                                                   kap_par$h^(-kap_par$k))/kap_par$k
          if (!is.na(thresh)) {
            quant <- quant[which(quant > thresh)]
          }
          data_kap <- rand.kappaLD(length(data_window),
                                   xi = kap_par$xi, alfa = kap_par$alfa, k = kap_par$k,
                                   h = kap_par$h)
          if (tolower(GoFtest) == "ks")
            p_vals[d] <- ks_test(data_window, data_kap)
          if (tolower(GoFtest) == "ad") {
            try_ad_test <- try(ad.test(data_window,
                                       F.kappa, xi = kap_par$xi, alfa = kap_par$alfa,
                                       k = kap_par$k, h = kap_par$h), silent = TRUE)
            if (length(try_ad_test) == 1) {
              p_vals[d] <- NA
            }
            else {
              p_vals[d] <- try_ad_test$p.value
            }
          }
        }
        else {
          if (d == 1) {
            p_vals[d] <- NA
            par_day[d, ] <- NA
          }
          else {
            p_vals[d] <- p_vals[d - 1]
            par_day[d, ] <- par_day[d - 1, ]
          }
        }
      }
      if (length(which(is.na(par_day[, 1]))) == 365) {
        marginal_list[[l]] <- "empirical"
      }
      else {
        if (length(which(is.na(par_day[, 1]))) > 0) {
          indices <- rev(which(is.na(par_day[, 1])))
          for (i in 1:length(indices)) {
            par_day[indices[i], ] <- par_day[indices[i] +
                                               1, ]
          }
        }
      }
      par_day_list[[l]] <- par_day
    }
    if (marginal != "kappa" & marginal != "empirical") {
      marginal_list[[l]] <- marginal
      p_vals <- numeric(365)
      par_day <- matrix(0, nrow = 365, ncol = n_par)
      for (d in c(1:365)) {
        win_length <- seq(1:15)
        before <- data[[l]]$index[d + 365 - win_length]
        after <- data[[l]]$index[d + 365 + win_length -
                                   1]
        ids <- c(before, after)
        data_window <- data[[l]]$Qobs[which(data[[l]]$index %in%
                                              ids)]
        theta <- CDF_fit(xdat = data_window, ...)
        data_random <- rCDF(n = length(data_window),
                            theta)
        if (tolower(GoFtest) == "ks") {
          p_vals[d] <- ks_test(data_window, data_random)
        }
        if (tolower(GoFtest) == "ad") {
          p_vals[d] <- ad.test(data_window, pCDF, theta)$p.value
        }
        par_day[d, ] <- theta
      }
      par_day_list[[l]] <- par_day
    }
  }
  for (l in 1:length(data)) {
    data[[l]]$norm <- data[[l]]$Qobs - mean(data[[l]]$Qobs,
                                            na.rm = T)
  }
  if (verbose)
    cat(paste0("Starting ", number_sim, " simulations:\n"))
  out_list <- list()
  for (l in 1:length(data)) {
    data_sim <- list()
    for (r in c(1:number_sim)) {
      wt_morlet <- wavCWT(x = data[[l]]$norm, wavelet = "morlet",
                          n.scale = n_wave)
      morlet_mat <- as.matrix(wt_morlet)
      modulus <- Mod(morlet_mat)
      phases <- Arg(morlet_mat)
      noise_mat <- noise_mat_r[[r]]
      phases_random <- Arg(noise_mat)
      mat_new <- matrix(complex(modulus = modulus, argument = phases_random),
                        ncol = ncol(phases_random))
      rec_orig = fun_icwt(x = morlet_mat) + mean(data[[l]]$Qobs)
      rec <- fun_icwt(x = mat_new)
      rec_random <- rec + mean(data[[l]]$Qobs)
      data_new <- data.frame(random = rec_random)
      data_new$MM <- data[[l]]$MM
      data_new$DD <- data[[l]]$DD
      data_new$YYYY <- data[[l]]$YYYY
      data_new$index <- data[[l]]$index
      data_new$seasonal <- data_new$random
      data_new$rank <- rank(data_new$seasonal)
      d <- 1
      data_new$simulated_seasonal <- NA
      for (d in c(1:365)) {
        data_day <- data[[l]][which(data[[l]]$index %in%
                                      c(d)), ]
        if (marginal_list[[l]] == "kappa") {
          colnames(par_day_list[[l]]) <- names(kap_par)
          data_day$kappa <- rand.kappaLD(length(data_day$Qobs),
                                         xi = par_day_list[[l]][d, "xi"], alfa = par_day_list[[l]][d,
                                                                                                   "alfa"], k = par_day_list[[l]][d, "k"],
                                         h = par_day_list[[l]][d, "h"])
          data_day$rank <- rank(data_day$kappa)
          data_new$rank <- rank(data_new$seasonal)
          data_new$rank[which(data[[l]]$index %in% c(d))] <- rank(data_new[which(data[[l]]$index %in%
                                                                                   c(d)), ]$seasonal)
          data_ordered <- data_day[order(data_day$rank),
                                   ]
          data_new$simulated_seasonal[which(data_new$index %in%
                                              c(d))] <- data_ordered$kappa[data_new$rank[which(data[[l]]$index %in%
                                                                                                 c(d))]]
          if (length(which(data_new$simulated_seasonal <
                           0)) > 0) {
            rep_value <- runif(n = 1, min = 0, max = min(data_day$Qobs))
            data_new$simulated_seasonal[which(data_new$simulated_seasonal <
                                                0)] <- rep_value
          }
        }
        if (marginal_list[[l]] == "empirical") {
          data_day$rank <- rank(data_day$Qobs)
          data_new$rank <- rank(data_new$seasonal)
          data_new$rank[which(data[[l]]$index %in% c(d))] <- rank(data_new[which(data[[l]]$index %in%
                                                                                   c(d)), ]$seasonal)
          data_ordered <- data_day[order(data_day$rank),
                                   ]
          data_new$simulated_seasonal[which(data_new$index %in%
                                              c(d))] <- data_ordered$Qobs[data_new$rank[which(data[[l]]$index %in%
                                                                                                c(d))]]
        }
        if (marginal_list[[l]] != "kappa" & marginal_list[[l]] !=
            "empirical") {
          data_day$cdf <- rCDF(n = length(data_day$Qobs),
                               par_day_list[[l]][d, ])
          data_day$rank <- rank(data_day$cdf)
          data_new$rank <- rank(data_new$seasonal)
          data_new$rank[which(data[[l]]$index %in% c(d))] <- rank(data_new[which(data[[l]]$index %in%
                                                                                   c(d)), ]$seasonal)
          data_ordered <- data_day[order(data_day$rank),
                                   ]
          data_new$simulated_seasonal[which(data_new$index %in%
                                              c(d))] <- data_ordered$cdf[data_new$rank[which(data[[l]]$index %in%
                                                                                               c(d))]]
        }
      }
      data_sim[[r]] <- data_new$simulated_seasonal
      if (verbose)
        cat(".")
    }
    if (verbose)
      cat("\nFinished.\n")
    data_sim <- as.data.frame(data_sim)
    names(data_sim) <- paste("r", seq(1:number_sim), sep = "")
    data_stoch <- data.frame(data[[l]][, c("YYYY", "MM",
                                           "DD", "timestamp", "Qobs")], data_sim)
    if (GoFtest == "NULL") {
      p_vals <- NULL
    }
    if (marginal != "empirical") {
      if (marginalpar) {
        out_list[[l]] <- list(simulation = data_stoch,
                              pars = par_day, p_val = p_vals)
      }
      else {
        out_list[[l]] <- list(simulation = data_stoch,
                              pars = NULL, p_val = p_vals)
      }
    }
    else {
      out_list[[l]] <- list(simulation = data_stoch, pars = NULL,
                            p_val = NULL)
    }
  }
  return(out_list)
}
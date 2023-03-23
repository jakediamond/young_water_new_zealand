# Data loading functions --------------------------------------------------
load_Rdata <- function(file) {
  # Results are saved in .Rdata file so need to load like this
  data <- load(file)
  return(get(data))
}

# Data cleaning/combining functions ---------------------------------------
# Clean raw data
clean_raw_fun <- function (rawdata) {
  rawdata %>%
    select(site = Site, date = Date.y,
           discharge = `Flow..m3.s.`, d18O, d2H) %>%
    mutate(date = lubridate::dmy(date),
           discharge = as.numeric(discharge)) %>%
    arrange(site, date) %>%
    mutate(date_simple = lubridate::round_date(date, unit = "month")) %>%
    filter(!(site == "NN1" & date == lubridate::ymd(20191015))) %>%
    drop_na() %>%
    distinct() %>%
    pivot_longer(cols = d18O, names_to = "variable", values_to = "value")
}

# Clean up data provided by Andy
clean_fun <- function (rivdata) { #, pptdata
  names(rivdata) <- tolower(names(rivdata))
  # names(pptdata) <- tolower(names(pptdata))
  # pdf <- mutate(pptdata,
  #               date = lubridate::make_date(year = year, month = month),
  #               type = "precip",
  #               variable = case_match(variable,
  #                                     dD ~ d2H,
  #                                     .default = variable)) %>%
  #   select(-site.number, -sampleid, -agent.number) %>%
  #   filter(variable %in% c("d18O", "d2H"))
  rdf <- select(rivdata, -`...1`) %>%
    rename(value = val)
  rdf
}

# Clean the already calculated sine model data from Andy
clean_fun_results <- function (data) {
  data %>%
    rename(site = Site,
           variable = `fitted quantity`,
           amp_riv_old = estimate,
           damp_riv_old = std.error,
           amp_ppt_unwt = Unweight.Avg.Precip.Amplitude,
           amp_ppt = Weight.Avg.Precip.Amplitude,
           amp_ppt_var = Weight.Avg.Precip.Variance,
           Fwy_old = YWF,
           Fwy_2.5_old = YWF.low,
           Fwy_97.5_old = YWF.hgh
           ) %>%
    select(-term) %>%
    mutate(damp_ppt = sqrt(amp_ppt_var))
}

# Functions to fit sine curve to isotope data -----------------------------
# Use multiple linear regression with iteratively reweighted least squares,
# which is the default for glm()
# d18O = A*cos(2*pi*f*t) + B*sin(2*pi*f*t)
sine_fit <- function (data) {
  # values to fit
  y <- data$value
  # Weights initialized to NULL
  wts <- NULL
  # Weights based on precipitation or discharge amount
  if ("discharge" %in% colnames(data)) {
    wts <- data$discharge
  }
  if ("precip" %in% colnames(data)) {
    wts <- data$precip
  }
  # Get date as a fraction of the year
  x <- lubridate::yday(data$date) / 366
  # cosine component
  xc <- cos(2*pi*x)
  # sine component
  xs <- sin(2*pi*x)
  # Fit the model
  fit <- glm(y ~ xc + xs, weights = wts)
}

# # Get the overall amplitude from the two components
# amp_fun <- function (fit) {
#   amp <- sqrt(coef(fit)[2]^2 + coef(fit)[3]^2)
# }
#
# # Get the overall phase
# phase_fun <- function (fit) {
#   phase <- -atan2(coef(fit)[2], coef(fit)[3])
# }
#
# # Get the overall offset
# offset_fun <- function (fit) {
#   offset <- coef(fit)[1]
# }

# Functions to map over all data ------------------------------------------
riv_analysis <- function(data) {
  # do the models
  df <- data %>%
    group_by(site, variable) %>%
    nest() %>%
    mutate(mods = map(data, sine_fit),
           tid = map(mods, broom::tidy)) %>%
    select(-data, -mods) %>%
    unnest(tid)

  # get in wider format for easier analysis
  df_w <- df %>%
    mutate(term = case_match(term, "(Intercept)" ~ "offset", .default = term)) %>%
    select(-p.value, -statistic) %>%
    rename(d = std.error) %>%
    pivot_wider(names_from = term, values_from = c(estimate, d),
                names_sep = "") %>%
    rename_with(str_remove, starts_with("estimate"), "estimate")

  # Get overall amplitude, phase, and offset
  df_res <- df_w %>%
    mutate_with_error(amp_riv ~ sqrt(xc^2 + xs^2)) %>%
    mutate_with_error(phase_riv ~ -atan(xc / xs))
}

# Function to get amplitude ratios
amp_ratio_fun <- function (rivdata, pptdata) {
  df_amp <- left_join(rivdata, pptdata) %>%
    mutate_with_error(Fwy_mean_Gauss ~ amp_riv / amp_ppt)
}

# Monte carlo functions
monte_carlo_fun <- function (data) {
  df <- data %>%
    filter(variable != "d2H") %>%
    group_by(site, variable) %>%
    mutate(
      # distributions of river amplitudes d18O amplitudes
      riv_dist = map2(amp_riv, damp_riv, ~rnorm(10000, .x, .y)),
      # distributions of precipitation d180 amplitudes
      ppt_dist = map2(amp_ppt, damp_ppt,
                    ~rnorm(10000, .x, .y)),
      # young water fraction distributions
      ywf_dist = map2(riv_dist, ppt_dist, ~.x / .y),
      # metrics based on distribution
      Fwy_mean_MC = map_dbl(ywf_dist, mean),
      Fwy_se_MC   = map_dbl(ywf_dist, sd) / sqrt(10000),
      Fwy_sd_MC   = map_dbl(ywf_dist, sd),
      Fwy_2.5_MC  = map_dbl(ywf_dist, quantile, 0.025),
      Fwy_med_MC  = map_dbl(ywf_dist, quantile, 0.5),
      Fwy_97.5_MC = map_dbl(ywf_dist, quantile, 0.975)) %>%
    select(-riv_dist, -ppt_dist, -ywf_dist)
}

# Functions to estimate gamma distribution parameters ---------------------
# gamma shape factor, alpha, using equation 11 from Kirchner 2016
# psi_riv - psi_ppt = alpha * arctan(sqrt((amp_riv / amp_ppt)^(-2 / alpha) -1))
shape_fun <- function (x, amp_ratio, phase_dif) {
  x * atan(sqrt(amp_ratio^(-2 / x) - 1)) - phase_dif
}

# function to iteratively solve for alpha
alpha_fun <- function (fun, amp_ratio, phase_dif, lims = c(0, 2)) {
  alpha <- uniroot(fun, lims, tol = 0.0001,
                   amp_ratio = amp_ratio,
                   phase_dif = phase_dif)$root
}

# gamma scale factor, beta, using equation 10 from Kirchner 2016
# beta = (1 / (2 * pi * f)) * sqrt((amp_riv / amp_ppt)^-(2 / alpha) - 1)
# assume yearly cycle, so f = 1 yr^-1
scale_fun <- function (alpha, amp_ratio, f = 1) {
  beta <- (1 / (2 * pi * f)) * sqrt((amp_ratio)^-(2 / alpha) - 1)
}

# Estimate the threshold age for young water, tau_yw with equation (14)
tau_fun <- function (alpha) {
  0.0949 + 0.1065 * alpha - 0.0126 * alpha^2
}

# Young water fraction function -------------------------------------------
# Estimate the young water fraction as Fyw = Gamma(tau, alpha, beta)
# Where Gamma is the lower incomplete gamma function with equation 13
Fyw_fun <- function (amp_ratio, phase_dif) {
  alpha <- alpha_fun(shape_fun, amp_ratio, phase_dif)
  beta <- scale_fun(alpha, amp_ratio)
  tau <- tau_fun(alpha)
  Fyw <- pgamma(tau, shape = alpha, scale = beta)
}

# Function to return weighted standard error of weighted mean ------------
# Based on Kirchner and Allen 2019
se_wt <- function (x, weights = NULL, normwt = FALSE,
                   na.rm = TRUE, method = c("unbiased", "ML"))
{
  method <- match.arg(method)
  if (!length(weights)) {
    if (na.rm)
      x <- x[!is.na(x)]
    return(var(x))
  }
  if (na.rm) {
    s <- !is.na(x + weights)
    x <- x[s]
    weights <- weights[s]
  }
  if (normwt)
    weights <- weights * length(x)/sum(weights)
  if (normwt || method == "ML")
    return(as.numeric(stats::cov.wt(cbind(x), weights, method = method)$cov))
  sw <- sum(weights)
  neff <- sw^2 / sum(weights^2)
  if (sw <= 1)
    warning("only one effective observation; variance estimate undefined")
  xbar <- sum(weights * x) / sw
  varw <- (sum(weights * ((x - xbar)^2)) / sw) * neff / (neff - 1)
  sqrt(varw / neff)
}
#
# df <- load_amp_data(here("src",
#                    "R",
#                    "sinusoidal model for river sites",
#                    "output",
#                    "river.site.est.Amplitude.Rdata"))
#

#
#
# x = df2$ywf_mean - df2$YWF
# y = x^2
# z = mean(sqrt(y))
#
#
#
#
#
#


#
#
# A function for propagating error with mutate function
# via Lee Pang (https://www.r-bloggers.com/2015/01/easy-error-propagation-in-r/)
mutate_with_error = function(.data, f) {
  exprs = list(
    # expression to compute new variable values
    deparse(f[[3]]),
    # expression to compute new variable errors
    sapply(all.vars(f[[3]]), function(v) {
      dfdp = deparse(D(f[[3]], v))
      sprintf('(d%s*(%s))^2', v, dfdp)
    }) %>%
      paste(collapse='+') %>%
      sprintf('sqrt(%s)', .)
  )
  names(exprs) = c(
    deparse(f[[2]]),
    sprintf('d%s', deparse(f[[2]]))
  )
  .data %>%
    # the standard evaluation alternative of mutate()
    mutate_(.dots=exprs)
}
#
#
# fx <- function(x1, x2, delta_x) -atan2(x1, x2)
# x <- c(1.27,0.806)
# delta_x <- c(0.470,0.499)
# fx(x[1], x, delta_x)
# gaussian_error_propagation <- function(x, fx, delta_x) {
#   # x: vector of input values
#   # fx: vector of output values corresponding to x
#   # delta_x: vector of input uncertainties
#
#   # Calculate the partial derivatives of fx with respect to x
#   grad_fx <- matrix(0, nrow = length(x), ncol = length(x))
#   for (i in 1:length(x)) {
#     grad_fx[,i] <- diff(fx(x[i], x, delta_x)) / diff(x)
#   }
#
#   # Calculate the covariance matrix of x
#   cov_x <- diag(delta_x^2)
#
#   # Calculate the covariance matrix of fx
#   cov_fx <- grad_fx %*% cov_x %*% t(grad_fx)
#
#   # Calculate the standard deviation of fx
#   sd_fx <- sqrt(diag(cov_fx))
#
#   return(sd_fx)
# }
# gaussian_error_propagation(c(1.27,0.806), fx, c(0.470,0.499))
#
# gaussian_error_propagation(c(1.27,0.806), fx, c(0.470,0.499))
#
# test <- tibble(x1 = 1.27, dx1 = 0.470, x2 = 0.806, dx2 = 0.499) %>%
#   mutate_with_error(y ~ sqrt(x1^2 + x2^2)) %>%
#   mutate_with_error(y2 ~ -atan(x1/x2))
# test
#
#
# gaussian_error_propagation2 <- function(f, x, dx) {
#   # Calculate the output of the function at the input values
#   y <- f(x)
#
#   num_deriv <- function(f, x, delta_x) {
#     f1 <- f(x + delta_x)
#     f2 <- f(x - delta_x)
#     return((f1 - f2) / (2 * delta_x))
#   }
#
#   # Calculate the partial derivatives of the function with respect to each input
#   df_dx <- sapply(seq_along(x), function(i) {
#     num_deriv(f, x[i], dx[i])
#   })
#
#   # Calculate the covariance matrix of the inputs
#   cov_mat <- diag(dx^2)
#
#   # Calculate the variance of the output using the partial derivatives and covariance matrix
#   var_y <- t(df_dx) %*% cov_mat %*% df_dx
#
#   # Return the output and its associated uncertainty as a list
#   return(list(y = y, dy = sqrt(var_y)))
# }
#
# # Example usage:
# f <- function(x) x^2 + 2*x + 1
# x <- c(1, 2, 3)
# dx <- c(0.1, 0.2, 0.1)
# result <- gaussian_error_propagation(f, x, dx)
# print(result)
# numericDeriv
#
# ?Hmisc::wtd.mean()
# Hmisc::wtd.var()

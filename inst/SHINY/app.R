#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

# Shiny app: Demonstrating MLE for multiple univariate distributions
# Save this file as app.R and run with: shiny::runApp('.') or open in RStudio and click Run App

library(shiny)
library(numDeriv) # for hessian when needed

# Negative log-likelihood functions --------------------------------------------------
negll_normal <- function(par, x){
  mu <- par[1]
  sigma <- exp(par[2])
  -sum(dnorm(x, mean = mu, sd = sigma, log = TRUE))
}

negll_exponential <- function(par, x){
  rate <- exp(par[1])
  -sum(dexp(x, rate = rate, log = TRUE))
}

negll_poisson <- function(par, x){
  lambda <- exp(par[1])
  -sum(dpois(x, lambda = lambda, log = TRUE))
}

negll_binomial <- function(par, x, size){
  p <- 1/(1+exp(-par[1])) # logistic transform to keep 0-1
  -sum(dbinom(x, size = size, prob = p, log = TRUE))
}

negll_gamma <- function(par, x){
  shape <- exp(par[1])
  rate <- exp(par[2])
  -sum(dgamma(x, shape = shape, rate = rate, log = TRUE))
}

negll_beta <- function(par, x){
  a <- exp(par[1])
  b <- exp(par[2])
  -sum(dbeta(x, shape1 = a, shape2 = b, log = TRUE))
}

# Generic MLE wrapper that returns estimates, se from inverse Hessian (when available)
fit_mle <- function(dist, x, extra = list()){
  res <- list()
  if(dist == "normal"){
    init <- c(mean(x), log(sd(x)))
    opt <- optim(init, negll_normal, x = x, hessian = TRUE, control = list(maxit=2000))
    mu_hat <- opt$par[1]
    sigma_hat <- exp(opt$par[2])
    se <- tryCatch({
      vcov <- solve(opt$hessian)
      c(mu = sqrt(vcov[1,1]), sigma = sqrt(vcov[2,2]) * sigma_hat) # delta for sigma
    }, error = function(e) c(mu = NA, sigma = NA))
    res$estimate <- c(mu = mu_hat, sigma = sigma_hat)
    res$se <- se
    res$negll <- opt$value
  } else if(dist == "exponential"){
    init <- log(1/mean(x))
    opt <- optim(init, negll_exponential, x = x, hessian = TRUE)
    rate_hat <- exp(opt$par[1])
    se <- tryCatch({
      vcov <- 1/opt$hessian
      se_rate <- sqrt(vcov) * rate_hat
      c(rate = se_rate)
    }, error = function(e) c(rate = NA))
    res$estimate <- c(rate = rate_hat)
    res$se <- se
    res$negll <- opt$value
  } else if(dist == "poisson"){
    init <- log(mean(x) + 1e-6)
    opt <- optim(init, negll_poisson, x = x, hessian = TRUE)
    lambda_hat <- exp(opt$par[1])
    se <- tryCatch({
      vcov <- 1/opt$hessian
      c(lambda = sqrt(vcov) * lambda_hat)
    }, error = function(e) c(lambda = NA))
    res$estimate <- c(lambda = lambda_hat)
    res$se <- se
    res$negll <- opt$value
  } else if(dist == "binomial"){
    size <- extra$size
    init <- 0
    opt <- optim(init, negll_binomial, x = x, size = size, hessian = TRUE)
    p_hat <- 1/(1+exp(-opt$par[1]))
    se <- tryCatch({
      vcov <- 1/opt$hessian
      se_logit <- sqrt(vcov)
      se_p <- se_logit * p_hat * (1 - p_hat) # delta method
      c(p = se_p)
    }, error = function(e) c(p = NA))
    res$estimate <- c(p = p_hat)
    res$se <- se
    res$negll <- opt$value
  } else if(dist == "gamma"){
    init <- c(log( (mean(x)^2)/var(x) ), log(mean(x)/var(x)))
    opt <- optim(init, negll_gamma, x = x, hessian = TRUE, control = list(maxit=2000))
    shape_hat <- exp(opt$par[1])
    rate_hat <- exp(opt$par[2])
    se <- tryCatch({
      vcov <- solve(opt$hessian)
      # approximate s.e. on original scale using delta method
      se_shape <- sqrt(vcov[1,1]) * shape_hat
      se_rate <- sqrt(vcov[2,2]) * rate_hat
      c(shape = se_shape, rate = se_rate)
    }, error = function(e) c(shape = NA, rate = NA))
    res$estimate <- c(shape = shape_hat, rate = rate_hat)
    res$se <- se
    res$negll <- opt$value
  } else if(dist == "beta"){
    init <- c(log(1), log(1))
    opt <- optim(init, negll_beta, x = x, hessian = TRUE, control = list(maxit=2000))
    a_hat <- exp(opt$par[1]); b_hat <- exp(opt$par[2])
    se <- tryCatch({
      vcov <- solve(opt$hessian)
      c(a = sqrt(vcov[1,1]) * a_hat, b = sqrt(vcov[2,2]) * b_hat)
    }, error = function(e) c(a = NA, b = NA))
    res$estimate <- c(a = a_hat, b = b_hat)
    res$se <- se
    res$negll <- opt$value
  }
  res
}

# UI -------------------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("MLE Demonstrator — Univariate Distributions"),
  sidebarLayout(
    sidebarPanel(
      selectInput("dist", "Choose distribution:",
                  choices = c("Normal" = "normal",
                              "Exponential" = "exponential",
                              "Poisson" = "poisson",
                              "Binomial" = "binomial",
                              "Gamma" = "gamma",
                              "Beta" = "beta")),
      numericInput("n", "Sample size (n):", value = 100, min = 1, step = 1),
      conditionalPanel("input.dist == 'normal'",
                       numericInput("true_mu", "True mu:", value = 0),
                       numericInput("true_sigma", "True sigma:", value = 1, min = 1e-6)),
      conditionalPanel("input.dist == 'exponential'",
                       numericInput("true_rate", "True rate:", value = 1, min = 1e-6)),
      conditionalPanel("input.dist == 'poisson'",
                       numericInput("true_lambda", "True lambda:", value = 3, min = 0)),
      conditionalPanel("input.dist == 'binomial'",
                       numericInput("true_p", "True p:", value = 0.4, min = 0, max = 1, step = 0.01),
                       numericInput("bin_size","Number of trials per obs (size):", value = 10, min = 1, step = 1)),
      conditionalPanel("input.dist == 'gamma'",
                       numericInput("true_shape", "True shape:", value = 2, min = 1e-6),
                       numericInput("true_rate_g", "True rate:", value = 1, min = 1e-6)),
      conditionalPanel("input.dist == 'beta'",
                       numericInput("true_a", "True alpha:", value = 2, min = 1e-6),
                       numericInput("true_b", "True beta:", value = 2, min = 1e-6)),
      actionButton("simulate","Simulate / Fit MLE"),
      hr(),
      numericInput("seed","Random seed (optional):", value = NA, min = 1),
      helpText("You can also paste your own sample values (comma-separated) into the box below."),
      textAreaInput("custom_data", "Custom sample (leave blank to use simulated sample)", rows = 3),
      width = 3
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Summary",
                 verbatimTextOutput("mle_summary"),
                 br(),
                 verbatimTextOutput("mle_se")),
        tabPanel("Plot",
                 plotOutput("distPlot", height = "500px")),
        tabPanel("Log-likelihood",
                 plotOutput("llPlot", height = "450px")),
        tabPanel("Data",
                 tableOutput("data_table"))
      ),
      width = 9
    )
  )
)

# Server ---------------------------------------------------------------------------
server <- function(input, output, session){
  # reactive sample
  sample_data <- eventReactive(input$simulate, {
    if(!is.na(input$seed)) set.seed(input$seed)
    # if user provided custom data
    if(nzchar(trimws(input$custom_data))){
      vals <- as.numeric(unlist(strsplit(input$custom_data, ",")))
      vals <- vals[!is.na(vals)]
      if(length(vals) < 1) showNotification("No valid numeric values in custom data.", type = "error")
      return(vals)
    }
    n <- input$n
    if(input$dist == "normal"){
      rnorm(n, mean = input$true_mu, sd = input$true_sigma)
    } else if(input$dist == "exponential"){
      rexp(n, rate = input$true_rate)
    } else if(input$dist == "poisson"){
      rpois(n, lambda = input$true_lambda)
    } else if(input$dist == "binomial"){
      rbinom(n, size = input$bin_size, prob = input$true_p)
    } else if(input$dist == "gamma"){
      rgamma(n, shape = input$true_shape, rate = input$true_rate_g)
    } else if(input$dist == "beta"){
      rbeta(n, shape1 = input$true_a, shape2 = input$true_b)
    }
  }, ignoreNULL = FALSE)
  
  fit_res <- reactive({
    x <- sample_data()
    if(is.null(x) || length(x) == 0) return(NULL)
    dist <- input$dist
    extra <- list()
    if(dist == "binomial") extra$size <- input$bin_size
    fit_mle(dist, x, extra = extra)
  })
  
  output$data_table <- renderTable({
    x <- sample_data()
    if(is.null(x)) return()
    head(data.frame(x = x), 50)
  })
  
  output$mle_summary <- renderPrint({
    res <- fit_res()
    if(is.null(res)) { cat("No results yet — simulate or provide data and click Simulate / Fit MLE") ; return() }
    cat("MLE estimates:\n")
    print(res$estimate)
    cat('\nNegative log-likelihood:', res$negll, '\n')
    cat('\nNotes: values are MLEs computed by numeric optimization (transforms used to respect parameter domains).')
  })
  
  output$mle_se <- renderPrint({
    res <- fit_res()
    if(is.null(res)) return()
    cat("Approx. standard errors (delta method where needed):\n")
    print(res$se)
  })
  
  output$distPlot <- renderPlot({
    x <- sample_data()
    if(is.null(x)) return()
    res <- fit_res()
    par(mfrow = c(1,1))
    dist <- input$dist
    if(dist %in% c("normal","exponential","gamma","beta")){
      hist(x, prob = TRUE, main = paste("Sample and fitted density —", dist), xlab = "x")
      xx <- seq(min(x), max(x), length.out = 400)
      if(dist == "normal") lines(xx, dnorm(xx, mean = res$estimate["mu"], sd = res$estimate["sigma"]), lwd = 2)
      if(dist == "exponential") lines(xx, dexp(xx, rate = res$estimate["rate"]), lwd = 2)
      if(dist == "gamma") lines(xx, dgamma(xx, shape = res$estimate["shape"], rate = res$estimate["rate"]), lwd = 2)
      if(dist == "beta") lines(xx, dbeta(xx, shape1 = res$estimate["a"], shape2 = res$estimate["b"]), lwd = 2)
    } else if(dist == "poisson"){
      tbl <- table(x)
      barplot(tbl/sum(tbl), main = "Empirical pmf vs fitted Poisson", xlab = "k")
      k <- as.numeric(names(tbl))
      points(k, dpois(k, lambda = res$estimate["lambda"]), pch = 19, col = "red")
      legend("topright", legend = c("empirical","fitted Poisson"), pch = c(15,19), col = c("grey","red"))
    } else if(dist == "binomial"){
      tbl <- table(x)
      barplot(tbl/sum(tbl), main = paste("Empirical pmf vs fitted Binomial (size=", input$bin_size, ")"), xlab = "k")
      k <- as.numeric(names(tbl))
      points(k, dbinom(k, size = input$bin_size, prob = res$estimate["p"]), pch = 19, col = "red")
      legend("topright", legend = c("empirical","fitted Binomial"), pch = c(15,19), col = c("grey","red"))
    }
  })
  
  output$llPlot <- renderPlot({
    x <- sample_data()
    if(is.null(x)) return()
    dist <- input$dist
    res <- fit_res()
    par(mfrow = c(1,1))
    # show log-likelihood across one parameter (fix others at MLE) when possible
    if(dist == "normal"){
      mu_seq <- seq(res$estimate["mu"] - 2*sd(x), res$estimate["mu"] + 2*sd(x), length.out = 200)
      ll <- sapply(mu_seq, function(m) -negll_normal(c(m, log(res$estimate["sigma"])), x))
      plot(mu_seq, ll, type = 'l', main = "Profile log-likelihood for mu", xlab = "mu", ylab = "log-likelihood")
      abline(v = res$estimate["mu"], col = 'red')
    } else if(dist == "exponential"){
      rate_seq <- seq(max(1e-6, res$estimate["rate"]/3), res$estimate["rate"]*3, length.out = 200)
      ll <- sapply(rate_seq, function(r) -negll_exponential(log(r), x))
      plot(rate_seq, ll, type = 'l', main = "Log-likelihood for rate", xlab = "rate", ylab = "log-likelihood")
      abline(v = res$estimate["rate"], col = 'red')
    } else if(dist == "poisson"){
      lam_seq <- seq(max(1e-6, res$estimate["lambda"]/3), res$estimate["lambda"]*3, length.out = 200)
      ll <- sapply(lam_seq, function(l) -negll_poisson(log(l), x))
      plot(lam_seq, ll, type = 'l', main = "Log-likelihood for lambda", xlab = "lambda", ylab = "log-likelihood")
      abline(v = res$estimate["lambda"], col = 'red')
    } else if(dist == "binomial"){
      p_seq <- seq(1e-4, 1-1e-4, length.out = 200)
      size <- input$bin_size
      ll <- sapply(p_seq, function(p) sum(dbinom(x, size = size, prob = p, log = TRUE)))
      plot(p_seq, ll, type = 'l', main = "Log-likelihood for p (Binomial)", xlab = "p", ylab = "log-likelihood")
      abline(v = res$estimate["p"], col = 'red')
    } else if(dist == "gamma"){
      shape_seq <- seq(max(1e-3, res$estimate["shape"]/3), res$estimate["shape"]*3, length.out = 200)
      ll <- sapply(shape_seq, function(s) -negll_gamma(c(log(s), log(res$estimate["rate"])), x))
      plot(shape_seq, ll, type = 'l', main = "Profile log-likelihood for shape", xlab = "shape", ylab = "log-likelihood")
      abline(v = res$estimate["shape"], col = 'red')
    } else if(dist == "beta"){
      a_seq <- seq(max(1e-3, res$estimate["a"]/3), res$estimate["a"]*3, length.out = 200)
      ll <- sapply(a_seq, function(a) -negll_beta(c(log(a), log(res$estimate["b"])), x))
      plot(a_seq, ll, type = 'l', main = "Profile log-likelihood for alpha", xlab = "alpha", ylab = "log-likelihood")
      abline(v = res$estimate["a"], col = 'red')
    }
  })
}

shinyApp(ui, server)

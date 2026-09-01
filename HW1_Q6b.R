set.seed(123)

p_true <- 23 * pi / 192
p_true

# sample size
Nmax <- 1e6

# Generate X,Y,Z
x <- runif(Nmax)
y <- runif(Nmax)
z <- runif(Nmax)

# Indicator for satisfying  inequalities
hit <- (x^2 + y^2 < z) & (z^2 > x * y)

# Estimate, sum(TRUE)/total
running_est <- cumsum(hit) / seq_len(Nmax)

# Sample sizes at which to plot the estimate
N_plot <- unique(round(10^seq(1, 6, length.out = 200)))

plot(
  N_plot,
  running_est[N_plot],
  type = "l",
  log = "x",
  xlab = "Sample size N",
  ylab = "Monte Carlo estimate",
  main = "Monte Carlo Estimate of Probability"
)


abline(h = p_true, lty = 2, lwd = 2)

legend(
  "topright",
  legend = c("Monte Carlo estimate", "Analytic value"),
  lty = c(1, 2),
  lwd = c(1, 2),
  bty = "n"
)

running_est[Nmax]
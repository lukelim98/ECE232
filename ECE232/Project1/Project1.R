library(igraph)
# 1a
# Probabilities
p_values <- c(0.002, 0.006, 0.012, 0.045, 0.1)
# Number of nodes
n <- 900
# Create graphs, calculate and plot degree distributions, and report metrics
for (p in p_values) {
  g <- erdos.renyi.game(n, p, type='gnp', directed=FALSE)
  deg <- degree(g)
  
  # Plotting the degree distribution
  hist(deg, main=paste('Degree Distribution p=', p), xlab='Degree', ylab='Frequency', breaks=50)
  
  # Calculating mean and variance
  mean_deg <- mean(deg)
  var_deg <- var(deg)
  
  cat("For p =", p, "\n")
  cat("Mean of degree distribution:", mean_deg, "\n")
  cat("Variance of degree distribution:", var_deg, "\n")
  
  # Theoretical mean and variance
  theoretical_mean <- p * (n-1)
  theoretical_var <- p * (1-p) * (n-1)
  
  cat("Theoretical mean:", theoretical_mean, "\n")
  cat("Theoretical variance:", theoretical_var, "\n\n")
}



# 1b
# Parameters
p_values <- c(0.002, 0.006, 0.012, 0.045, 0.1)
n <- 900
num_realizations <- 100 # Number of realizations for each p to estimate connectivity probability

# Loop through each p value
for (p in p_values) {
  connected_count <- 0
  
  for (i in 1:num_realizations) {
    g <- erdos.renyi.game(n, p, type = "gnp", directed = FALSE)
     
    # Check if the network is connected
    if (is_connected(g)) {
      connected_count <- connected_count + 1
    } else if (i == 1) { # For the first non-connected instance, find GCC and its diameter
      # Identifying the GCC
      cl <- clusters(g)
      gcc <- induced.subgraph(g, which(cl$membership == which.max(cl$csize)))
      
      # Calculating the diameter of the GCC
      gcc_diameter <- diameter(gcc)
      cat("For p =", p, ", the diameter of the GCC in the first non-connected instance is:", gcc_diameter, "\n")
    }
  }
  
  # Estimate the probability of connectivity
  connectivity_prob <- connected_count / num_realizations
  cat("For p =", p, ", the estimated probability of a network being connected is:", connectivity_prob, "\n\n")
}



# 1c
n <- 900
p_values <- seq(0, log(n)/n, length.out = 100) # Incrementing p up to pmax
num_realizations <- 100 # Number of realizations for each p

# Initialize a list to store the normalized sizes of the GCC for each p
gcc_sizes <- vector("list", length(p_values))
names(gcc_sizes) <- as.character(p_values)

# Sweep over values of p
for (p in p_values) {
  sizes <- numeric(num_realizations) # Store sizes for each realization
  
  for (i in 1:num_realizations) {
    g <- erdos.renyi.game(n, p, type = "gnp", directed = FALSE)
    cl <- clusters(g)
    # Size of the largest component divided by total number of nodes
    sizes[i] <- max(cl$csize) / n
  }
  
  gcc_sizes[[as.character(p)]] <- sizes
}

# Calculate average normalized GCC sizes for each p
average_sizes <- sapply(gcc_sizes, mean)

# Plotting
plot(p_values, average_sizes, type = "l", col = "red", xlab = "p", ylab = "Normalized GCC Size", main = "Normalized GCC Size vs p")
points(rep(p_values, each = num_realizations), unlist(gcc_sizes), col = "blue", cex = 0.5)

# Empirical estimation of emergence and dominance of GCC
# Criterion for emergence: We define "emergence" as the first p where the average size jumps significantly
emergence_p <- p_values[which(diff(average_sizes) > 0.05)[1]] # Example criterion
# Criterion for dominance: First p where average size exceeds 0.99
dominance_p <- p_values[which(average_sizes > 0.99)[1]]

cat("Emergence of GCC starts at p =", emergence_p, "\n")
cat("GCC dominates over 99% of the nodes at p =", dominance_p, "\n")


# 1d.i
n_values <- seq(100, 10000, by=100)
c <- 0.5
gcc_sizes <- numeric(length(n_values))

for (i in 1:length(n_values)) {
  n <- n_values[1]
  p <- c/n
  g <- erdos.renyi.game(n, p, type='gnp')
  cl <- clusters(g)
  gcc_sizes[i] <- max(cl$csize) / n # Normalized GCC size
}

plot(n_values, gcc_sizes, type = "l", xlab = "Number of Nodes (n)", ylab = "Normalized GCC Size", main = "GCC Size vs n for c=0.5")



# 1d.ii
n_values <- seq(100, 10000, by=100)
c <- 1
gcc_sizes <- numeric(length(n_values))

for (i in 1:length(n_values)) {
  n <- n_values[1]
  p <- c/n
  g <- erdos.renyi.game(n, p, type='gnp')
  cl <- clusters(g)
  gcc_sizes[i] <- max(cl$csize) / n # Normalized GCC size
}

plot(n_values, gcc_sizes, type = "l", xlab = "Number of Nodes (n)", ylab = "Normalized GCC Size", main = "GCC Size vs n for c=1")



# 1d.iii
c_values <- c(1.15, 1.25, 1.35)
colors <- c("red", "green", "blue")
plot(NULL, xlim = c(min(n_values), max(n_values)), ylim = c(0, 1), xlab = "Number of Nodes (n)", ylab = "Normalized GCC Size", main = "GCC Size vs n for Various c")

for (j in 1:length(c_values)) {
  c <- c_values[j]
  gcc_sizes <- numeric(length(n_values))
  
  for (i in 1:length(n_values)) {
    n <- n_values[i]
    p <- c / n
    g <- erdos.renyi.game(n, p, type = "gnp")
    cl <- clusters(g)
    gcc_sizes[i] <- max(cl$csize) / n # Normalized GCC size
  }
  
  lines(n_values, gcc_sizes, col = colors[j], type = "l", lwd = 2)
}

legend("topright", legend = paste("c =", c_values), col = colors, lwd = 2)







# 2a
n <- 1050 # number of nodes
m <- 1 # Each new node attaches to m old nodes

# Creating the network
g <- barabasi.game(n, m=m, directed=FALSE)

# Check if the network is connected
is_connected <- is_connected(g)

cat("Is the network alwasy connected?", is_connected, "\n")



# 2b
# Community detection with fast greedy method
communities <- cluster_fast_greedy(g)

# Measure modularity
mod <- modularity(communities)
cat("Modularty of the network:", mod, "\n")

# Compute the degree assortativity 
assortativity  <- assortativity_degree(g)
cat("Degree Assortativity of the network:", assortativity , "\n")

################################################################################
######################                SETUP              #######################   
################################################################################

# Set working directory, packages and data -- only on macOS
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# Load libraries
library(ggplot2)
library(patchwork)
library(nleqslv)
library(data.table)

# Set seed for random observations
set.seed(234259203)

################################################################################
######################            DEFINE GRID            #######################   
################################################################################

# The first step consists in creating a grid with a set of locations that we will
# use as units of observations

# Choose number of points to make the grid
n <- 20

# Number of locations 
N <- n * n

# Create coordinates of points in the grid
x_coords <- seq(1, n, by = 1)
y_coords <- seq(1, n, by = 1)

grid_data <- expand.grid(x = x_coords, y = y_coords)

# We will assume that there is an homogeneous constant speed across our city
v <- 30         #Average commuting speed (km/h)
c <- 0.25       #Average commuting cost (euro/km)

# Initialize distance matrix and travel time matrix
distance_matrix <- matrix(0, nrow = N, ncol = N)
travel_time <- matrix(NA, nrow = N, ncol = N)
travel_cost <- matrix(NA, nrow = N, ncol = N)

# Calculate pairwise distances and travel times
for (i in 1:N) {
  for (j in 1:N){
    distance_matrix[i,j] <- abs(grid_data$x[i] - grid_data$x[j])+abs(grid_data$y[i] - grid_data$y[j])     #Barcelona distance
    travel_time[i,j] <- distance_matrix[i,j] / v
    travel_cost[i,j] <- distance_matrix[i,j] * c
  }
}


for (i in 1:N) {
  travel_time[i,i] <- min(travel_time[i,-i]) 
  travel_cost[i,i] <- min(travel_cost[i,-i]) 
}

#Find central location (min distance to all other locations)
loc = which.min(rowSums(distance_matrix))


################################################################################
###################       RANDOM GRID CITY OBSERVATIONS      ###################   
################################################################################

# Draw a random matrix for commuting trips
R_ij <- matrix(NA, nrow = N, ncol = N)

for (i in 1:N) {
  for (j in 1:N) {
    R_ij[i,j] <- rlnorm(1, meanlog = 8 - 0.5 * ((grid_data$x[j]-grid_data$x[loc])^2 + (grid_data$y[j]-grid_data$y[loc])^2)^(1/2), sdlog = 0.5)  
  }
}

# Calculate population in every cell
R_i <- apply(R_ij,1,sum)
grid_data$population <- R_i

p1 <- ggplot(grid_data, aes(x = x, y = y, fill = population)) +
              geom_tile() +
              scale_fill_gradient(low = "white", high = "blue") +
              theme_minimal() +
              labs(title = "Population distribution", fill = "Population",  x = NULL, y = NULL)

# Calculate employment in every cell
R_j <- apply(R_ij,2,sum)
grid_data$employment <- R_j

p2 <- ggplot(grid_data, aes(x = x, y = y, fill = employment)) +
              geom_tile() +
              scale_fill_gradient(low = "white", high = "red") +
              theme_minimal() +
              labs(title = "Employment distribution", fill = "Employment",  x = NULL, y = NULL)

# Draw a random wage vector
w_j <- matrix(NA,  nrow = N, ncol = 1)

for (i in 1:N) {
  w_j[i] <- rlnorm(1, meanlog = 3 - 0.1 * (abs(grid_data$x[i]-grid_data$x[loc]) + abs(grid_data$y[i]-grid_data$y[loc])), sdlog = 0.1)
}

grid_data$wages <- w_j

p3 <- ggplot(grid_data, aes(x = x, y = y, fill = wages)) +
             geom_tile() +
             scale_fill_gradient(low = "white", high = "green") +
             theme_minimal() +
             labs(title = "Wage distribution", fill = "Wage",  x = NULL, y = NULL)


# Draw a random rents vector
q_i <- matrix(NA,  nrow = N, ncol = 1)

for (i in 1:N) {
  q_i[i] <- rlnorm(1, meanlog = 1 - 0.1 * (abs(grid_data$x[i]-grid_data$x[loc]) + abs(grid_data$y[i]-grid_data$y[loc])), sdlog = 0.1)
}

grid_data$rents <- q_i

p4 <- ggplot(grid_data, aes(x = x, y = y, fill = rents)) +
             geom_tile() +
             scale_fill_gradient(low = "white", high = "purple") +
             theme_minimal() +
             labs(title = "Rent distribution", fill = "Rents",  x = NULL, y = NULL)

wrap_plots(p1, p2, p3, p4)



# Draw a random array for non-commuting trips
      grid_data <- as.data.table(grid_data)
      # Step 1: Precompute the mean for each k
      loc_half_x <- grid_data$x[loc / 2]
      loc_y <- grid_data$y[loc]
      
      grid_data[, mean_k := 1 - 0.5 * sqrt((x - loc_half_x)^2 + (y - loc_y)^2)]
      
      # Step 2: Create the long table of all (i, j, k) combinations
      ijk <- as.data.table(expand.grid(i = 1:N, j = 1:N, k = 1:N))
      
      # Step 3: Join mean_k values by k
      ijk[, mean := grid_data$mean_k[k]]
      
      # Step 4: Generate random values and apply truncation at 0
      set.seed(234259203)  # optional: for reproducibility
      ijk[, n_ijk := pmax(rnorm(.N, mean = mean, sd = 1), 0)]
      
      # Step 5: Convert back to 3D array
      n_ijk <- array(ijk$n_ijk, dim = c(N, N, N))

#Aggregate them by destination
n_k <- apply(n_ijk, 3, sum)

grid_data$non_commuting_trips <- n_k

ggplot(grid_data, aes(x = x, y = y, fill = non_commuting_trips)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  labs(title = "Non-commuting trips distribution", fill = "Non-commuting trips",  x = NULL, y = NULL)

################################################################################
######################       DEFINE CITY PARAMETERS      #######################   
################################################################################

#delta <- 3        #Non-commuting trips generalized cost elasticility
#nu <- 8
#eta <- 3

alpha <- 0.3      # Diminishing returns time per visit
beta <- 0.3       # Share of time spent working
gamma <- 0.75     # Share of income spent in consumption
rho <- 0.7        # Elasticity of substitution among locations
theta <- 5        # Commuting elasticity parameter
mu <- 0.7         # Share of revenue spent in labour


T <- 24 * 7       # Total time endowment (h)
W <- 8            # Daily working time (h)


#Composite parameters. These are the parameters that we can actually estimate
delta <- (1-alpha*rho)/(1-rho)
nu <- ((1-alpha) * (1-beta) + beta) * theta
eta <- (1-rho)/rho * (1-beta) * theta



################################################################################
######################        DERIVED OBSERSATIONS       #######################   
################################################################################

# Step 1: Compute vot_ij
ij <- as.data.table(expand.grid(i = 1:N, j = 1:N))
ij[, w_j := w_j[j]]
ij[, travel_cost := as.vector(travel_cost)]
ij[, travel_time := as.vector(travel_time)]
ij[, vot_ij := pmax(w_j * W - travel_cost, 0) / (W + travel_time)]

# Step 2: Compute tau_ijk = travel_cost[i,k] + vot_ij[i,j] * travel_time[i,k]
ijk <- as.data.table(expand.grid(i = 1:N, j = 1:N, k = 1:N))
ijk[, travel_cost_ik := as.vector(travel_cost)[(i - 1) * N + k]]
ijk[, travel_time_ik := as.vector(travel_time)[(i - 1) * N + k]]
ijk[, vot_ij := ij[.SD, on = .(i, j), vot_ij]]
ijk[, tau_ijk := travel_cost_ik + vot_ij * travel_time_ik]

# Step 3: Compute totalTravelTime[i,j,k] = n_ijk * travel_time[i,k]
ijk[, n_ijk := as.vector(n_ijk)]
ijk[, total_time := n_ijk * travel_time_ik]
ijk[, total_cost := n_ijk * travel_cost_ik]

# Step 4: Aggregate over k
ij_totals <- ijk[, .(
  totalTravelTime = sum(total_time, na.rm = TRUE),
  totalTravelCost = sum(total_cost, na.rm = TRUE)
), by = .(i, j)]

# Step 5: Compute x_ij
ij_totals[, travel_time_ij := ij[.SD, on = .(i, j), travel_time]]
ij_totals[, travel_cost_ij := ij[.SD, on = .(i, j), travel_cost]]
ij_totals[, w_j := ij[.SD, on = .(i, j), w_j]]

ij_totals[, x_ij := 
            beta / (alpha * (1 - beta) + beta) * (T - totalTravelTime) / (W + travel_time_ij) +
            alpha * (1 - beta) / (alpha * (1 - beta) + beta) * totalTravelCost / (w_j * W - travel_cost_ij)
]

# Step 6: Compute N_ij = x_ij * R_ij * W
ij_totals[, R_ij := as.vector(R_ij)]
ij_totals[, N_ij := x_ij * R_ij * W]

# Step 7: Aggregate N_j = sum over i
N_j <- ij_totals[, .(N_j = sum(N_ij, na.rm = TRUE)), by = j][order(j), N_j]

# Step 8: Compute A_j
A_j <- w_j / (mu * N_j^(mu - 1))

tau_ijk <- array(ijk$tau_ijk, dim = c(N, N, N))


################################################################################
##########################         FUNCTIONS         ###########################
################################################################################

inversionAmenities <- function(a) {

      # Initial value in the positive simplex
      # a <- a / sum(a)
                      
      # Step 1: Reshape n_ijk and tau_ijk to long format
      ijk <- as.data.table(expand.grid(i = 1:N, j = 1:N, k = 1:N))
      ijk[, n_ijk := as.vector(n_ijk)]
      ijk[, tau_ijk := as.vector(tau_ijk)]
      
      # Step 2: Compute n_o = n_ijk * tau_ijk
      ijk[, n_o := n_ijk * tau_ijk]
      
      # Step 3: Compute N_o[k]
      N_o_dt <- ijk[, .(N_o = sum(n_o, na.rm = TRUE)), by = k]
      
      # Step 4: Compute n_ij[i,j]
      n_ij_dt <- ijk[, .(n_ij = sum(n_o, na.rm = TRUE)), by = .(i, j)]
      
      # Step 5: Precompute a² and tau^(−(delta−1))
      a_sq <- a^2
      tau_pow <- ijk[, tau_pow := tau_ijk^(-(delta - 1))]
      
      # Step 6: Compute weighted share for each (i,j,k)
      ijk[, numerator := a_sq[k] * tau_pow]
      ijk[, denom := sum(a_sq * tau_pow), by = .(i, j)]
      ijk[, n_d := numerator / denom * n_ij_dt[.SD, on = .(i, j), n_ij]]
      
      # Step 7: Aggregate N_d[k]
      N_d_dt <- ijk[, .(N_d = sum(n_d, na.rm = TRUE)), by = k]
      
      # Step 8: Final vector
      N_dt <- merge(N_o_dt, N_d_dt, by = "k")
      N <- N_dt[, N_o - N_d]
                      
return (N)}

inversionWorkplace <- function(lambda_j) {
      
  # Create long table of (i, j) pairs
      dt <- as.data.table(expand.grid(i = 1:N, j = 1:N))
      dt[, lambda_j := lambda_j[j]]
      dt[, A_ij := as.vector(A_ij)]
      dt[, vot_ij := as.vector(vot_ij)]
      
      # Calculate unnormalised probabilities
      dt[, pi_j_i := lambda_j * A_ij^eta * vot_ij^nu]
      
      # Normalise by row i
      dt[, denom := sum(pi_j_i), by = i]
      dt[, pi_j_i := pi_j_i / denom]
      
      # Merge R_i for weighted sum
      dt[, R_i := R_i[i]]
      
      # Aggregate to get predicted R_j
      R_j_p <- dt[, .(R_j_p = sum(pi_j_i * R_i)), by = j][order(j), R_j_p]
      
      # Return residual
      return(R_j - R_j_p)
      
}

inversionResidence <- function(lambda_i) {
      # Create long table of (i, j) pairs
      dt <- as.data.table(expand.grid(i = 1:N, j = 1:N))
      dt[, lambda_i := lambda_i[i]]
      dt[, q_i := q_i[i]]
      dt[, A_ij := as.vector(A_ij)]
      dt[, vot_ij := as.vector(vot_ij)]
      
      # Calculate unnormalised probabilities
      dt[, pi_i_j := lambda_i * q_i * A_ij^eta * vot_ij^nu]
      
      # Normalise by column j
      dt[, denom := sum(pi_i_j), by = j]
      dt[, pi_i_j := pi_i_j / denom]
      
      # Merge R_j for weighted sum
      dt[, R_j := R_j[j]]
      
      # Aggregate to get predicted R_i
      R_i_p <- dt[, .(R_i_p = sum(pi_i_j * R_j)), by = i][order(i), R_i_p]
      
      # Return residual
      return(R_i - R_i_p)
}


################################################################################
##########################         INVERSION         ###########################
################################################################################

#D raw a random amenities vector
a_init <- matrix(1/N,  nrow = N, ncol = 1)

solution <- nleqslv(a_init, inversionAmenities, method = "Broyden", control = list(allowSingular = TRUE))

grid_data$amenities <- (solution$x)^2 / sum((solution$x)^2)

p5 <- ggplot(grid_data, aes(x = x, y = y, fill = amenities)) +
             geom_tile() +
             scale_fill_gradient(low = "white", high = "red") +
             theme_minimal() +
             labs(title = "Amenities Distribution", fill = "Amenities",  x = NULL, y = NULL)


A_ij <- matrix(NA,  nrow = N, ncol = N)

# Step 1: Flatten tau_ijk into a long table
ijk <- as.data.table(expand.grid(i = 1:N, j = 1:N, k = 1:N))
ijk[, tau := as.vector(tau_ijk)]

# Step 2: Precompute weights
a2 <- (solution$x)^2
ijk[, weight := a2[k]]

# Step 3: Compute tau^(-(delta - 1))
ijk[, tau_power := tau^(-(delta - 1))]

# Step 4: Compute A_ij = sum_k weight[k] * tau_power[i,j,k]
A_dt <- ijk[, .(A_val = sum(weight * tau_power, na.rm = TRUE)), by = .(i, j)]

# Step 5: Convert to wide matrix format
A_ij <- dcast(A_dt, i ~ j, value.var = "A_val")[, -1] |> as.matrix()

lambda_j_init <- matrix(1/N,  nrow = N, ncol = 1)

solution <- nleqslv(lambda_j_init, inversionWorkplace, method = "Broyden", control = list(allowSingular = TRUE))

grid_data$workplace_fe <- solution$x / sum(solution$x)

p6 <- ggplot(grid_data, aes(x = x, y = y, fill = workplace_fe)) +
             geom_tile() +
             scale_fill_gradient(low = "white", high = "green") +
             theme_minimal() +
             labs(title = "Workplace FE Distribution", fill = "Workplace FE",  x = NULL, y = NULL)


lambda_i_init <- matrix(1/N,  nrow = N, ncol = 1)

solution <- nleqslv(lambda_i_init, inversionResidence, method = "Broyden", control = list(allowSingular = TRUE))

grid_data$residence_fe <- solution$x / sum((solution$x))

p7 <- ggplot(grid_data, aes(x = x, y = y, fill = residence_fe)) +
             geom_tile() +
             scale_fill_gradient(low = "white", high = "blue") +
             theme_minimal() +
             labs(title = "Residence Fixed Effect Distribution", fill = "Residence FE",  x = NULL, y = NULL)

grid_data$productivity <- A_j

p8 <- ggplot(grid_data, aes(x = x, y = y, fill = productivity)) +
             geom_tile() +
             scale_fill_gradient(low = "white", high = "purple") +
             theme_minimal() +
             labs(title = "Productivity Distribution", fill = "Productivity",  x = NULL, y = NULL)

wrap_plots(p5, p6, p7, p8)

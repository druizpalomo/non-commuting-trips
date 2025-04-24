################################################################################
######################                SETUP              #######################   
################################################################################

#Load libraries
library(ggplot2)
library(patchwork)
library(nleqslv)

#Set seed for random observations
set.seed(234259203)

################################################################################
######################            DEFINE GRID            #######################   
################################################################################

#The first step consists in creating a grid with a set of locations that we will
#use as units of observations

#Choose number of points to make the grid
n <- 5

#Number of locations 
N <- n * n

#Create coordinates of points in the grid
x_coords <- seq(1, n, by = 1)
y_coords <- seq(1, n, by = 1)

grid_data <- expand.grid(x = x_coords, y = y_coords)

#We will assume that there is an homogeneous constant speed across our city
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

#Draw a random matrix for commuting trips
R_ij <- matrix(NA, nrow = N, ncol = N)

for (i in 1:N) {
  for (j in 1:N) {
    R_ij[i,j] <- rlnorm(1, meanlog = 8 - 0.5 * ((grid_data$x[j]-grid_data$x[loc])^2 + (grid_data$y[j]-grid_data$y[loc])^2)^(1/2), sdlog = 0.5)  
  }
}

#Calculate population in every cell
R_i <- apply(R_ij,1,sum)
grid_data$population <- R_i

p1 <- ggplot(grid_data, aes(x = x, y = y, fill = population)) +
              geom_tile() +
              scale_fill_gradient(low = "white", high = "blue") +
              theme_minimal() +
              labs(title = "Population distribution", fill = "Population",  x = NULL, y = NULL)

#Calculate employment in every cell
R_j <- apply(R_ij,2,sum)
grid_data$employment <- R_j

p2 <- ggplot(grid_data, aes(x = x, y = y, fill = employment)) +
              geom_tile() +
              scale_fill_gradient(low = "white", high = "red") +
              theme_minimal() +
              labs(title = "Employment distribution", fill = "Employment",  x = NULL, y = NULL)

#Draw a random wage vector
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


#Draw a random rents vector
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

#Draw a random array for non-commuting trips
n_ijk <- array(NA, dim = c(N, N, N))

for (i in 1:N) {
  for (j in 1:N) {
    for (k in 1:N) {
      n_ijk[i,j,k] <- max(rnorm(1, mean = 1 - 0.5 * ((grid_data$x[k]-grid_data$x[loc/2])^2 + (grid_data$y[k]-grid_data$y[loc])^2)^(1/2), sd = 1),0)
    }
  }
}

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

alpha <- 0.3      #Diminishing returns time per visit
beta <- 0.3       #Share of time spent working
gamma <- 0.75     #Share of income spent in consumption
rho <- 0.7        #Elasticity of substitution among locations
theta <- 5        #Commuting elasticity parameter
mu <- 0.7         #Share of revenue spent in labour


T <- 24 * 7       #Total time endowment (h)
W <- 8            #Daily working time (h)


#Composite parameters. These are the parameters that we can actually estimate
delta <- (1-alpha*rho)/(1-rho)
nu <- ((1-alpha) * (1-beta) + beta) * theta
eta <- (1-rho)/rho * (1-beta) * theta



################################################################################
######################        DERIVED OBSERSATIONS       #######################   
################################################################################

vot_ij <- matrix(NA, nrow = N, ncol = N)
vot_ij <- pmax(matrix(rep(w_j, each = N), nrow = N) * W - travel_cost, 0) / (W + travel_time)

tau_ijk <- array(NA, dim = c(N, N, N))

for (i in 1:N) {          
  for (j in 1:N) {        
    for (k in 1:N) {      
      tau_ijk[i,j,k] <- travel_cost[i,k] + vot_ij[i,j] * travel_time[i,k]  
    }
  }
}

totalTravelTime <- array(NA, dim = c(N, N, N))
totalTravelCost <- array(NA, dim = c(N, N, N))

for (i in 1:N) {
  for (j in 1:N) {
    for (k in 1:N) {
      totalTravelTime[i,j,k] <- n_ijk[i,j,k] * travel_time[i,k]
      totalTravelCost[i,j,k] <- n_ijk[i,j,k] * travel_cost[i,k]
    }
  }
}

totalTravelTime <- apply(totalTravelTime, c(1,2), sum)
totalTravelCost <- apply(totalTravelCost, c(1,2), sum)

x_ij <- matrix(NA, nrow = N, ncol = N)

for (i in 1:N) {
  for (j in 1:N) {
    x_ij[i,j] <- beta/(alpha*(1-beta)+beta) * (T - totalTravelTime[i,j]) / (W+travel_time[i,j]) + 
      alpha*(1-beta)/(alpha * (1-beta) + beta) * totalTravelCost[i,j] / (w_j[j]*W-travel_cost[i,j])
  }
}

N_ij <- matrix(NA, nrow = N, ncol = N) 

for (i in 1:N) {
  for (j in 1:N) {
    N_ij[i,j] <- x_ij[i,j] * R_ij[i,j] * W
  }
}

N_j <- apply(N_ij, 2, sum)

A_j <- w_j /(mu * N_j^(mu-1))

################################################################################
##########################         FUNCTIONS         ###########################
################################################################################

inversionAmenities <- function(a) {

          
                      #Initial value in the positive simplex
                      #a <- a / sum(a)
                      
                      
                      n_o <- array(NA, dim = c(N, N, N))
                      for (i in 1:N) {          
                        for (j in 1:N) {     
                          for (k in 1:N) {
                            n_o[i,j,k] = n_ijk[i,j,k] * tau_ijk[i,j,k] 
                          }
                        }
                      }
                      
                      N_o <- apply(n_o,3, sum)
                      
                      n_ij <- apply(n_o, c(1,2), sum)
                      
                      n_d <- array(NA, dim = c(N, N, N))
                      
                      for (i in 1:N) {          
                        for (j in 1:N) {        
                          for (k in 1:N) {      
                            n_d[i,j,k] <- (a[k]^2 * tau_ijk[i,j,k]^(-(delta-1))) / (tau_ijk[i,j,]^(-(delta-1)) %*% a^2) * n_ij[i,j]
                          }
                        }
                      }
  
                      N_d <- apply(n_d,3, sum)
                      
                      N <- N_o - N_d
                      
                      return (N)}

inversionWorkplace <- function(lambda_j){
                      
                      pi_j_i <- matrix(NA,  nrow = N, ncol = N)
                      for (i in 1:N) {
                        for (j in 1:N) {
                          pi_j_i[i,j] <- lambda_j[j] * A_ij[i,j]^eta * vot_ij[i,j]^nu  
                        }
                      }
                      
                      pi_i <- apply(pi_j_i, 1, sum)
                      
                      for (i in 1:N) {
                        for (j in 1:N) {
                          pi_j_i[i,j] <- pi_j_i[i,j] / pi_i[i]  
                        }
                      }
                      
                      R_j_p <- matrix(NA,  nrow = N, ncol = 1)
                      
                      for (j in 1:N) { 
                        R_j_p[j] <- pi_j_i[,j] %*% R_i 
                      }
                      
                      R_j_o <- R_j
                      
                      R_j <- R_j_o - R_j_p
                      
                      return(R_j)}

inversionResidence <- function(lambda_i){
                      
                      pi_i_j <- matrix(NA, nrow = N, ncol = N)
                      
                      for (i in 1:N) {
                        for (j in 1:N) {
                          pi_i_j[i,j] <- lambda_i[i] * q_i[i] * A_ij[i,j]^eta * vot_ij[i,j]^nu
                        }
                      }
                      
                      pi_j <- apply(pi_i_j, 2, sum)
                      
                      for (i in 1:N) {
                        for (j in 1:N) {
                          pi_i_j[i,j] <- pi_i_j[i,j] / pi_j[j]
                        }
                      }
                      
                      R_i_p <- matrix(NA,  nrow = N, ncol = 1)
                      
                      for (i in 1:N) { 
                        R_i_p[i] <- pi_i_j[i,] %*% R_j 
                      }
                      
                      R_i_o <- R_i
                      
                      R_i <- R_i_o - R_i_p
                      
                      return(R_i)}


################################################################################
##########################         INVERSION         ###########################
################################################################################

#Draw a random amenities vector
a_init <- matrix(1/N,  nrow = N, ncol = 1)

solution <- nleqslv(a_init, inversionAmenities, method = "Broyden", control = list(allowSingular = TRUE))

grid_data$amenities <- (solution$x)^2 / sum((solution$x)^2)

p5 <- ggplot(grid_data, aes(x = x, y = y, fill = amenities)) +
             geom_tile() +
             scale_fill_gradient(low = "white", high = "red") +
             theme_minimal() +
             labs(title = "Amenities Distribution", fill = "Amenities",  x = NULL, y = NULL)


A_ij <- matrix(NA,  nrow = N, ncol = N)

for (i in 1:N) {
  for (j in 1:N) {
    A_ij[i,j] <- (solution$x)^2 %*% tau_ijk[i,j,]^(-(delta-1))
  }
}

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
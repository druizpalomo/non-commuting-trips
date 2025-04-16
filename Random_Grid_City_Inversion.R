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
n <- 7

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

#Draw a random wage vector
w_j <- matrix(NA,  nrow = N, ncol = 1)

for (i in 1:N) {
  w_j[i] <- rlnorm(1, meanlog = 3 - 0.1 * sqrt((grid_data$x[i]-grid_data$x[loc])^2 + (grid_data$y[i]-grid_data$y[loc])^2), sdlog = 0.1)
}

grid_data$wages <- w_j

ggplot(grid_data, aes(x = x, y = y, fill = wages)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  labs(title = "Wage Distribution", fill = "Wage",  x = NULL, y = NULL)

#Draw a random matrix for commuting trips
#R_ij <- matrix(NA, nrow = N, ncol = N)

#for (i in 1:N) {
#  for (j in 1:N) {
#    R_ij[i,j] <- rlnorm(1, meanlog = 12 - sqrt((grid_data$x[j]-grid_data$x[loc])^2 + (grid_data$y[j]-grid_data$y[loc])^2), sdlog = 0.5)  
#  }
#}

#Draw a random array for non-commuting trips
n_ijk <- array(NA, dim = c(N, N, N))

for (i in 1:N) {
  for (j in 1:N) {
    for (k in 1:N) {
      n_ijk[i,j,k] <- max(rnorm(1, mean = 1 - 0.5 * (abs(grid_data$x[k]-grid_data$x[loc]) + abs(grid_data$y[k]-grid_data$y[loc])), sd = 1),0)
    }
  }
}


################################################################################
######################       DEFINE CITY PARAMETERS      #######################   
################################################################################

#alpha <- 0.3      #Share of revenues spent in labour
beta <- 0.3       #Share of time spent working
gamma <- 0.75     #Share of income spent in consumption

rho <- 0.8        #Elasticity of substitution among locations
alpha <- 0.9      #Diminishing returns time per visit

theta <- 5        #Commuting elasticity parameter

T <- 24 * 7       #Total time endowment (h)
W <- 8            #Daily working time (h)


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

################################################################################
##########################         FUNCTIONS         ###########################
################################################################################

inversionAmenities <- function(a) {

          
                      #Initial value in the positive simplex
                      a <- a / sum(a)
                      
                      
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
                            n_d[i,j,k] <- (a[k] * tau_ijk[i,j,k]^(-rho*(1-alpha)/(1-rho))) / (tau_ijk[i,j,]^(-rho*(1-alpha)/(1-rho)) %*% a) * n_ij[i,j]
                          }
                        }
                      }
  
                      N_d <- apply(n_d,3, sum)
                      
                      N <- N_o - N_d
                      
                      return (N)}


################################################################################
##########################         INVERSION         ###########################
################################################################################

#Draw a random amenities vector
a_init <- matrix(1/N,  nrow = N, ncol = 1)

solution <- nleqslv(a_init, inversionAmenities, method = "Broyden", control = list(allowSingular = TRUE))

grid_data$amenities <- solution$x / sum(solution$x)

ggplot(grid_data, aes(x = x, y = y, fill = amenities)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  labs(title = "Amenities Distribution", fill = "Amenities",  x = NULL, y = NULL)

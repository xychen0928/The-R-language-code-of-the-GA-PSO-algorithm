#The goal is f6

popmin <- c(10,10,50,20.5)            # Define lower bounds for each variable (4-dimensional problem)
popmax <- c(30,30,70,22.5)            # Define upper bounds for each variable (4-dimensional problem)
lenchrom <- 4                         # Number of variables/chromosome length (problem dimension = 4)
pc <- 0.9                             # Initial crossover probability (can also be dynamic or adaptive)
# The initial crossover probability is set to 0.9 (later overwritten)
pm <- 0.1                             # Initial mutation probability (can also be dynamic or adaptive)
# The initial mutation probability is set to 0.1 (later overwritten)
#w1=0.7,w2=0.3                        # Commented example: possible weighting coefficients for hybrid algorithms
pc <- 0.8                             # Overwrite: final crossover probability = 0.8
pm <- 0.2                             # Overwrite: final mutation probability = 0.2
maxgen <- 200                         # Number of generations (evolution iterations)
popsize <- 100                        # Population size (number of individuals)

# Particle velocity limits
Vmax <- 1                             # Upper limit for velocity components
Vmin <- -1                            # Lower limit for velocity components

# Dynamic learning factor functions (remain unchanged)
c1fun <- function(t) 2*sin(pi/2*(1 - t/maxgen))^2
# Define time-varying learning factor c1(t)
# Formula: 2 * [ sin( (π/2) * (1 - t/maxgen) ) ]^2
# As t increases from 0 to maxgen, c1 changes according to this expression

c2fun <- function(t) 2*sin(pi*t/(2*maxgen))^2
# Define time-varying learning factor c2(t)
# Formula: 2 * [ sin( π * t / (2*maxgen) ) ]^2
# c2 changes complementarily to c1 to balance individual and global learning in PSO

# Inertia weight parameters
wmax <- 0.9                           # Maximum inertia weight
wmin <- 0.4                           # Minimum inertia weight

# Initialize population (4 columns)
par_num <- lenchrom                   # Store number of parameters for convenience
Z <- matrix(0, popsize, par_num)      # Initialize popsize x par_num matrix Z for chaotic sequence generation

set.seed(0)                           # Set random seed for reproducibility
Z[1, ] <- runif(par_num)              # First row: uniform random numbers in (0,1)
for(i in 2:popsize) {
  Z[i, ] <- 4.0 * Z[i-1, ] * (1 - Z[i-1, ])
}
# Use Logistic map to generate a chaotic sequence:
# For each subsequent row: Z[i,] = 4 * Z[i-1,] * (1 - Z[i-1,])
# This is a common method for producing pseudo-random samples with high diversity

pop <- matrix(0, popsize, par_num)    # Initialize population position matrix pop (each row = one individual)
V <- matrix(0, popsize, par_num)      # Initialize velocity matrix V (same dimensions as pop)
fitness <- numeric(popsize)           # Initialize fitness vector

set.seed(0)                           # Reset random seed to ensure reproducibility
for(i in 1:popsize) {
  pop[i, ] <- popmin + (popmax - popmin) * Z[i, ]
  # Map Z from [0,1] to [popmin, popmax] for each dimension
  V[i, ] <- (Vmax - Vmin)*runif(par_num, -1, 1) + Vmin
  # Initialize velocity uniformly in [Vmin, Vmax]
  fitness[i] <- f6(pop[i, ])          # Evaluate fitness using user-defined objective function f6 (4D version)
  # Note: f6 should accept a 1x4 vector and return a scalar (lower is better)
}

# Initialize best values
bestindex <- which.min(fitness)       # Index of individual with minimum fitness (best)
zbest <- pop[bestindex, ]             # Global best position
gbest <- pop                          # Personal best positions (initialize as current population)
fitnessgbest <- fitness               # Personal best fitness values
fitnesszbest <- fitness[bestindex]    # Global best fitness value

# Record fitness history
yy <- numeric(maxgen)                 # Store global best fitness per generation (for plotting or analysis)

set.seed(0)                           # Reset random seed again for reproducibility in main loop
# Main optimization loop (automatically handles 4D case)
for(ii in 1:maxgen) {
  c1 <- c1fun(ii)                     # Compute current learning factor c1
  c2 <- c2fun(ii)                     # Compute current learning factor c2
  
  # Compute average and best fitness
  favg <- mean(fitness)
  fmin <- min(fitness)
  
  for(j in 1:popsize) {
    # Dynamic inertia weight
    if(fitness[j] <= favg) {
      w <- wmin + (fitness[j] - fmin)*(wmax - wmin)/(favg - fmin)
    } else {
      w <- wmax
    }
    # Explanation:
    # - If fitness[j] <= favg, assign smaller inertia (closer to wmin)
    #   so better individuals exploit more locally.
    # - Otherwise assign wmax, encouraging worse individuals to explore more broadly.
    # Note: if favg == fmin, division by zero could occur (no explicit safeguard here)
    
    # Update velocity (4D automatically)
    V[j, ] <- w*V[j, ] + c1*runif(1)*(gbest[j, ] - pop[j, ]) + 
      c2*runif(1)*(zbest - pop[j, ])
    # Velocity update components:
    # - w * V[j, ]: inertia term
    # - c1 * r1 * (gbest[j,] - pop[j,]): cognitive term (individual learning)
    # - c2 * r2 * (zbest - pop[j,]): social term (global learning)
    # Here runif(1) produces a scalar, so the same random factor applies to all dimensions
    
    V[j, ] <- pmin(pmax(V[j, ], Vmin), Vmax)
    # Clamp velocity to [Vmin, Vmax] to prevent instability
    
    # Update position
    pop[j, ] <- pop[j, ] + 0.5*V[j, ]
    # Step size scaled by 0.5 to avoid excessive movement
    
    pop[j, ] <- pmin(pmax(pop[j, ], popmin), popmax)
    # Project positions back to [popmin, popmax]
    
    # Evaluate fitness
    fitness[j] <- f6(pop[j, ])
    
    # Update personal best
    if(fitness[j] < fitnessgbest[j]) {
      gbest[j, ] <- pop[j, ]
      fitnessgbest[j] <- fitness[j]
    }
    
    # Update global best
    if(fitness[j] < fitnesszbest) {
      zbest <- pop[j, ]
      fitnesszbest <- fitness[j]
    }
  }
  
  # Genetic algorithm operators (handle 4D automatically)
  # Crossover
  sortindexbest <- order(fitnessgbest)
  # Sort individuals by personal best fitness (ascending)
  
  numPool <- round(pc * popsize)
  # Number of individuals selected for crossover
  
  Poolx <- pop[sortindexbest[1:numPool], ]
  # Best numPool individuals (positions)
  
  PoolVx <- V[sortindexbest[1:numPool], ]
  # Corresponding velocities
  
  for(i in 1:numPool) {
    seed1 <- sample(numPool, 1)
    seed2 <- sample(numPool, 1)
    while(seed1 == seed2) seed2 <- sample(numPool, 1)
    # Randomly select two distinct parents
    
    pb <- runif(1)
    childx1 <- pb*Poolx[seed1, ] + (1-pb)*Poolx[seed2, ]
    # Arithmetic crossover: linear interpolation
    # child = pb * parent1 + (1-pb) * parent2
    
    childv1 <- (PoolVx[seed1, ] + PoolVx[seed2, ]) * 
      norm(as.matrix(PoolVx[seed1, ]), "2") / 
      norm(as.matrix(PoolVx[seed1, ] + PoolVx[seed2, ]), "2")
    # Child velocity:
    # - Combine parent velocities (vector sum)
    # - Scale to have same 2-norm as parent1’s velocity
    #   (heuristic to preserve momentum direction and magnitude)
    
    if(f6(pop[i, ]) > f6(childx1)) {
      pop[i, ] <- childx1
      V[i, ] <- childv1
    }
    # Replace individual with child if child has better (lower) fitness
  }
  
  # Mutation (4D automatically)
  mutationpool <- round(pm * popsize)
  # Number of individuals subject to mutation
  
  for(i in 1:mutationpool) {
    seed3 <- sample(popsize, 1)
    mutationchild <- pop[seed3, ] * (1 + rnorm(par_num))
    # Mutation strategy:
    # Multiply each dimension by (1 + Gaussian noise N(0,1))
    mutationchild <- pmin(pmax(mutationchild, popmin), popmax)
    # Clamp mutated values to [popmin, popmax]
    
    if(f6(pop[i, ]) > f6(mutationchild)) {
      pop[i, ] <- mutationchild
    }
    # If the mutated individual is better, replace current individual i
    # Note: uses seed3 for mutation, but replaces index i — may be intentional or a minor logic inconsistency
  }
  
  # Update fitness
  for(q in 1:popsize) {
    fitness[q] <- f6(pop[q, ])
    if(fitness[q] < fitnessgbest[q]) {
      gbest[q, ] <- pop[q, ]
      fitnessgbest[q] <- fitness[q]
    }
    if(fitness[q] < fitnesszbest) {
      zbest <- pop[q, ]
      fitnesszbest <- fitness[q]
    }
  }
  # Recalculate all fitness values and update personal/global bests for consistency
  
  yy[ii] <- fitnesszbest
  # Record the global best fitness at this generation
}

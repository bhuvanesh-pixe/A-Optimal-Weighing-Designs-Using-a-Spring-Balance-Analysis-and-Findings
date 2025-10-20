# Deriving the A-Optimal Weighing Design for n=11
# Based on principles in the paper by Pena Pardo and Sarkar

# The paper mentions that D-optimal designs for n=11 are connected to 
# the symmetric BIBD(11,5,2) derived from the Paley biplane
# Let's implement this using the principles described in the paper

# Function to generate a symmetric BIBD(11,5,2) based on Paley biplane principles
generate_bibd_11_5_2 <- function() {
  # According to the paper, the symmetric BIBD(11,5,2) can be constructed 
  # from the Paley biplane
  
  # We'll implement this using the cyclic construction method mentioned in the paper
  # where designs can be found by rotation about the center
  
  # The first row of the incidence matrix (one "base block")
  # Note: This follows the structure described in the paper for the BIBD(11,5,2)
  base_block <- c(1,1,1,1,1,0,0,0,0,0,0)
  
  # Generate the cyclic design
  design <- matrix(0, nrow=11, ncol=11)
  for(i in 1:11) {
    # Rotate the base block i-1 positions
    design[i,] <- c(base_block[((i-1):10 %% 11) + 1], base_block[1:(i-1)])
  }
  
  return(design)
}

# Alternative implementation based on quadratic residues
# as the paper mentions connections to Paley construction
generate_paley_bibd <- function() {
  # For a prime p = 11, quadratic residues are important in Paley constructions
  # Quadratic residues modulo 11: 1, 3, 4, 5, 9
  q_residues <- c(1, 3, 4, 5, 9)
  
  # Initialize the design matrix
  design <- matrix(0, nrow=11, ncol=11)
  
  # Construct the design based on differences
  for(i in 0:10) {
    row <- rep(0, 11)
    for(qr in q_residues) {
      # Calculate position using modular arithmetic
      pos <- (i + qr) %% 11 + 1
      row[pos] <- 1
    }
    design[i+1,] <- row
  }
  
  return(design)
}

# Function to implement the incidence matrix for n=11 as shown in the paper
generate_paper_design <- function() {
  # This reproduces the exact design from the paper
  X <- matrix(c(
    1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1,
    1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0,
    0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1,
    1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1,
    1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1,
    1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0,
    0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0,
    0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0,
    0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1,
    1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0,
    0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1
  ), nrow=11, byrow=TRUE)
  
  return(X)
}

# We'll also generate the complement of the BIBD
get_complement_design <- function(design) {
  return(1 - design)
}

# Calculate information matrix
calculate_info_matrix <- function(X) {
  t(X) %*% X
}

# Calculate the trace of the inverse of the information matrix (A-optimality criterion)
calculate_trace <- function(X) {
  info_matrix <- calculate_info_matrix(X)
  inv_info <- solve(info_matrix)
  return(sum(diag(inv_info)))
}

# Calculate the determinant of the information matrix (D-optimality criterion)
calculate_determinant <- function(X) {
  info_matrix <- calculate_info_matrix(X)
  return(det(info_matrix))
}

# Generate all possible designs by row and column permutations
# For illustration, we'll generate a few variations
generate_variations <- function(base_design, num_variations = 10) {
  designs <- list()
  designs[[1]] <- base_design
  
  n <- nrow(base_design)
  
  # Generate variations by random row permutations
  for(i in 2:num_variations) {
    row_perm <- sample(n)
    designs[[i]] <- base_design[row_perm,]
  }
  
  return(designs)
}

# Evaluate all designs and find the one with minimum trace
find_min_trace_design <- function(designs) {
  min_trace <- Inf
  min_idx <- 1
  traces <- numeric(length(designs))
  
  for(i in 1:length(designs)) {
    traces[i] <- calculate_trace(designs[[i]])
    if(traces[i] < min_trace) {
      min_trace <- traces[i]
      min_idx <- i
    }
  }
  
  cat("Trace values for all designs:\n")
  print(traces)
  
  return(list(design = designs[[min_idx]], trace = min_trace, index = min_idx))
}

# Main execution
cat("Generating designs based on principles in the paper...\n")

# Get the design directly from the paper
paper_design <- generate_paper_design()
cat("Design directly from the paper:\n")
print(paper_design)

# Calculate its trace (A-optimality)
paper_trace <- calculate_trace(paper_design)
cat("\nTrace of (X'X)^-1 for the paper design:", paper_trace, "\n")

# Calculate its determinant (D-optimality)
paper_det <- calculate_determinant(paper_design)
cat("Determinant of X'X for the paper design:", paper_det, "\n")

# Now let's try to generate BIBD-based designs as per principles
cat("\nGenerating BIBD-based design...\n")
bibd_design <- generate_paley_bibd()
cat("BIBD(11,5,2) design:\n")
print(bibd_design)

# Get complement as per the paper
cat("\nComplement design (BIBD(11,6,3)):\n")
complement_design <- get_complement_design(bibd_design)
print(complement_design)

# Calculate trace for the complement design
complement_trace <- calculate_trace(complement_design)
cat("\nTrace of (X'X)^-1 for the complement design:", complement_trace, "\n")

# Generate variations to find minimum trace
cat("\nGenerating variations to find minimum trace...\n")
variations <- generate_variations(paper_design, 20)
min_trace_result <- find_min_trace_design(variations)

cat("\nDesign with minimum trace (", min_trace_result$trace, "):\n")
print(min_trace_result$design)

# Compare information matrices
cat("\nInformation matrix (X'X) for minimum trace design:\n")
info_min <- calculate_info_matrix(min_trace_result$design)
print(info_min)

cat("\nInverse information matrix (X'X)^-1 for minimum trace design:\n")
inv_info_min <- solve(info_min)
print(inv_info_min)

# Display the trace components (diagonal elements of the inverse info matrix)
cat("\nDiagonal elements of (X'X)^-1 (contributing to the trace):\n")
print(diag(inv_info_min))

# Display actual weighing sets for the min-trace design
cat("\nWeighing sets for the minimum trace design:\n")
for(i in 1:nrow(min_trace_result$design)) {
  items <- which(min_trace_result$design[i,] == 1)
  cat("Set", i, ":", paste(items, collapse=", "), "\n")
}
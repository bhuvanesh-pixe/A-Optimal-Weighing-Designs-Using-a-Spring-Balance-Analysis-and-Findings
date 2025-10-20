# D-Optimal Weighing Design for n=8
# Based on the principles from the paper by Pena Pardo and Sarkar

# From Table 1 in the paper:
# For n=8, max det(X) = 56, and the number of D-optimal matrices is 195,955,200

# The paper mentions how designs for n=8 can sometimes be constructed 
# by extending designs for n=7. This is shown in Section 6.

# Function to generate D-optimal matrix for n=8 based on n=7 extension
generate_n8_from_n7 <- function() {
  # First, create the D-optimal matrix for n=7 (as shown in Section 5.1)
  X7 <- matrix(c(
    1, 1, 1, 1, 0, 0, 0,
    1, 1, 0, 0, 1, 1, 0,
    1, 0, 1, 0, 1, 0, 1,
    1, 0, 0, 1, 0, 1, 1,
    0, 1, 1, 0, 0, 1, 1,
    0, 1, 0, 1, 1, 0, 1,
    0, 0, 1, 1, 1, 1, 0
  ), nrow=7, byrow=TRUE)
  
  # Add a row and column as shown in Section 6 of the paper
  X8 <- matrix(0, nrow=8, ncol=8)
  X8[1:7, 1:7] <- X7
  X8[1:7, 8] <- 1  # Add a column of 1s
  X8[8, 1:7] <- 1  # Add a row of 1s
  X8[8, 8] <- 0    # The diagonal entry is 0
  
  return(X8)
}

# Alternative approach: Generate a D-optimal matrix for n=8 using the 
# construction principles inspired by Hadamard matrices
generate_n8_hadamard_based <- function() {
  # Create a base pattern that satisfies the properties
  # of a D-optimal weighing design for n=8
  X8 <- matrix(c(
    1, 1, 1, 1, 0, 0, 0, 1,
    1, 1, 0, 0, 1, 1, 0, 1,
    1, 0, 1, 0, 1, 0, 1, 1,
    1, 0, 0, 1, 0, 1, 1, 1,
    0, 1, 1, 0, 0, 1, 1, 1,
    0, 1, 0, 1, 1, 0, 1, 1,
    0, 0, 1, 1, 1, 1, 0, 1,
    1, 1, 1, 1, 1, 1, 1, 0
  ), nrow=8, byrow=TRUE)
  
  return(X8)
}

# Function to calculate the determinant of a matrix
calculate_determinant <- function(X) {
  return(det(X))
}

# Function to calculate the information matrix
calculate_info_matrix <- function(X) {
  return(t(X) %*% X)
}

# Function to calculate the trace of the inverse of the information matrix
calculate_trace <- function(X) {
  info_matrix <- calculate_info_matrix(X)
  inv_info <- solve(info_matrix)
  return(sum(diag(inv_info)))
}

# Function to verify the D-optimality of the design
verify_optimality <- function(X) {
  det_value <- calculate_determinant(X)
  # According to Table 1, the maximum determinant for n=8 is 56
  expected_max_det <- 56
  
  cat("Determinant of X:", det_value, "\n")
  if(abs(det_value) == expected_max_det) {
    cat("This design is D-optimal (maximum determinant achieved).\n")
  } else {
    cat("This design is NOT D-optimal (maximum determinant not achieved).\n")
  }
}

# Function to display the design properties
display_design_properties <- function(X) {
  n <- nrow(X)
  cat("\nWeighing Design Properties for n =", n, "\n")
  cat("--------------------------------\n")
  
  # Calculate and display the determinant
  det_value <- calculate_determinant(X)
  cat("Determinant of X:", det_value, "\n")
  
  # Calculate and display the information matrix
  info_matrix <- calculate_info_matrix(X)
  cat("\nInformation Matrix (X'X):\n")
  print(info_matrix)
  
  # Calculate and display the inverse of the information matrix
  inv_info_matrix <- solve(info_matrix)
  cat("\nInverse of Information Matrix (X'X)^-1:\n")
  print(inv_info_matrix)
  
  # Calculate and display the trace (A-optimality measure)
  trace_value <- sum(diag(inv_info_matrix))
  cat("\nTrace of (X'X)^-1 (A-optimality measure):", trace_value, "\n")
  
  # Generate the geometric representation described in section 6
  cat("\nGeometric Representation (1s are represented by 'X', 0s by '.'):\n")
  for(i in 1:n) {
    row_str <- ""
    for(j in 1:n) {
      if(X[i,j] == 1) {
        row_str <- paste0(row_str, "X ")
      } else {
        row_str <- paste0(row_str, ". ")
      }
    }
    cat(row_str, "\n")
  }
}

# Function to generate the weighing sets from the design matrix
generate_weighing_sets <- function(X) {
  n <- nrow(X)
  weighing_sets <- list()
  
  for(i in 1:n) {
    items <- which(X[i,] == 1)
    weighing_sets[[i]] <- items
  }
  
  return(weighing_sets)
}

# Function to display the weighing sets
display_weighing_sets <- function(weighing_sets) {
  cat("\nWeighing Sets (which items to weigh together):\n")
  for(i in 1:length(weighing_sets)) {
    cat("Set", i, ":", paste(weighing_sets[[i]], collapse=", "), "\n")
  }
}

# Main execution
cat("Generating D-optimal Weighing Design for n=8\n")
cat("Based on principles from the paper by Pena Pardo and Sarkar\n\n")

# Method 1: Extension from n=7
cat("Method 1: Extending n=7 design to n=8 as described in Section 6\n")
X8_method1 <- generate_n8_from_n7()
cat("Matrix X for n=8 (Method 1):\n")
print(X8_method1)
verify_optimality(X8_method1)
display_design_properties(X8_method1)
weighing_sets_method1 <- generate_weighing_sets(X8_method1)
display_weighing_sets(weighing_sets_method1)

# Method 2: Hadamard-based construction
cat("\n\nMethod 2: Construction based on Hadamard matrix principles\n")
X8_method2 <- generate_n8_hadamard_based()
cat("Matrix X for n=8 (Method 2):\n")
print(X8_method2)
verify_optimality(X8_method2)
display_design_properties(X8_method2)
weighing_sets_method2 <- generate_weighing_sets(X8_method2)
display_weighing_sets(weighing_sets_method2)

# Compare the trace (A-optimality) of both methods
trace_method1 <- calculate_trace(X8_method1)
trace_method2 <- calculate_trace(X8_method2)

cat("\n\nComparing A-optimality (trace of (X'X)^-1):\n")
cat("Method 1 trace:", trace_method1, "\n")
cat("Method 2 trace:", trace_method2, "\n")

if(trace_method1 < trace_method2) {
  cat("Method 1 is more A-optimal (has lower trace)\n")
} else if(trace_method2 < trace_method1) {
  cat("Method 2 is more A-optimal (has lower trace)\n")
} else {
  cat("Both methods have the same A-optimality\n")
}

# Generate variations to search for better A-optimality
generate_variations <- function(base_design, num_variations = 20) {
  designs <- list()
  designs[[1]] <- base_design
  
  n <- nrow(base_design)
  
  # Generate variations by random row and column permutations
  for(i in 2:num_variations) {
    row_perm <- sample(n)
    designs[[i]] <- base_design[row_perm, row_perm]
  }
  
  return(designs)
}

# Search for the most A-optimal design among variations
cat("\n\nSearching for the most A-optimal design among variations...\n")
variations <- generate_variations(X8_method1, 20)
traces <- sapply(variations, calculate_trace)
best_idx <- which.min(traces)

cat("Best trace found:", traces[best_idx], "\n")
cat("Matrix with best trace:\n")
print(variations[[best_idx]])

# Final analysis of the best design
cat("\n\nFinal Analysis of the Best Design:\n")
display_design_properties(variations[[best_idx]])
best_weighing_sets <- generate_weighing_sets(variations[[best_idx]])
display_weighing_sets(best_weighing_sets)
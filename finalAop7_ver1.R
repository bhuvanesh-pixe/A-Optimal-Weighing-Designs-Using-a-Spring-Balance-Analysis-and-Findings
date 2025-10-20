# Step 1: Define the 7 lines of the Fano plane using indices 1 to 7
fano_lines <- list(
  c(1, 2, 3),
  c(1, 4, 5),
  c(1, 6, 7),
  c(2, 4, 6),
  c(2, 5, 7),
  c(3, 4, 7),
  c(3, 5, 6)
)

# Step 2: Generate all 30 labelings (via all permutations of 7 labels, then de-duplicate)
library(gtools)

labels <- c("A", "B", "C", "D", "E", "F", "G")
perms <- permutations(n = 7, r = 7, v = labels)
perms <- perms[!duplicated(apply(perms, 1, paste, collapse = "")), ] # Remove dupes

# Only take 30 evenly spaced permutations to represent the 30 DWDs
set.seed(1)
selected_perms <- perms[sample(1:nrow(perms), 30), ]

# Step 3: Function to convert a labeling into a complemented design matrix
get_design_matrix <- function(labeling) {
  # Get the Fano plane blocks with actual labels
  blocks <- lapply(fano_lines, function(line) labeling[line])
  
  # Initialize incidence matrix (7x7)
  item_labels <- labeling
  incidence <- matrix(0, nrow = 7, ncol = 7)
  colnames(incidence) <- item_labels
  
  for (i in 1:7) {
    subset <- blocks[[i]]
    incidence[i, subset] <- 1
  }
  
  # Complement incidence matrix (i.e., flip 1s and 0s)
  comp_incidence <- 1 - incidence
  
  return(comp_incidence)
}

# Step 4: Function to compute trace of (X^T X)^-1
compute_trace <- function(X) {
  M <- t(X) %*% X
  if (det(M) == 0) return(Inf)  # Non-invertible matrix
  return(sum(diag(solve(M))))
}

# Step 5: Evaluate all 30 designs
results <- data.frame(index = 1:30, trace = NA)

for (i in 1:30) {
  labeling <- selected_perms[i, ]
  X <- get_design_matrix(labeling)
  results$trace[i] <- compute_trace(X)
}

# Step 6: Display the best A-optimal DWD(s)
best_index <- which.min(results$trace)
cat("Minimum trace:", results$trace[best_index], "\n")
cat("Best design index:", best_index, "\n")
cat("Best labeling:", selected_perms[best_index, ], "\n")

# Optional: View all trace values
print(results)
# Step 7: Show the design matrix with the minimum trace
best_labeling <- selected_perms[best_index, ]
best_matrix <- get_design_matrix(best_labeling)

cat("\nBest D-optimal design matrix (complemented Fano incidence):\n")
print(best_matrix)


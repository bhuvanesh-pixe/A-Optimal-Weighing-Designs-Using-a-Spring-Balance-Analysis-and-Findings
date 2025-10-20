library(Matrix)
library(gtools)
library(parallel)

# 1. Generate all relevant binary row vectors of length n (natural order)
generateRelevantSubsets <- function(n) {
  bits <- lapply(1:(2^n - 1), function(i) intToBits(i)[1:n])
  lapply(bits, function(x) as.integer(rev(x)))
}

# 2. Generate full-rank matrices without row permutations
generateFullRankDesigns <- function(n) {
  subs <- generateRelevantSubsets(n)
  combos <- combn(length(subs), n, simplify = FALSE)
  designMatrices <- list()
  counter <- 1
  
  for (idx in combos) {
    mat <- do.call(rbind, subs[idx])
    if (qr(mat)$rank == n) {
      sig <- paste0(c(mat), collapse = "")  # Unique signature for exact row order
      designMatrices[[counter]] <- list(matrix = mat, row_signature = sig)
      counter <- counter + 1
    }
  }
  designMatrices[1:(counter - 1)]
}

# 3. A-optimality: trace of inverse of XtX
compute_A_optimality <- function(X) {
  XtX <- tcrossprod(X)
  if (det(XtX) < 1e-8) return(Inf)
  sum(diag(solve(XtX)))
}

# 4. Process single matrix
process_matrix_base <- function(item) {
  trace_val <- compute_A_optimality(item$matrix)
  if (!is.finite(trace_val)) return(NULL)
  list(matrix = item$matrix, trace = trace_val, sig = item$row_signature)
}

# 5. Main driver
find_strict_base_A_optimal_designs <- function(n) {
  matrices <- generateFullRankDesigns(n)
  if (length(matrices) == 0) return(NULL)
  
  cl <- makeCluster(detectCores() - 1)
  clusterExport(cl, c("compute_A_optimality", "process_matrix_base"))
  clusterEvalQ(cl, library(Matrix))
  
  base_results <- parLapply(cl, matrices, process_matrix_base)
  stopCluster(cl)
  
  traces <- vapply(base_results, function(res) if (!is.null(res)) res$trace else Inf, numeric(1))
  minTrace <- min(traces)
  
  seen <- new.env(hash = TRUE, parent = emptyenv())
  unique_bases <- list()
  trace_match_count <- 0
  
  for (res in base_results) {
    if (!is.null(res) && abs(res$trace - minTrace) < 1e-8) {
      trace_match_count <- trace_match_count + 1
      if (!exists(res$sig, envir = seen)) {
        assign(res$sig, TRUE, envir = seen)
        unique_bases[[length(unique_bases) + 1]] <- res$matrix
      }
    }
  }
  
  list(
    "Min Trace" = minTrace,
    "Optimal Matrices" = unique_bases,
    "Count All Matrices with Same Trace" = trace_match_count
  )
}

# --- Run ---
n <- 5
cat("\nFinding strict base A-optimal designs for n =", n, "\n")
system.time({
  result <- find_strict_base_A_optimal_designs(n)
})

# --- Output ---
if (!is.null(result)) {
  cat("\nMinimum Trace:", result$`Min Trace`, "\n")
  cat("Number of Unique A-Optimal Base Matrices:", length(result$`Optimal Matrices`), "\n")
  cat("Total Number of Matrices with Same Trace:", result$`Count All Matrices with Same Trace`, "\n")
  
  for (i in seq_along(result$`Optimal Matrices`)) {
    cat("\nMatrix", i, ":\n")
    print(result$`Optimal Matrices`[[i]])
  }
}

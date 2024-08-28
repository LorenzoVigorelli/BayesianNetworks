library(dplyr)
library(readr)
library(bnlearn)
library(bnstruct)
library(Rgraphviz)
library(foreach)
library(doParallel)
library(ggplot2)
library(tidyr)
library(purrr)


########################################################
#               K2 Scoring Function(s)                 #
########################################################

## Scoring function for exact results
scoring_function <- function(data, x_i, parents){
  # Find the possible values of the attribute i
  r_i <- data |> distinct(data[[x_i]]) |> nrow()

  if (length(parents) == 0){
    N <- nrow(data)

    num1 <- data |> 
      group_by(data[[x_i]]) |> 
      count() |> 
      mutate(a = factorial(n)) |> 
      ungroup() |> 
      pull(a) |>
      prod()

    # To do consistently with the other case
  } else {
    alpha <- data |> group_by(data[c(x_i, parents)]) |> count()

    N <- alpha |> group_by(alpha[parents]) |> 
      summarise(N = sum(n), .groups = "drop") |> 
      select(N)

    num1 <- alpha |>
      group_by(alpha[parents]) |>
      summarise(alpha = prod(factorial(n)), .groups = "drop") |>
      pull(alpha)
  }

  # Calculate the score
  den <- sapply(N, function(x, r) factorial(x + r - 1), r=r_i)
  num2 <- factorial(r_i - 1)

  return(prod(num2 / den * num1))
}

## Log scoring function to avoid numerical issues
log_scoring_function <- function(data, x_i, parents){
  # Find the possible values of the attribute i
  r_i <- data |> distinct(data[[x_i]]) |> nrow()

  if (length(parents) == 0){
    N <- nrow(data)

    num1 <- data |> 
      group_by(data[[x_i]]) |> 
      count() |> 
      mutate(a = lfactorial(n)) |> 
      ungroup() |> 
      pull(a) |>
      sum()

    # To do consistently with the other case
  } else {
    alpha <- data |> group_by(data[c(x_i, parents)]) |> count()

    N <- alpha |> group_by(alpha[parents]) |> 
      summarise(N = sum(n), .groups = "drop") |> 
      select(N)

    num1 <- alpha |>
      group_by(alpha[parents]) |>
      summarise(alpha = sum(lfactorial(n)), .groups = "drop") |>
      pull(alpha)
  }

  # Calculate the score
  den <- sapply(N, function(x, r) lfactorial(x + r - 1), r=r_i)
  num2 <- lfactorial(r_i - 1)

  return(sum(num1 - den + num2))
}



########################################################
#               K2 Algorithm Implementation            #
########################################################


K2_algorithm <- function(data, max_parents){
  names <- colnames(data)
  results <- c()

  for (i in 1:ncol(data)){
    x_i <- names[i]
    parents <- c()
    p_old <- log_scoring_function(data, x_i, parents)
    proceed <- TRUE
    
    while (proceed){
      # Check if the maximum number of parents has been reached
      if (length(parents) >= max_parents){
        break
      }

      # Compute the predecessors
      predecessors <- setdiff(names[0:(i-1)], parents)
      if (length(predecessors) == 0){
        break
      }

      # Try adding a new parent
      scores <- sapply(predecessors, function(z) log_scoring_function(data, x_i, c(z, parents)))
      p_new <- max(scores)
      
      # If the score increases, add the parent
      if (p_new > p_old){
        p_old <- p_new
        parents <- c(parents, names[which.max(scores)])
      } else {
        proceed <- FALSE
      }
    } # end while
    
    results[[x_i]] <- parents
  } # end for

  return(list(names=names, parents_list=results))
}



########################################################
#               K2 Pipeline Implementation             #
########################################################


# Convert the parent-child relationships to a DAG
get_dag <- function(names, parents_list){
  dag <- empty.graph(names)

  # Add arcs based on the parent-child relationships
  for (child in names) {
    parents <- parents_list[[child]]
    if (length(parents) > 0) {
      for (parent in parents) {
        dag <- set.arc(dag, from = parent, to = child)
      }
    }
  }
  return(dag)
}


K2_to_dag <- function(data, max_parents){
  # Run the K2 algorithm
  results <- K2_algorithm(data, max_parents)
  names <- results$names
  parents_list <- results$parents_list

  # Convert the parent-child relationships to a DAG
  dag <- get_dag(names, parents_list)

  # Get the score of the DAG
  score <- score(dag, data)

  return(list(dag=dag, score=score))
}


K2_pipeline <- function(data, max_parents, max_iter, mode="local", n_cores=-1, return_history=FALSE){
  # Check if the mode is valid
  if (mode != "local" && mode != "parallel") {
    stop("Invalid mode. Please use 'local' or 'parallel'.")
  }

  # Return history only works in local mode
  if (return_history && mode == "parallel") {
    stop("Return history only works in local mode.")
  }
  
  # If the data does not have column names, assign them
  if (is.null(colnames(data))) {
    colnames(data) <- paste0("X", 1:ncol(data))
  }

  # Randomly shuffle data rows
  data <- data[sample(nrow(data)), ]

  # Try different random orders of the columns
  if (mode == "local") {
    history <- array(data = NA, dim = c(max_iter))

    # Initialize the best score and DAG
    score_best <- -Inf
    dag_best <- NULL

    # Try different random orders of the columns
    for (i in 1:max_iter){
      if (i > 1) {
        data <- data |> select(sample(colnames(data)))
      }
      result <- K2_to_dag(data, max_parents)
      score <- result$score

      # Update the best DAG
      if (score > score_best){
        score_best <- score
        dag_best <- result$dag
      }

      # Save the history
      history[i] <- score_best

    } 

    if (return_history) {
      return(list(dag=dag_best, score=score_best, history=history))
    }
    return(list(dag=dag_best, score=score_best))

  } else if (mode == "parallel") {
      # Setup parallel processing
      if (n_cores == -1) {
        n_cores <- detectCores()
      }
      cl <- makeCluster(n_cores)
      registerDoParallel(cl)

      # Try different random orders of the columns in parallel
      results <- mclapply(1:max_iter, function(i) {
        data_sampled <- data |> select(sample(colnames(data)))
        result <- K2_to_dag(data_sampled, max_parents)
        return(result)
      }, mc.cores = n_cores)

      # Initialize best score and DAG
      score_best <- -Inf
      dag_best <- NULL

      # Find the best result
      for (res in results) {
        if (res$score > score_best) {
          score_best <- res$score
          dag_best <- res$dag
        }
      }

      return(list(dag=dag_best, score=score_best))
    }

  return(NULL)
} 

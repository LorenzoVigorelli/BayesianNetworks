########################################################
#                  Introduction                        #
########################################################
# Authors: Paolo Lapo Cerni, Lorenzo Vigorelli, Arman Singh Bains
# Date: 28/09/2024
# Description: This script contains the functions used to implement the K2 algorithm.
# The K2 algorithm is a scoring-based algorithm used to learn the structure of a Bayesian network from data.


########################################################
#               Load Required Libraries                #
########################################################

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


########################################################
#               BNstruct function Function(s)          #
########################################################

learning <- function (data, maxParent, algo = "k2", percentage = 1, plot = F) {
    # Check if the percentage is within the valid range
    if (percentage < 0 || percentage > 1) {
        stop("Percentage must be between 0 and 1.")
    }

    # Subset the data based on the given percentage
    # Shuffle the rows of the data randomly
    num_rows <- nrow(data)
    data <- data[1:as.integer(percentage * num_rows), ]

    # Determine the starting value for the dataset
    minValue <- min(data, na.rm = TRUE)
    startsFrom <- ifelse(minValue == 0, 0, 1)

    # Calculate the number of unique values for each column in the data
    sizes <- sapply(data, function(x) length(unique(x)))
    sizes <- as.numeric(sizes)

    # Create a BNDataset object with the given data
    dataset <- BNDataset(data = data,
                         discreteness = rep(TRUE, ncol(data)),
                         variables = colnames(data), 
                         starts.from = startsFrom, 
                         node.sizes = sizes)
                    
    # Learn the network structure using the specified algorithm
    dag <- learn.network(algo = algo, x = dataset, max.parents = maxParent)

    # Create an empty graph and set its adjacency matrix
    net = empty.graph(names(data))
    amat(net) <- dag(dag)

    # Convert data columns to factors
    for (i in 1:length(names(data))) {
        name = names(data)[i]
        data[, name] = as.factor(as.character(data[, name]))
    }

    # Calculate the score of the network
    score <- score(net, data = data)

    # Plot the DAG if the plot parameter is TRUE
    if (plot) {
        plot(dag)
    }
        
    # Return the DAG and its score as a list
    return(list(dag = net, score = score))
}


########################################################
#             Compute shd and plots Function(s)        #
########################################################

# Function to compute the Structural Hamming Distance (SHD) between two networks
computeShdSingle <- function(theor, empir, plot=FALSE) {
  
  # Convert the theoretical and empirical networks to adjacency matrices
  DAG1 <- amat(theor)
  DAG2 <- amat(empir)
  
  # Compute the Structural Hamming Distance between the two adjacency matrices
  shd <- shd(DAG1, DAG2)
  
  # If plot is TRUE, visualize the comparison between the two networks
  if (plot) {
    graphviz.compare(theor, empir, shape="rectangle")
  }
  
  # Return the computed SHD
  return(shd)
}


########################################################
#     Compute shd multiples and plots Function(s)      #
########################################################

# Function to compute the Structural Hamming Distance (SHD) between a theoretical model and a list of empirical models
computeShd <- function(theor, empirList, plot=FALSE) {
  # Convert the theoretical model to an adjacency matrix
  DAG1 <- amat(theor)
  
  # Compute the SHD between the theoretical model and each empirical model in the list
  shd_values <- sapply(empirList, function(empir) shd(DAG1, amat(empir)))
  
  # If plot is TRUE, visualize the comparison between the models
  if (plot) {
    # Create a list of graphs to visualize, including the theoretical model and the empirical models
    plot_list <- c(list(theor), empirList)
    
    # Number of empirical models
    num_empirical <- length(empirList)
    
    # Main titles for the graphs
    titles <- c("THEORETICAL MODEL", paste("EMPIRICAL MODEL", seq_len(num_empirical)))
    # Subtitles for the graphs, including the SHD values
    subtitles <- c(paste("SHD =", "0"), paste("SHD =", shd_values))
    
    # Visualize the comparative graphs using graphviz.compare
    do.call(graphviz.compare, c(list(plot_list[[1]]), plot_list[-1], 
                                list(shape = "rectangle", 
                                     main = titles, 
                                     sub = subtitles,
                                     diff.args = list(tp.lwd = 2, tp.col = "green", fn.col = "orange"))))
  }
  
  # Return the computed SHD values
  return(shd_values)
}


########################################################
#             Compute scores and DAGs                  #
########################################################


# Function to compute scores and DAGs for a given algorithm, data, and percentages
compute_scores_and_dags <- function(algo, data, percentages) {
  # Apply the learning function to each percentage and store the results
  results <- lapply(percentages, function(p) learning(data = data, algo = algo, maxParent = 3, percentage = p))
  
  # Extract the scores from the results
  scores <- sapply(results, function(res) res$score)
  
  # Extract the DAGs from the results
  dags <- lapply(results, function(res) res$dag)
  
  # Normalize the scores
  min_score <- min(scores)
  max_score <- max(scores)
  scores_normalized <- (scores - min_score) / (max_score - min_score)
  
  # Return a list containing the normalized scores and the DAGs
  return(list(scores = scores_normalized, dags = dags))
}

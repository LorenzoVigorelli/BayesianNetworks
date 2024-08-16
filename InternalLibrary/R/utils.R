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
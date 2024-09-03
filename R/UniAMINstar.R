#' Asymptotic Unimax test for testing equality means against umbrella-ordered alternatives with unknown peak in one-way ANOVA
#' @export
#' @param sample_data list
#' @param significance_level numeric
#' @return Peak numeric
#' @return Critical value numeric
#' @return Test statistic value numeric
#' @return Result Character
#' @details Testing of H_0:mu_1 = mu_2 = ... = mu_k vs H_1:mu_1 <=.....<= mu_(h-1)<= mu_h >= mu_(h+1)>=....>= mu_k (at least one strict inequality, h is unknown), where mu_i represents the population means of the i-th treatment. The input consists of two variables: sample_data and significance_level. The output consists of the peak of the umbrella, the critical value, the UniAMAX* test statistic value, and the result, which indicates whether to reject or not reject the null hypothesis.
#' @importFrom Iso pava
#' @importFrom MASS mvrnorm
#' @importFrom stats qchisq
#' @importFrom stats quantile rnorm var
#' @author Subha Halder

UniAMINstar <- function(sample_data, significance_level){
  set.seed(456)
  R_MLE <- function(X, n) {
    X1 <- X[-1]
    n1 <- n[-1]
    sorted_indices <- order(X1)
    X1_sorted <- X1[sorted_indices]
    n1_sorted <- n1[sorted_indices]
    A <- numeric(length(X1_sorted))
    for (j in 2:length(X)) {
      A[j-1] <- (n[1] * X[1] + sum(n1_sorted[1:(j - 1)] * X1_sorted[1:(j - 1)])) /
        (n[1] + sum(n1_sorted[1:(j - 1)]))
    }
    if (all(X1 >= X[1])) {
      new_X <- X
    } else if (A[length(X)-2] >= X1_sorted[length(X)-1]) {
      X <- rep(A[length(X)-1], length(X))
      new_X <- X
    } else {
      comparisons <- logical(length(X1_sorted) - 1)
      comparisons1 <- logical(length(X1_sorted) - 1)
      stored_values <- numeric(0)
      for (k in 1:(length(X1_sorted) - 1)) {
        comparisons[k] <- A[k] < X1_sorted[k + 1]
        if(comparisons1[k] <- A[k] < X1_sorted[k + 1]) {
          for (s in 1:k) {
            stored_values[s] <- X1_sorted[s]
          }
          break
        }
      }
      selected_A_values <- A[comparisons]
      X[1] <- selected_A_values[1]
      for (l in 2:length(X)) {
        if (X[l] %in% stored_values) {
          X[l] <- selected_A_values[1]
        }
      }
      new_X <- X
    }
    return(new_X)
  }

  unimodal <- function(X, w, lmode) {
    if(lmode==length(X)){
      new_X <- Iso::pava(X,w)
    } else if (lmode==1){
      new_X <- -(Iso::pava(-X,w))
    } else{
      X1 <- X[1:(lmode-1)]
      X2 <- X[(lmode+1):length(X)]
      w1 <- w[1:(lmode-1)]
      w2 <- w[(lmode+1):length(X)]
      newX1 <- Iso::pava(X1,w1)
      newX2 <- -(Iso::pava(-X2,w2))
      Y <- c(X[lmode],newX1,newX2)
      w_new<- c(w[lmode],w1,w2)
      v <- -(R_MLE(-Y,w_new))
      max_v <- v[-1]
      new_X1 <- max_v[1:length(X1)]
      new_X2 <- max_v[(length(X1)+1):length(max_v)]
      new_X <- c(new_X1,v[1],new_X2)
    }
    return(new_X)
  }
  unimod_peak <- function(sample_data_list) {
    n <- sapply(sample_data_list, length)
    mu0 <- sapply(sample_data_list, mean)
    var0 <- sapply(1:length(sample_data_list), function(i) sum((sample_data_list[[i]] - mu0[i])^2) / n[i])
    w0 <- n / var0
    Uni_star <- list()
    res <- numeric()
    for (p in 1:length(mu0)) {
      repeat {
        new_mu0 <- unimodal(X = sapply(sample_data_list, mean), w = w0, lmode = p)
        new_var0 <- sapply(1:length(sample_data_list), function(i) sum((sample_data_list[[i]] - new_mu0[i])^2) / n[i])
        new_w0 <- n / new_var0
        if (max(abs(new_mu0 - mu0)) <= 0.0000001) {
          break  # Exit the loop if the difference is less than epsilon
        }
        w0 <- new_w0
        mu0 <- new_mu0
        var0 <- new_var0
      }
      Uni_star[[p]] <- new_mu0
      res[p] <- sum(new_w0*(sapply(sample_data_list, mean)-Uni_star[[p]])^2)
    }
    return(which.min(res))
  }
  sample_data <- lapply(sample_data, function(x) x[!is.na(x)])
  num_samples = 20000
  num_datasets <- length(sample_data)
  n <- sapply(sample_data, length)
  proportions <- n / sum(n)
  var_data <- sapply(1:num_datasets, function(j) var(sample_data[[j]]))
  h <- unimod_peak(sample_data)
  s <- num_datasets - 1
  b_sq <- numeric(s)
  for (p in 1:s) {
    b_sq[p] <- proportions[p] * var_data[p+1] + proportions[p+1] * var_data[p]
  }
  B <- diag(sqrt(b_sq))
  diag_elements <- b_sq
  if (h==1||h==length(num_datasets)){
    upper_diag <- numeric(s - 1)
    lower_diag <- numeric(s - 1)
    for (i in 1:(s-1)) {
      upper_diag[i] <- -sqrt(proportions[i] * proportions[i + 2]) * var_data[i + 1]
      lower_diag[i] <- upper_diag[i]
    }
    tridiag_matrix <- matrix(0, nrow = s, ncol = s)
    tridiag_matrix[row(tridiag_matrix) == col(tridiag_matrix)] <- diag_elements
    tridiag_matrix[row(tridiag_matrix) == col(tridiag_matrix) - 1] <- upper_diag
    tridiag_matrix[row(tridiag_matrix) == col(tridiag_matrix) + 1] <- lower_diag
  } else {
    if(num_datasets == 3){
      upper_diag <- numeric(s - 1)
      lower_diag <- numeric(s - 1)
      upper_diag[h - 1] <- sqrt(proportions[h - 1] * proportions[h + 1]) * var_data[h]
      lower_diag[h - 1] <- upper_diag[h - 1]
      tridiag_matrix <- matrix(0, nrow = s, ncol = s)
      tridiag_matrix[row(tridiag_matrix) == col(tridiag_matrix)] <- diag_elements
      tridiag_matrix[row(tridiag_matrix) == col(tridiag_matrix) - 1] <- upper_diag
      tridiag_matrix[row(tridiag_matrix) == col(tridiag_matrix) + 1] <- lower_diag
    } else if (num_datasets == 4){
      if (h==2) {
        upper_diag[1] <- sqrt(proportions[1] * proportions[3]) * var_data[2]
        lower_diag[1] <- upper_diag[1]
        upper_diag[2] <- -sqrt(proportions[2] * proportions[4]) * var_data[3]
        lower_diag[2] <- upper_diag[2]
        tridiag_matrix <- matrix(0, nrow = s, ncol = s)
        tridiag_matrix[row(tridiag_matrix) == col(tridiag_matrix)] <- diag_elements
        tridiag_matrix[row(tridiag_matrix) == col(tridiag_matrix) - 1] <- upper_diag
        tridiag_matrix[row(tridiag_matrix) == col(tridiag_matrix) + 1] <- lower_diag
      } else{
        upper_diag[1] <- -sqrt(proportions[1] * proportions[3]) * var_data[2]
        lower_diag[1] <- upper_diag[1]
        upper_diag[2] <- sqrt(proportions[2] * proportions[4]) * var_data[3]
        lower_diag[2] <- upper_diag[2]
        tridiag_matrix <- matrix(0, nrow = s, ncol = s)
        tridiag_matrix[row(tridiag_matrix) == col(tridiag_matrix)] <- diag_elements
        tridiag_matrix[row(tridiag_matrix) == col(tridiag_matrix) - 1] <- upper_diag
        tridiag_matrix[row(tridiag_matrix) == col(tridiag_matrix) + 1] <- lower_diag
      }
    } else{
      upper_diag <- numeric(s - 1)
      lower_diag <- numeric(s - 1)
      if (h==2){
        upper_diag[h-1] <- sqrt(proportions[h-1] * proportions[h+1]) * var_data[h]
        lower_diag[h-1] <- upper_diag[h-1]
        for (i in (h+1):(num_datasets-1)){
          upper_diag[i-1] <- -sqrt(proportions[i-1] * proportions[i + 1]) * var_data[i]
          lower_diag[i-1] <- upper_diag[i-1]
        }
      } else if (h==num_datasets-1){
        upper_diag[h-1] <- sqrt(proportions[h-1] * proportions[h+1]) * var_data[h]
        lower_diag[h-1] <- upper_diag[h-1]
        for (i in 1:(num_datasets-3)){
          upper_diag[i] <- -sqrt(proportions[i] * proportions[i + 2]) * var_data[i+1]
          lower_diag[i] <- upper_diag[i]
        }
      } else{
        for (i in 1:(h-2)){
          upper_diag[i] <- -sqrt(proportions[i] * proportions[i + 2]) * var_data[i+1]
          lower_diag[i] <- upper_diag[i]
        }
        upper_diag[h-1] <- sqrt(proportions[h-1] * proportions[h+1]) * var_data[h]
        lower_diag[h-1] <- upper_diag[h-1]
        for (i in h:(num_datasets-2)){
          upper_diag[i] <- -sqrt(proportions[i] * proportions[i + 2]) * var_data[i+1]
          lower_diag[i] <- upper_diag[i]
        }
      }
      tridiag_matrix <- matrix(0, nrow = s, ncol = s)
      tridiag_matrix[row(tridiag_matrix) == col(tridiag_matrix)] <- diag_elements
      tridiag_matrix[row(tridiag_matrix) == col(tridiag_matrix) - 1] <- upper_diag
      tridiag_matrix[row(tridiag_matrix) == col(tridiag_matrix) + 1] <- lower_diag
    }
  }
  D <- solve(B) %*% tridiag_matrix %*% solve(B)
  D_star_amin_values <- numeric(num_samples)
  for (k in 1:num_samples) {
    bootstrap_s <- MASS::mvrnorm(mu = rep(0, s), Sigma = D)
    D_star_amin_values[k] <- min(bootstrap_s)
  }
  quantile_value <- quantile(sort(D_star_amin_values), probs = 1-significance_level)
  if (h==1){
    Aunimin <- min(c(
      sapply(1:(num_datasets-1), function(i) {
        (mean(sample_data[[i]]) - mean(sample_data[[i + 1]])) /
          sqrt(
            (var(sample_data[[i]]) / length(sample_data[[i]])) +
              (var(sample_data[[i + 1]]) / length(sample_data[[i + 1]]))
          )
      })))
  } else if (h==num_datasets){
    Aunimin <- min(c(
      sapply(1:(num_datasets - 1), function(i) {
        (mean(sample_data[[i+1]]) - mean(sample_data[[i]])) /
          sqrt(
            (var(sample_data[[i+1]]) / length(sample_data[[i+1]])) +
              (var(sample_data[[i]]) / length(sample_data[[i]]))
          )
      })))
  } else {
    Aunimin <- min(c(
      sapply(2:h, function(i) {
        (mean(sample_data[[i]]) - mean(sample_data[[i - 1]])) /
          sqrt(
            (var(sample_data[[i]]) / length(sample_data[[i]])) +
              (var(sample_data[[i - 1]]) / length(sample_data[[i - 1]]))
          )
      }),
      sapply((h + 1):num_datasets, function(i) {
        (-mean(sample_data[[i]]) + mean(sample_data[[i - 1]])) /
          sqrt(
            (var(sample_data[[i]]) / length(sample_data[[i]])) +
              (var(sample_data[[i - 1]]) / length(sample_data[[i - 1]]))
          )
      })
    ))
  }
  if (Aunimin > quantile_value) {
    result <- "Reject null hypothesis"
  } else {
    result <- "Do not reject null hypothesis"
  }
  return(paste("Peak:", h,"; UniAMIN* Critical value:", quantile_value, "; UniAMIN* Test statistic:", Aunimin, "; Result:", result))
}

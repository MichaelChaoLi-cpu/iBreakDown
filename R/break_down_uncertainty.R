#' Explanation Level Uncertainty of Sequential Variable Attribution
#'
#' This function calculates the break down algorithm for \code{B} random orderings.
#' Then it calculates the distribution of attributions for these different orderings.
#' Note that the \code{shap()} function is just a simplified interface to the \code{break_down_uncertainty()} function
#' with a default value set to \code{B=25}.
#'
#' @param x an explainer created with function \code{\link[DALEX]{explain}} or a model.
#' @param data validation dataset, will be extracted from \code{x} if it is an explainer.
#' @param predict_function predict function, will be extracted from \code{x} if it is an explainer.
#' @param new_observation a new observation with columns that correspond to variables used in the model.
#' @param ... other parameters.
#' @param B number of random paths
#' @param keep_distributions if \code{TRUE} then we will keep distribution for predicted values. It's needed by the describe function.
#' @param path if specified, then this path will be highlighed on the plot. Use \code{average} in order to show an average effect
#' @param label name of the model. By default it's extracted from the 'class' attribute of the model.
#'
#' @return an object of the \code{break_down_uncertainty} class.
#' @importFrom utils head
#'
#' @seealso \code{\link{break_down}}, \code{\link{local_attributions}}
#'
#' @references Explanatory Model Analysis. Explore, Explain and Examine Predictive Models. \url{https://ema.drwhy.ai}
#'
#' @examples
#' library("DALEX")
#' library("iBreakDown")
#' set.seed(1313)
#' model_titanic_glm <- glm(survived ~ gender + age + fare,
#'                        data = titanic_imputed, family = "binomial")
#' explain_titanic_glm <- explain(model_titanic_glm,
#'                                data = titanic_imputed,
#'                                y = titanic_imputed$survived,
#'                            label = "glm")
#'
#' # there is no explanation level uncertanity linked with additive models
#' bd_glm <- break_down_uncertainty(explain_titanic_glm, titanic_imputed[1, ])
#' bd_glm
#' plot(bd_glm)
#'
#' \dontrun{
#' ## Not run:
#' library("randomForest")
#' set.seed(1313)
#' model <- randomForest(status ~ . , data = HR)
#' new_observation <- HR_test[1,]
#'
#' explainer_rf <- explain(model,
#'                         data = HR[1:1000, 1:5])
#'
#' bd_rf <- break_down_uncertainty(explainer_rf,
#'                            new_observation)
#' bd_rf
#' plot(bd_rf)
#'
#' # example for regression - apartment prices
#' # here we do not have intreactions
#' model <- randomForest(m2.price ~ . , data = apartments)
#' explainer_rf <- explain(model,
#'                         data = apartments_test[1:1000, 2:6],
#'                         y = apartments_test$m2.price[1:1000])
#'
#' bd_rf <- break_down_uncertainty(explainer_rf, apartments_test[1,])
#' bd_rf
#' plot(bd_rf)
#'
#' bd_rf <- break_down_uncertainty(explainer_rf, apartments_test[1,], path = 1:5)
#' plot(bd_rf)
#'
#' bd_rf <- break_down_uncertainty(explainer_rf,
#'                                      apartments_test[1,],
#'                                      path = c("floor", "no.rooms", "district",
#'                                          "construction.year", "surface"))
#' plot(bd_rf)
#'
#' bd <- break_down(explainer_rf,
#'                     apartments_test[1,])
#' plot(bd)
#'
#' s <- shap(explainer_rf,
#'                    apartments_test[1,])
#' plot(s)
#' }
#' @export
#' @rdname break_down_uncertainty
break_down_uncertainty <- function(x, ...,
                                   keep_distributions = TRUE,
                                   B = 10, clusterNumber = 2)
  UseMethod("break_down_uncertainty")

#' @export
#' @rdname break_down_uncertainty
break_down_uncertainty.explainer <- function(x, new_observation,
                       ...,
                       keep_distributions = TRUE,
                       B = 10, clusterNumber = 2) {
  # extracts model, data and predict function from the explainer
  model <- x$model
  data <- x$data
  predict_function <- x$predict_function
  label <- x$label

  break_down_uncertainty.default(model, data, predict_function,
                     new_observation = new_observation,
                     label = label,
                     ...,
                     keep_distributions = keep_distributions,
                     B = B, clusterNumber = clusterNumber)
}

get_single_random_path <- function(x, data, predict_function, new_observation,
                                   label, random_path, clusterNumber) {
  ###NOTE:
  ###ML: reduction in loops is extremely hard,
  ###we move foreach part to this function to reduce the time

  # if predict_function returns a single vector, conrvet it to a data frame
  #if (length(unlist(predict_function(x, new_observation))) > 1) {
  #  predict_function_df <- predict_function
  #} else {
  #  predict_function_df <- function(...) {
  #    tmp <- as.data.frame(predict_function(...))
  #    colnames(tmp) = label
  #    tmp
  #  }
  #} ###ML: test regressor or classifier; useless; drop
  predict_function_df <- predict_function

  vnames <- colnames(data)
  names(vnames) <- vnames
  current_data <- data

  yhats <- list()
  yhats[[1]] <- mean(predict_function_df(x, current_data))
  ###ML: revise colMeans to mean
  cl <- parallel::makeCluster(clusterNumber)
  doSNOW::registerDoSNOW(cl)
  opts <- list(progress = progress_fun)
  yhats.foreach <- foreach::foreach(i = seq_along(random_path), .combine = 'c',
                                    .packages = 'randomForest') %dopar% {
    candidate <- random_path[i]
    current_data[,candidate] <- new_observation[,candidate]
    mean(predict_function_df(x, current_data))
  } ###ML: numerous prediction
  parallel::stopCluster(cl)
  for(j in yhats.foreach) {
    yhats <- append(yhats, j)
  }

  diffs <- apply(do.call(rbind, yhats), 2, diff)
  if (is.vector(diffs)) { #93
    diffs <- t(diffs)
  }
  #76
  new_observation_vec <- sapply(as.data.frame(new_observation), nice_format) # same as in BD

  single_cols <- lapply(1:ncol(diffs), function(col) {

    variable_names <- vnames[random_path]
    data.frame(
      contribution = diffs[,col],
      variable_name = variable_names
    )
  })

  do.call(rbind,single_cols)
}

#' @export
#' @rdname break_down_uncertainty
break_down_uncertainty.default <- function(x, data, predict_function = predict,
                               new_observation,
                               label = class(x)[1],
                               ...,
                               path = NULL,
                               keep_distributions = TRUE,
                               B = 10, clusterNumber = 2) {
  # here one can add model and data and new observation
  # just in case only some variables are specified
  # this will work only for data.frames
  if ("data.frame" %in% class(data)) {
    common_variables <- intersect(colnames(new_observation), colnames(data))
    new_observation <- new_observation[1, common_variables, drop = FALSE]
    data <- data[,common_variables, drop = FALSE]
  } ###ML: keep only one observation; need to revise

  # Now we know the path, so we can calculate contributions
  # set variable indicators
  # start random path
  p <- ncol(data)
  previous_paths <- list()

  #### original version
  ###ML10.13: here we go back to original version

  result <- lapply(1:B, function(b) {
    # ensure each path was previously unused
    needs_new_path <- TRUE
    while (needs_new_path) {
      random_path <- sample(1:p)
      # was the path previously used?
      needs_new_path <- any(
        sapply(previous_paths, function(one_prev_path) {
          identical(one_prev_path, random_path)}))
    }
    previous_paths <- c(previous_paths, list(random_path))
    cat(b, " ")
    tmp <- get_single_random_path(x, data, predict_function,
                                  new_observation, label, random_path, clusterNumber)
    tmp
  })

  result <- do.call(rbind, result)

  result <- result[,c("contribution", "variable_name")]
  result <- aggregate(result$contribution, list(result$variable_name), FUN = mean)
  colnames(result) <- c("variable_name", "contribution")

  result <- as.data.frame(t(result))
  colnames(result) <- as.character(unlist(result[1,]))
  result <- result[-1,]

#  if(clusterNumber > 1){
#    #### foreach and doSNOW
#    # current version is for random forest only
#    cl <- parallel::makeCluster(clusterNumber)
#    doSNOW::registerDoSNOW(cl)
#    opts <- list(progress=progress_fun)
#    result <- foreach::foreach(b = seq(1,B, 1), .combine = 'rbind',
#                               .packages='randomForest', .options.snow=opts) %dopar% {
#                                 # ensure each path was previously unused
#                                 needs_new_path <- TRUE
#                                 while (needs_new_path) {
#                                   random_path <- sample(1:p)
#                                   # was the path previously used?
#                                   needs_new_path <- any(
#                                     sapply(previous_paths, function(one_prev_path) {
#                                       identical(one_prev_path, random_path)}))
#                                 }
#                                 previous_paths <- c(previous_paths, list(random_path))
#                                 tmp <- do.call(get_single_random_path,
#                                                list(x, data, predict_function,
#                                                     new_observation, label, random_path))
#                                 tmp$B <- b
#                                 tmp
#                               }
#    parallel::stopCluster(cl)
#
#    result <- result[,c("contribution", "variable_name")]
#    result <- aggregate(result$contribution, list(result$variable_name), FUN = mean)
#    colnames(result) <- c("variable_name", "contribution")
#
#    result <- as.data.frame(t(result))
#    colnames(result) <- as.character(unlist(result[1,]))
#    result <- result[-1,]
#  }

#  class(result) <- c("break_down_uncertainty", "data.frame")

#  if (keep_distributions) {
#    ## this yhats is not calculated like in breakDown
#    yhats <- list(NULL)
#
#    yhats_distribution <- calculate_yhats_distribution(x, data, predict_function, label, yhats)
#
#    attr(result, "yhats_distribution") <- yhats_distribution
#  }

#  target_yhat <- predict_function(x, new_observation)
#  yhatpred <- as.data.frame(predict_function(x, data))
#  baseline_yhat <- colMeans(yhatpred)

#  attr(result, "prediction") <- as.numeric(target_yhat)
#  attr(result, "intercept") <- as.numeric(baseline_yhat)

  result
}

#' @export
#' @rdname break_down_uncertainty
shap <- function(x, ..., B = 25, clusterNumber = 4) {
  ret <- break_down_uncertainty(x, ..., B = B, path = "average", clusterNumber = clusterNumber)

  class(ret) <- c("shap", class(ret))

  ret
}

progress_fun <- function(n){
  cat(n, ' ')
  if (n%%100==0){
    cat('\n')
  }
}

# Note: currently we just want to get the instance result.

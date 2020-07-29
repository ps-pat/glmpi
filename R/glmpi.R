modified_residuals <- function(object, center = TRUE)
    scale(residuals(object, type = "response") / sqrt(1 - hatvalues(object)),
          center = center, scale = FALSE)

res_boot <- function(object){
    ## Extract needed components from the model.
    model_frame <- model.frame(object)
    n <- nrow(model_frame)
    offset <- model.offset(model_frame)
    weights <- model.weights(model_frame)

    object_family <- family(object)
    link_inv <- object_family$linkinv

    model_matrix <- model.matrix(object)
    intercept_index <- which(attr(model_matrix, "assign") == 0)
    X <- if(!length(intercept_index))
             model_matrix
         else
             model_matrix[, -intercept_index, drop = FALSE]

    ## Sample error.
    epsilon_star <- sample(modified_residuals(object), replace = TRUE)

    ## Compute new response.
    y_star <- link_inv(predict(object, type = "link") + epsilon_star)

    ## Perform regression and extract coefficients.
    new_formula <- formula(paste("y_star",
                                 paste(formula(delete.response(terms(object))),
                                       collapse = " ")))
    glm(new_formula,
        data = object$data,
        offset = offset,
        weights = weights,
        family = object_family)
}

pred_boot <- function(object, newdata, R = 1000, M = 1){
    link_inv <- family(object)$linkinv

    boot_out <- res_boot(object)
    pred_out <- predict(boot_out, newdata, type = "response")

    pred_hat <- predict(object, newdata, type = "link")

    res <- replicate(n = R, simplify = FALSE, expr = {
        replicate(n = M, expr = {
            epsilon_star <- sample(modified_residuals(res_boot(boot_out)),
                                   size = nrow(newdata),
                                   replace = TRUE)
             pred_out - link_inv(pred_hat + epsilon_star)
        })
    })

    do.call(cbind, res)
}

#' @export
glmpi <- function(object, newdata, alpha = 0.05, R = 1000, M = 1){
    prediction <- predict(object, newdata, type = "response")
    predictions_boot<- pred_boot(object, newdata, R, M)

    lb <- alpha / 2
    ub <- 1 - lb
    quantiles <- apply(predictions_boot, 1, quantile,
                       probs = c(alpha / 2, 1 - alpha / 2),
                       names = FALSE)

    res <- t(matrix(prediction, 1) %x% matrix(c(1, 1), 2) + quantiles)
    colnames(res) <- paste("ic", c(lb, ub), sep = "")

    data.frame(pred = prediction, res)
}

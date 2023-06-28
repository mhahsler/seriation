#######################################################################
# seriation - Infrastructure for seriation
# Copyright (C) 2015 Michael Hahsler, Christian Buchta and Kurt Hornik
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


#' Register Seriation Based on a Variational Autoencoder
#'
#' Embeds the data or distance information in one dimension using a
#' variational autoencoder and returns the ordered points.
#'
#' Registers the method `"vae"` for [seriate()]. This method learns a
#' 1D embedding by training a variational autoencoder
#' (see Kingma and Welling, 2013) that
#' minimizes the reconstruction loss plus a penalty for
#' the Kullback–Leibler divergence of the
#' embedding space to a standard normal distribution.
#'
#' The method can embed a data matrix or a distance matrix.
#'
#' The chosen network topology for input data of size \eqn{n} is:
#'
#' **Encoder**
#'  - Input layer of size \eqn{n} accepts data \eqn{x}.
#'  - 1 shared dense layer of size \eqn{floor(n/2)} with relu activation.
#'      The two output layers are directly connected to this layer.
#'  - first dense output layer of size 1 which produces the embedding mean
#'      \eqn{\mu_z} with linear activation.
#'  - second dense output layer of size 1 for the log of the variance of \eqn{z},
#'      \eqn{log(\sigma_z)} used to sample
#'      with linear activation.
#'
#' **Sampler**
#'  - Samples a value \eqn{z \sim N(\mu_z, \sigma_z)}
#'    (optional: \eqn{\sigma_z} multiplied by `epsilon_std`).
#'
#' **Decoder**
#'  - input layer of size 1 (takes the embedding \eqn{z}).
#'  - 1 dense layers of size \eqn{floor(n/2)} and relu activation.
#'  - dense output layer of size \eqn{n} with activation specified in `output_activation`
#'    that reconstructs \eqn{x}.
#'    The default activation function is `"sigmoid"` which is suitable when all data
#'    in \eqn{x} is between 0 and 1.
#'    Use `"linear"` if the data can fall outside the \eqn{[0, 1]} range.
#'
#' **Loss function**
#'
#' The loss function is
#'
#' \deqn{L = loss(x, \hat{x}) + \beta\ KL(N(\mu_z,\sigma_z) || N(0,1)),}
#'
#' where the loss function is specified via the option `reconstruction_loss`.
#' The default for \eqn{\beta} is 1, but can be increased to create a beta-VAE.
#'
#' The encoder learns \eqn{p(z|x)}, the distribution in the embedding space given
#' the object in the input space. The decoder learns \eqn{p(x|z)} to reconstruct
#' \eqn{x} given a point \eqn{z} from the embedding space.
#' The prior for the embedding space is chosen to be \eqn{z \sim N(1,0)}. This
#' distribution is enforced by the KL divergence penalty in the loss function.
#'
#' Smoothness in the embedding space is created by the sampling process which
#' makes similar points in the embedding space decode as similar data and the
#' KL divergence which makes sure that points in the embedding space
#' are pushed together. These
#' properties create an embedding score whose order represents an approximate linear
#' order which can be used as a seriation.
#'
#' **Number of Epochs**
#' The number of epoch defaults to three times the number of observations. Early
#' stopping is used. This behavior can be changed via the `control` parameters.
#'
#'
#'
#' \bold{Note:} Package \pkg{keras} needs to be installed.
#'
#' @aliases register_vae vae
#' @family seriation
#' @returns Nothing.
#' @references
#' D. P. Kingma and M. Welling (2013). Auto-encoding variational bayes. ICLR.
#' \doi{10.48550/arXiv.1312.6114}
#'
#' Irina Higgins, Loic Matthey, Arka Pal, Christopher Burgess,
#' Xavier Glorot, Matthew Botvinick, Shakir Mohamed, Alexander Lerchner (2017).
#' beta-VAE: Learning Basic Visual Concepts with a Constrained
#' Variational Framework, International Conference on Learning Representations.
#'
#' @keywords optimize cluster
#' @examples
#' \dontrun{
#' register_vae()
#' get_seriation_method("matrix", "vae")
#'
#' data("Zoo")
#' Zoo[,"legs"] <- (Zoo[,"legs"] > 0)
#' x <- as.matrix(Zoo[,-17])
#' label <- rownames(Zoo)
#' class <- Zoo$class
#'
#' ### embed the rows in the data matrix
#' o <- seriate(x, method = "vae")
#'
#' pimage(x, o, prop = FALSE)
#'
#' # look at the class of the animals after ordering
#' class[get_order(o, 1)]
#'
#' # look at the embedding of the first dimension (rows)
#' attr(o[[1]], "vae")
#' attr(o[[1]], "vae")$get_embedding()
#' hist(attr(o[[1]], "vae")$get_embedding())
#'
#' ### embed the rows in the distance matrix (i.e.,
#' ###   the distance relationship between objects)
#' d <- dist(x, method = "binary")
#'
#' o <- seriate(d, method = "vae")
#' class[get_order(o)]
#' pimage(d, order = get_order(o))
#'
#' # get the embedding
#' sc <- attr(o[[1]], "vae")$get_embedding()
#' sc
#'
#' orderplot(sc, col = class)
#' legend("topright", legend = levels(class), col  = seq_along(levels(class)), pch = 16)
#' }
#' @export
register_vae <- function() {
  check_installed("keras")

  ## the code follows the example from:
  ##   https://github.com/rstudio/keras/blob/main/vignettes/examples/variational_autoencoder.R
  ## Maybe we can move to TF2:
  ##   https://blogs.rstudio.com/ai/posts/2018-10-22-mmd-vae/
  if (tensorflow::tf$executing_eagerly())
    tensorflow::tf$compat$v1$disable_eager_execution()

  .contr <- list(
    batch_size = 32L,
    epochs = NULL,
    optimizer = "sgd",
    epsilon_std = 1,
    scale_input = TRUE,
    output_activation = "sigmoid",
    reconstruction_loss = keras::loss_binary_crossentropy,
    beta = 1,
    early_stopping = TRUE,
    return_model = FALSE
  )

  .seriate_vae_rows <- function(x, control) {
    control <- .get_parameters(control, .contr)
    batch_size <- control$batch_size
    epochs <- control$epochs
    optimizer <- control$optimizer
    epsilon_std <- control$epsilon_std
    beta <- control$beta
    loss_fun <- control$reconstruction_loss
    output_activation <- control$output_activation

    x_train <- as.matrix(x)

    if (is.null(epochs))
      epochs <- 3 * nrow(x_train)

    if (control$scale_input)
      x_train <- x_train/max(x_train, na.rm = TRUE)

    # Parameters --------------------------------------------------------------
    original_dim <- ncol(x_train)
    latent_dim <- 1L
    intermediate_dim <- floor(original_dim / 2)

    # Model definition --------------------------------------------------------
    # input -> dense -> split dense layers for mean and log_var
    x <- keras::layer_input(shape = c(original_dim))
    h <- keras::layer_dense(x, intermediate_dim, activation = "relu")
    z_mean <- keras::layer_dense(h, latent_dim, activation = "linear")
    z_log_var <- keras::layer_dense(h, latent_dim, activation = "linear")

    # create a sample
    sampling <- function(arg) {
      z_mean <- arg[, 1:(latent_dim)]
      z_log_var <- arg[, (latent_dim + 1):(2 * latent_dim)]

      epsilon <- keras::k_random_normal(shape = c(keras::k_shape(z_mean)[[1]]),
                                 mean = 0.,
                                 stddev = epsilon_std)

      z_mean + keras::k_exp(z_log_var / 2) * epsilon
    }

    # note that "output_shape" isn't necessary with the TensorFlow backend
    # add sampling layer

    # pipes are only in R since 4.1.0
    #z <- keras::layer_concatenate(list(z_mean, z_log_var)) |>
    #  keras::layer_lambda(sampling)
    z <- keras::layer_lambda(keras::layer_concatenate(list(z_mean, z_log_var)), sampling)

    # we instantiate these layers separately so as to reuse them later
    decoder_h <-
      keras::layer_dense(units = intermediate_dim, activation = "relu")
    decoder_mean <-
      keras::layer_dense(units = original_dim, activation = output_activation)
    h_decoded <- decoder_h(z)
    x_decoded_mean <- decoder_mean(h_decoded)

    # end-to-end autoencoder
    vae <- keras::keras_model(x, x_decoded_mean, name = "VAE_end_to_end")
    # summary(vae)

    # encoder, from inputs to latent space
    encoder <- keras::keras_model(x, z_mean, name = "VAE_encoder")

    # decoder, from latent space to reconstructed inputs
    decoder_input <- keras::layer_input(shape = latent_dim)
    h_decoded_2 <- decoder_h(decoder_input)
    x_decoded_mean_2 <- decoder_mean(h_decoded_2)
    decoder <- keras::keras_model(decoder_input, x_decoded_mean_2, name = "VAE_decoder")

    vae_loss <- function(x, x_decoded_mean) {
      xent_loss <- (original_dim / 1.0) * loss_fun(x, x_decoded_mean)

      ## D_KL as the relative entropy between a diagonal multivariate normal, and a standard normal distribution
      ## see: https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
      ##  D_\text{KL}\left(
      ## \mathcal{N}\left(\left(\mu_1, \ldots, \mu_k\right)^\mathsf{T}, \operatorname{diag} \left(\sigma_1^2, \ldots, \sigma_k^2\right)\right) \parallel
      ## \mathcal{N}\left(\mathbf{0}, \mathbf{I}\right)
      ## \right) =
      ## {1 \over 2} \sum_{i=1}^k \left(\sigma_i^2 + \mu_i^2 - 1 - \ln\left(\sigma_i^2\right)\right).

      kl_loss <-
        .5 * keras::k_mean((keras::k_exp(z_log_var) +
                              keras::k_square(z_mean) -
                              1 - z_log_var), axis = -1L)

      xent_loss +  beta * kl_loss
    }

    keras::compile(vae, optimizer = optimizer, loss = vae_loss)

    # Model training ----------------------------------------------------------
    callback_early_stopping <- NULL
    if (control$early_stopping)
      callback_early_stopping <-
      list(
        keras::callback_early_stopping(
          monitor = "val_loss",
          patience = 10,
          restore_best_weights = TRUE,
          verbose = 1
        )
      )

    fit_vae <-
      function(epochs = 50L,
               view_metrics = TRUE,
               validation_split = .2) {
        keras::fit(
          vae,
          x_train,
          x_train,
          shuffle = TRUE,
          epochs = epochs,
          batch_size = batch_size,
          validation_split = validation_split,
          view_metrics = view_metrics,
          callbacks = callback_early_stopping
        )
      }

    orig_fit_result <- fit_vae(
      epochs = epochs,
      view_metric = "auto",
      validation_split = .2
    )

    get_embedding <- function() {
      sc <- drop(stats::predict(encoder, x_train))
      names(sc) <- rownames(x_train)

      sc
    }

    model <- list(
      end_to_end = vae,
      encoder = encoder,
      decoder = decoder,
      get_embedding = get_embedding,
      orig_fit_result = orig_fit_result
      #run_more_epochs = fit_vae
    )

    embedding <- as.vector(get_embedding())
    o <- order(embedding)
    attr(o, "embedding") <- embedding

    if (control$return_model)
      attr(o, "model") <- model

    # this only embeds rows
    o
  }

  .seriate_vae_matrix <- function(x, control, margin = seq_along(dim(x))) {
    if (1L %in% margin)
      row <- .seriate_vae_rows(x, control)
    else
      row <- NA

    if (2L %in% margin)
      col <- .seriate_vae_rows(t(x), control)
    else col <- NA

    list(row, col)
  }


  set_seriation_method(
    "matrix",
    "vae",
    .seriate_vae_matrix,
    "Use 1D variational autoencoder embedding of rows in a data matrix to create an order",
    .contr,
    randomized = TRUE
  )

  set_seriation_method(
    "dist",
    "vae",
    .seriate_vae_rows,
    "Use 1D variational autoencoder embedding of a distance matrix to create an order",
    .contr,
    randomized = TRUE
  )

}


generate_rhythm <- function (amps, phases, thetime, noise,
                             means = rep(1, length(amps))) {
  ## matrices with colums as time points
  ## cos(wt - phi) = cos(phi)cos(wt) + sin(phi)sin(wt)
  pure_sine_mat <-
    cos(matrix(phases/12*pi, ncol=1)) %*% cos(matrix(thetime/12*pi, nrow=1)) +
    sin(matrix(phases/12*pi, ncol=1)) %*% sin(matrix(thetime/12*pi, nrow=1))
  therun <-
    means + amps*pure_sine_mat +
    noise*matrix(rnorm(length(thetime)*length(amps)), nrow=length(amps))
  ## translation to base line (we want no negative values)
  replace(therun, therun < 0, 0)
}

generate_rhythm_mult <- function (amps, phases, thetime, noise,
                                  means = rep(1, length(amps))) {
  ## matrices with colums as time points
  ## cos(wt - phi) = cos(phi)cos(wt) + sin(phi)sin(wt)
  pure_sine_mat <-
    cos(matrix(phases/12*pi, ncol=1)) %*% cos(matrix(thetime/12*pi, nrow=1)) +
    sin(matrix(phases/12*pi, ncol=1)) %*% sin(matrix(thetime/12*pi, nrow=1))
  ## Exploit column-wise vector-matrix element-wise multiplication with
  ## recycling; e.g., (1:3) * matrix(1, 3, 6). Same for addition
  therun <-
    (means + amps*pure_sine_mat) *
    (1 + noise*matrix(rnorm(length(thetime)*length(amps)), nrow=length(amps)))
  ## translation to base line (we want no negative values)
  replace(therun, therun < 0, 0)
}

generate_rhythm_mult_outl <- function (amps, phases, thetime, noise, frac_outl,
                                       means = rep(1, length(amps))) {
  ## matrices with colums as time points
  ## cos(wt - phi) = cos(phi)cos(wt) + sin(phi)sin(wt)
  pure_sine_mat <-
    cos(matrix(phases/12*pi, ncol=1)) %*% cos(matrix(thetime, nrow=1)/12*pi) +
    sin(matrix(phases/12*pi, ncol=1)) %*% sin(matrix(thetime, nrow=1)/12*pi)
  ## exploit column-wise multiplication, e.g., (1:3) * matrix(1, 3, 6)
  ## add outliers with probab of 'frac_outl', half pos. and half neg. outliers
  ## outliers drawn from normal distr. with mean=of 10*noise
  therun <-
    (means + amps * pure_sine_mat) *
    (1 + matrix(rnorm(length(thetime)*length(amps),mean=0,sd=noise),
                nrow=length(amps))) +
    matrix(sample(c(0, -1, 1), prob=c((1 - frac_outl),
                                      frac_outl/2, frac_outl/2),
                  size=length(thetime)*length(amps), replace=TRUE) *
             rnorm(length(thetime)*length(amps), mean = 10*noise, sd = noise),
           ncol = length(thetime))
  ## translation to base line (we want no negative values)
  replace(therun, therun < 0, 0)
}

generate_fourier_wave_mult <- function (amps, phases, thetime, noise,
                                        means = rep(1, length(amps)),
                                        comps = 2/pi*(-1)^(1:10)/(1:10)) {

  ## Default sawtooth example with n=10:
  ## (2*A/pi) * sum_k=1^n (-1)^k*sin(k*w*t)/k
  ## Note that we use sine waves; to bring the phase to the same definition
  ## (time of max) as with the common cosine formulation, add a phase shift 6
  time_t <- t(t(matrix(rep(thetime, length(phases)), ncol = length(phases))) -
                phases) + 6
  ## outer() produces a 3D array, with component number on 1st, time on 2nd and
  ## sample on 3rd dimensions.
  ## colSums() slices along the 3rd dimension, produces colums sums of the
  ## 1+2dim slices, places these colum sums as column vectors in return matrix.
  ## Thus, the return matrix will have 1 dimension of input 2nd, and 2 dimension
  ## of input 3 dimension (1st dimension is summed away.). So, time becomes 1d,
  ## and sample becomes 2d.
  ## Finally, multiplication of 3D array with vectors works first along 1d,
  ## then along 2d, then along 3d, as expected.
  ## Final transpose gets time on 2d instead of 1d.
  waves <- t(colSums(comps*sin(pi/12*outer((1:length(comps)), time_t))))

  therun <- (means + amps*waves) *
    (1 + noise*matrix(rnorm(length(thetime)*length(amps)), nrow=length(amps)))

  replace(therun, therun < 0, 0)

}

generate_fourier_wave_mult_outl <- function (amps, phases, thetime,
                                             noise, frac_outl,
                                             means = rep(1, length(amps)),
                                             comps = 2/pi*(-1)^(1:10)/(1:10)) {

  ## Default sawtooth example with n=10:
  ## (2*A/pi) * sum_k=1^n (-1)^k*sin(k*w*t)/k
  ## Note that we use sine waves; to bring the phase to the same definition
  ## (time of max) as with the common cosine formulation, add a phase shift 6
  time_t <- t(t(matrix(rep(thetime, length(phases)), ncol = length(phases))) -
                phases) + 6
  ## outer() produces a 3D array, with component number on 1st, time on 2nd and
  ## sample on 3rd dimensions.
  ## colSums() slices along the 3rd dimension, produces colums sums of the
  ## 1+2dim slices, places these colum sums as column vectors in return matrix.
  ## Thus, the return matrix will have 1 dimension of input 2nd, and 2 dimension
  ## of input 3 dimension (1st dimension is summed away.). So, time becomes 1d,
  ## and sample becomes 2d.
  ## Finally, multiplication of 3D array with vectors works first along 1d,
  ## then along 2d, then along 3d, as expected.
  ## Final transpose gets time on 2d instead of 1d.
  waves <- t(colSums(comps*sin(pi/12*outer((1:length(comps)), time_t))))

  therun <- (means + amps*waves) *
    (1 + noise*matrix(rnorm(length(thetime)*length(amps)), nrow=length(amps))) +

    matrix(sample(c(0, -1, 1), prob=c((1 - frac_outl),
                                      frac_outl/2, frac_outl/2),
                  size=length(thetime)*length(amps), replace=TRUE) *
             rnorm(length(thetime)*length(amps), mean = 10*noise, sd = noise),
           ncol = length(thetime))

  replace(therun, therun < 0, 0)

}






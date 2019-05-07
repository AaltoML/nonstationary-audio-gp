# End-to-End Probabilistic Inference for Nonstationary Audio Analysis

https://arxiv.org/abs/1901.11436

Paper accepted to International Conference on Machine Learning (ICML) 2019.

We present fully probabilistic joint inference in the Gaussian time-frequency analysis and NMF model (GT-NMF). This is a nonstationary spectral mixture Gaussian process model, with a GP hyper-prior over the amplitude of each frequency channel output.

- We perform inference via expectation propagation in the Kalman filter setting, and compare this method to the extended Kalman filter.

- We implement the infinite-horizon GP method, which scales linearly in time and quadratically in the state dimensionality, and frees up significant amounts of memory allowing us to process signals with hundreds of thousands of data points.

- We show how the generative model can be applied to missing data synthesis, denoising and source separation without modification.




matlab/ folder contains the code and example scripts.


matlab/experiments/ folder allows you to rerun the experiments from the paper and produce the plots.





#### Reference:
```
@inproceedings{wilkinson2019end,
       title = {End-to-End Probabilistic Inference for Nonstationary Audio Analysis},
      author = {Wilkinson, William J. and Andersen, Michael Riis and Reiss, Joshua D. and Stowell, Dan and Solin, Arno},
        year = {2019},
   booktitle = {International Conference on Machine Learning (ICML)},
      volume = {97},
      series = {PMLR}
}
```



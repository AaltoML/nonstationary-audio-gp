# End-to-End Probabilistic Inference for Nonstationary Audio Analysis

Codes for the paper:

* William J. Wilkinson, Michael Riis Andersen, Joshua D. Reiss, Dan Stowell, and Arno Solin (2019). **End-to-End Probabilistic Inference for Nonstationary Audio Analysis**. Accepted for publication in *International Conference on Machine Learning (ICML)*. [[arXiv]](https://arxiv.org/abs/1901.11436)

We present fully probabilistic joint inference in the Gaussian time-frequency analysis and NMF model (GT-NMF). This is a nonstationary spectral mixture Gaussian process model, with a GP hyper-prior over the amplitude of each frequency channel output.

- We perform inference via expectation propagation in the Kalman filter setting, and compare this method to the extended Kalman filter.

- We implement the infinite-horizon GP method, which scales linearly in time and quadratically in the state dimensionality, and frees up significant amounts of memory allowing us to process signals with hundreds of thousands of data points.

- We show how the generative model can be applied to missing data synthesis, denoising and source separation without modification.


## Getting started

* The `matlab/` folder contains the code and example scripts.

* The `matlab/experiments/` folder allows you to rerun the experiments from the paper and produce the plots.


## Reference:
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

## License

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

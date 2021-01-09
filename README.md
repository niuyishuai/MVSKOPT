# MVSKOPT
 A MATLAB toolbox of DC programming approaches for solving the higher-order moment Mean-Variance-Skewness-Kurtosis (MVSK) portfolio optimization model

This project is supported by the National Natural Science Foundation of China (Grant No: 11601327).

## Installation
  1. Install POLYLAB toolbox (see [Polylab](https://github.com/niuyishuai/Polylab)) and DCAM toolbox (see [DCAM](https://github.com/niuyishuai/DCAM))
  2. Download the package to a local folder or by running:
```console
git clone https://github.com/niuyishuai/MVSKOPT
```
  3. Run Matlab and navigate to the code folder, then run `install.m` script to install the package.

## Examples
  See example `test_udca_ubdca.m` for DCA with commonly used universal DC decomposition and the associated boosted-DCA.

## Citation

```
@article{niu2011efficient,
  title={An efficient DC programming approach for portfolio decision with higher moments},
  author={Pham, Dinh Tao and Niu, Yi-Shuai},
  journal={Computational Optimization and Applications},
  volume={50},
  number={3},
  pages={525--554},
  year={2011},
  publisher={Springer}
}

@article{niu2020higherorder,
      title={Higher-order Moment Portfolio Optimization via The Difference-of-Convex Programming and Sums-of-Squares}, 
      author={Yi-Shuai Niu and Ya-Juan Wang},
      year={2020},
      eprint={1906.01509},
      archivePrefix={arXiv},
      primaryClass={math.OC}
}
```  

## License

Released under MIT license

## Contact

niuyishuai@sjtu.edu.cn

# MaxSq-basis


Features
============

* The reconstruction of broadband modes with maximal possible squeezing from the given covariance matrix. For details see <https://arxiv.org/abs/2403.05259>.


Run example:
============
To run example:
```    $ python3 run_example.py ```


Documentation:
============

* The $\hbar=2$ convention is used (vacuum quadrature variance $=1$). 

* The function *build_MSq_basis* in **msq_funcs.py** returns the basis of maximal possible squeezing. The order of sorting corresponds to the increasing minimal quadrature variance.

* The shortcoming of the function *build_MSq_basis* in **msq_funcs.py** is the order of vectors: the occupied but unsqueezed modes are last calculated. For large matrices, it is convenient and efficient to change the order of modes and find these unsqueezed modes earlier.
This is realized in the function build_MSq_basis_swap, where the parameter swap_cutoff determines the criteria for the sorting change.



Citing:
============

If you are using the decomposition for broadband modes with maximal possible squeezing, please cite the <https://arxiv.org/abs/2403.05259>:

Denis A. Kopylov, Torsten Meier, Polina R. Sharapova; Theory of Multimode Squeezed Light Generation in Lossy Media, arXiv:2403.05259 (2024)


QmeQ: Quantum master equation for Quantum dot transport calculations
====================================================================

QmeQ is an open-source Python package for calculations of transport through
quantum  dot devices. The so-called Anderson-type models are used to describe
the quantum dot device, where quantum dots are coupled to the leads by
tunneling. QmeQ can calculate the stationary state **particle** and
**energy currents** using various approximate density matrix approaches. As for
now we have implemented the following first-order methods

* Pauli (classical) master equation
* Lindblad approach
* Redfield approach
* First order von Neumann (1vN) approach

which can describe the effect of Coulomb blockade. Additionally, there is a
possibility to model electron-phonon interaction inside a quantum dot using
the first-order approaches. QmeQ also has two second-order methods

* Second order von Neumann (2vN) approach
* Real Time Diagrammatic (RTD) approach

2vN and RTD approaches can address the effects of cotunneling and pair tunneling.
Additionally, the 2vN approach can describe broadening of quantum dot states.
The advantage of RTD approach is that it requires a lot less memory and
computation time resources.

Physics disclaimer
------------------

All the methods in QmeQ are approximate so depending on parameter regime they
**can fail**, and a good knowledge of the method is required whether to trust
the result or not. For example, Redfield, 1vN, and 2vN approaches can **violate
positivity** of the reduced density matrix and lead to **currents flowing against
the bias**. We still think it is important to have a package where a user can
duplicate existing calculations, check applicability of different methods, or
simply discover new kind of physics using different approximate master equations.

Prokop Notes
------------

Activate environment (because Ubuntu 24 LTS does not allow to use pip directly):
```
source /home/prokop/venvs/ML/bin/activate
```
Compile like this ( using modified setup.py from Vlado ):
```
python setup.py install --root=/home/prokop/bin/
```

anyway the package is installed to :
```
/home/prokop/bin/home/prokop/venvs/ML/lib/python3.12/site-packages/qmeq/
```

you must update this path inside the script `spinless_trimer1.py`
```
path.insert(0, '/home/prokop/bin/home/prokop/venvs/ML/lib/python3.12/site-packages/qmeq/')
```

it does not run with cython (fallback to numpy), this is actually good for debugging
```
(ML) prokop@GTX3090:~/git_SW/qmeq$ python spinless_trimer1.py
WARNING: Cannot import Cython compiled modules for the special functions (specfunc.__init__.py).
WARNING: Cannot import Cython compiled modules for the approaches (builder_base.py).
WARNING: Cannot import Cython compiled modules for the approaches (builder_elph.py).
# spinless_trimer1.py start, Thu Dec 19 18:29:39 2024
# eps1: -10.0 eps2: -10.0 eps3: -10.0 t: 0.0 U: 220.0 W: 20.0
# GammaS: 0.2 GammaT: 0.05 VS: 0.2523 VT: 0.1262
# muS: 0.0000 muT: 0.0000 Temp: 0.2240 DBand: 1000.0
# VBiasMin: 0.0 VBiasMax: 60.0 dVBias: 0.100 NPoints: 600

```


NOTES:
------

#### Where to start navidation the code:

* when we call `qmeq.Builder()` we go to `BuilderBase.__init__()` in `/home/prokop/git_SW/qmeq/qmeq/builder/builder_base.py`


When we run it we get following flow:
```
DEBUG: Approach.solve()
DEBUG: Approach.prepare_solver()
DEBUG: ApproachPauli.generate_fct()
DEBUG: ApproachPauli.generate_kern() ncharge: 4  statesdm: [[0], [1, 2, 4], [3, 5, 6], [7], []]
DEBUG: ApproachPauli.generate_coupling_terms() b: 0  bp: 0  bcharge: 0
DEBUG: ApproachPauli.generate_coupling_terms() b: 1  bp: 1  bcharge: 1
DEBUG: ApproachPauli.generate_coupling_terms() b: 2  bp: 2  bcharge: 1
DEBUG: ApproachPauli.generate_coupling_terms() b: 4  bp: 4  bcharge: 1
DEBUG: ApproachPauli.generate_coupling_terms() b: 3  bp: 3  bcharge: 2
DEBUG: ApproachPauli.generate_coupling_terms() b: 5  bp: 5  bcharge: 2
DEBUG: ApproachPauli.generate_coupling_terms() b: 6  bp: 6  bcharge: 2
DEBUG: ApproachPauli.generate_coupling_terms() b: 7  bp: 7  bcharge: 3
DEBUG: ApproachPauli.generate_kern() kh.kern:
 [[-1.318  0.     0.     0.     0.     0.     0.     0.   ]
 [ 0.5   -0.018  0.     0.     0.4    0.4    0.     0.   ]
 [ 0.409  0.    -0.509  0.     0.     0.     0.4    0.   ]
 [ 0.409  0.     0.    -0.509  0.     0.     0.4    0.   ]
 [ 0.     0.009  0.     0.5   -0.409  0.     0.     0.4  ]
 [ 0.     0.009  0.5    0.     0.    -0.409  0.     0.4  ]
 [ 0.     0.     0.009  0.009  0.     0.    -0.9    0.4  ]
 [ 0.     0.     0.     0.     0.009  0.009  0.1   -1.2  ]]
DEBUG: Approach.solve_kern()
DEBUG: ApproachPauli.generate_current()
```


Installation
------------


For installation instructions see [INSTALL.md](INSTALL.md).

Tutorial & Examples
-------------------

For an introduction to QmeQ see this [tutorial][tutorial]
and various [examples][examples].

License
-------

QmeQ has [The BSD 2-Clause License][license] and it can be found
in [LICENSE.md](LICENSE.md).

Citing QmeQ
-----------

Please consider citing QmeQ if the use of this project gives results which lead
to scientific publication:

G. Kiršanskas, J. N. Pedersen, O. Karlström, M. Leijnse, and A. Wacker,
*QmeQ 1.0: An open-source Python package for calculations of transport through
quantum dot devices*, [Comput. Phys. Commun. 221, 317 (2017)][qmeqdoi].

The preprint version of the paper can be found on the
[arXiv.org][qmeqarxiv] server.

[tutorial]: https://github.com/gedaskir/qmeq-examples/tree/master/tutorial/tutorial.ipynb
[examples]: https://github.com/gedaskir/qmeq-examples
[license]: https://opensource.org/licenses/BSD-2-Clause
[qmeqdoi]: https://dx.doi.org/10.1016/j.cpc.2017.07.024
[qmeqarxiv]: https://arxiv.org/abs/1706.10104

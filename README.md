# MetabEngine

[![Generic badge](https://img.shields.io/badge/metabengine-ver_0.0.15-%3CCOLOR%3E.svg)](https://github.com/Waddlessss/metabengine/)
![Maintainer](https://img.shields.io/badge/maintainer-Huaxu_Yu-blue)
[![PyPI Downloads](https://img.shields.io/pypi/dm/bago.svg?label=PyPI%20downloads)](https://pypi.org/project/metabengine/)

**metabengine** is an integrated Python package for liquid chromatography-mass spectrometry (LC-MS) data processing.

* **Documentation:** https://metabengine.readthedocs.io/en/latest/
* **Source code:** https://github.com/Waddlessss/metabengine/
* **Bug reports:** https://github.com/Waddlessss/metabengine/issues/

It provides:

* Ion-identity-based feature detection (i.e. peak picking).
* Peak quality evaluation via artificial neural network.
* Accurate annotation of isotopes, adducts, and in-source fragments.
* Ultra-fast annotation of MS/MS spectra supported by [Flash Entropy Search](https://github.com/YuanyueLi/MSEntropy)

## Installation

```sh
pip install metabengine
```
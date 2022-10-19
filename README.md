# eXplain Your MOLecule (XYMOL)

[![Python](https://img.shields.io/badge/Python-3.7+-blue.svg)](https://www.python.org)
[![github ci](https://github.com/smu-tao-group/xymol/actions/workflows/ci.yml/badge.svg)](https://github.com/smu-tao-group/xymol/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/smu-tao-group/xymol/branch/main/graph/badge.svg?token=fl4kUOywR3)](https://codecov.io/gh/smu-tao-group/xymol)

XYMOL: A Python package to understand and explain atom/bond contributions of small molecules in machine learning models.

## Install

Using `pip` to install:

```
# for release (stable) version
pip install xymol

# for the latest version
pip install git+https://github.com/smu-tao-group/xymol.git
```

## Usage

The easiest way to use XYMOL is to input your featurizer (function to featurize SMILES) and the trained ML model through `create_map` function.

```python
from xymol import XYMOL

SMILES = "CCC" # replace with your SMILES
xymol = XYMOL(SMILES)
xymol.create_map(FEATURIZER, TRAINED_MODEL)
```

One example is displayed below. Green color means dropping this atom would lead to an increase in prediction, and vice versa.

<img src="./examples/imp.png" width="400">

## License

Apache-2.0 license

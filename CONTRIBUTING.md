<!--- This file based on https://github.com/PyGithub/PyGithub/blob/master/CONTRIBUTING.md --->
# Contributing

`PyMassSpec` uses `tox` to automate testing, autoformatting and packaging.

For more information on installing `tox` see https://tox.readthedocs.io/en/latest/install.html

## Coding style

`PyMassSpec` uses `yapf` for code formatting and `isort` to sort imports.

`yapf` and `isort` can be run via `tox`:
```bash
tox -e yapf
tox -e isort
```

To run the complete autoformatting suite run:
```bash
tox -e qa
```

## Automated tests

Tests are run with `tox` and `pytest`. To run tests for a specific Python version, such as Python 3.6, run:
```bash
tox -e py36
```

To run tests for all Python versions, simply run:
```bash
tox
```

## Type Annotations

Type annotations are checked using `mypy`. Run `mypy` using `tox`:
```bash
tox -e mypy
```

## Build documentation locally

Documentation is powered by Sphinx. A local copy of the documentation can be built with `tox`:
```bash
tox -e docs
```

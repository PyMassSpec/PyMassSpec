Overview
---------

.. This file based on https://github.com/PyGithub/PyGithub/blob/master/CONTRIBUTING.md

``PyMassSpec`` uses `tox <https://tox.readthedocs.io>`_ to automate testing and packaging,
and `pre-commit <https://pre-commit.com>`_ to maintain code quality.

Install ``pre-commit`` with ``pip`` and install the git hook:

.. prompt:: bash

	python -m pip install pre-commit
	pre-commit install


Coding style
--------------

`formate <https://formate.readthedocs.io>`_ is used for code formatting.

It can be run manually via ``pre-commit``:

.. prompt:: bash

	pre-commit run formate -a


Or, to run the complete autoformatting suite:

.. prompt:: bash

	pre-commit run -a


Automated tests
-------------------

Tests are run with ``tox`` and ``pytest``.
To run tests for a specific Python version, such as Python 3.10:

.. prompt:: bash

	tox -e py310


To run tests for all Python versions, simply run:

.. prompt:: bash

	tox

A series of reference images for ``test_Display.py`` are in the "tests/baseline" directory.
If these files need to be regenerated, run the following command:

.. prompt:: bash

	pytest --mpl-generate-path="tests/baseline" tests/test_Display.py


Type Annotations
-------------------

Type annotations are checked using ``mypy``. Run ``mypy`` using ``tox``:

.. prompt:: bash

	tox -e mypy



Build documentation locally
------------------------------

The documentation is powered by Sphinx. A local copy of the documentation can be built with ``tox``:

.. prompt:: bash

	tox -e docs

#!/usr/bin/python3

# This file is managed by `git_helper`. Don't edit it directly

# stdlib
import platform
import sys
from io import StringIO

# 3rd party
import rst2txt
from docutils.core import publish_file

# this package
from __pkginfo__ import *  # pylint: disable=wildcard-import

recipe_dir = repo_root / "conda"

if not recipe_dir.exists():
	recipe_dir.mkdir()

# TODO: entry_points, manifest

all_requirements = install_requires[:]

if isinstance(extras_require, dict):
	for requires in extras_require.values():
		all_requirements += requires

all_requirements = set(x.replace(" ", '') for x in set(all_requirements))
requirements_block = "\n".join(f"    - {req}" for req in all_requirements if req)

# txt_readme = publish_file(source=StringIO(long_description), writer=rst2txt.Writer())
# description_block = "\n".join([line.replace('"', '\\"') for line in txt_readme.split("\n")])
description_block = conda_description.replace('"', '\\"')

with open(recipe_dir / "meta.yaml", "w") as fp:
	fp.write(f"""{{% set name = "{pypi_name}" %}}
{{% set version = "{__version__}" %}}

package:
  name: "{{{{ name|lower }}}}"
  version: "{{{{ version }}}}"

source:
  url: "https://pypi.io/packages/source/{{{{ name[0] }}}}/{{{{ name }}}}/{{{{ name }}}}-{{{{ version }}}}.tar.gz"

build:
#  entry_points:
#    - {import_name} = {import_name}:main
#  skip_compile_pyc:
#    - "*/templates/*.py"          # These should not (and cannot) be compiled
  noarch: python
  script: "{{{{ PYTHON }}}} -m pip install . -vv"

requirements:
  build:
    - python
    - setuptools
    - wheel
  host:
    - pip
    - python
{requirements_block}
  run:
    - python
{requirements_block}

test:
  imports:
    - {import_name}

about:
  home: "{web}"
  license: "{__license__}"
  # license_family: LGPL
  # license_file: requirements.txt
  summary: "{short_desc}"
  description: "{description_block}"
  doc_url: {project_urls["Documentation"]}
  dev_url: {project_urls["Source Code"]}

extra:
  maintainers:
    - {author}
    - github.com/{github_username}

""")

print(f"Wrote recipe to {recipe_dir / 'meta.yaml'}")
#
# plat = platform.system().lower()
# arch = platform.architecture()[0][:2]
#
# if plat == "linux":
# 	conda_arch = f"linux-{arch}"
# elif plat == "windows":
# 	conda_arch = f"win-{arch}"
# elif plat == "darwin":
# 	conda_arch = f"osx-{arch}"
# else:
# 	sys.exit(1)
#
# with open(recipe_dir / "conda_arch.sh", "w") as fp:
# 	fp.write(f'#!/bin/bash\necho "{conda_arch}"')

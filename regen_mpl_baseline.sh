#!/bin/bash

tox -e py36-display -- --mpl-generate-path=tests/baseline

tox -e py36-display -- --mpl-generate-hash-library=image_hashes_36.json
tox -e py37-display -- --mpl-generate-hash-library=image_hashes_37.json
tox -e py38-display -- --mpl-generate-hash-library=image_hashes_38.json
tox -e py39-display -- --mpl-generate-hash-library=image_hashes.json
tox -e py310-display -- --mpl-generate-hash-library=image_hashes_313.json

python3 sort_hashes_json.py

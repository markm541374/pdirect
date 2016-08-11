#!/usr/bin/env bash

(cd pdirect/core; cython -a *.pyx)

pip install --upgrade --user -e .

mv *.so pdirect/core/
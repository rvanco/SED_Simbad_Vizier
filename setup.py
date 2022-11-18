#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 16:45:58 2022

@author: rvancoellie
"""

from setuptools import setup
from Cython.Build import cythonize


setup(
    name = "function", # project module
    ext_modules=cythonize("cython_code.pyx", compiler_directives={"language_level": "3"}), # cython file to compile
)
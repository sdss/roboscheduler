# encoding: utf-8
#
# setup.py
#

from setuptools import setup, Extension, find_packages

import os
import argparse
import sys


class getPybindInclude(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked.
    https://github.com/pybind/python_example/blob/master/setup.py
    """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)


def getIncludes():
    return [
        'include',
        getPybindInclude(),

    ]


sources = [
    'cextern/cCadenceCore.cpp',
    'cextern/cadenceCore.cpp'
]

extra_compile_args = ["--std=c++11", "-fPIC", "-v", "-O3"]
extra_link_args = None
if sys.platform == 'darwin':
    extra_compile_args += ['-stdlib=libc++', '-mmacosx-version-min=10.9']
    extra_link_args = ["-v", '-mmacosx-version-min=10.9']

module = Extension(
    'roboscheduler.cCadenceCore',
    include_dirs=getIncludes(),
    extra_compile_args=extra_compile_args,
    extra_link_args=extra_link_args,
    sources=sources
)

setup(ext_modules=[module])

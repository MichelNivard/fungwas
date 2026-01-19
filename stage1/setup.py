"""
Build script for Stage 1 C++ extension using pybind11.

Usage:
    python setup.py build_ext --inplace
    
Or for development:
    pip install -e .
"""

import os
import sys
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

class get_pybind_include:
    """Helper class to determine pybind11 include path."""
    def __str__(self):
        import pybind11
        return pybind11.get_include()

def get_armadillo_include():
    """Find Armadillo include path from conda."""
    conda_prefix = os.environ.get('CONDA_PREFIX', '')
    if conda_prefix:
        arma_path = os.path.join(conda_prefix, 'include')
        if os.path.exists(os.path.join(arma_path, 'armadillo')):
            return arma_path
    # Fallback to common locations
    for path in ['/usr/include', '/usr/local/include']:
        if os.path.exists(os.path.join(path, 'armadillo')):
            return path
    return None

def get_armadillo_lib():
    """Find Armadillo library path from conda."""
    conda_prefix = os.environ.get('CONDA_PREFIX', '')
    if conda_prefix:
        lib_path = os.path.join(conda_prefix, 'lib')
        if os.path.exists(lib_path):
            return lib_path
    return None

ext_modules = [
    Extension(
        'fungwas_stage1._stage1_cpp',
        sources=['src/stage1_pybind.cpp'],
        include_dirs=[
            get_pybind_include(),
            get_armadillo_include(),
        ],
        library_dirs=[get_armadillo_lib()] if get_armadillo_lib() else [],
        runtime_library_dirs=[get_armadillo_lib()] if get_armadillo_lib() else [],
        libraries=['armadillo'],
        extra_compile_args=['-std=c++14', '-O3', '-fopenmp', '-fPIC'],
        extra_link_args=['-fopenmp'],
        language='c++'
    ),
]

class BuildExt(build_ext):
    """Custom build_ext command with compiler customization."""
    def build_extensions(self):
        # Remove -Wstrict-prototypes (not valid for C++)
        if '-Wstrict-prototypes' in self.compiler.compiler_so:
            self.compiler.compiler_so.remove('-Wstrict-prototypes')
        super().build_extensions()

setup(
    name='fungwas-stage1',
    version='0.1.0',
    ext_modules=ext_modules,
    cmdclass={'build_ext': BuildExt},
    packages=['fungwas_stage1'],
    python_requires='>=3.8',
    install_requires=[
        'numpy>=1.20',
        'scipy',
    ],
    extras_require={
        'bgen': ['bgen-reader>=4.0'],
    },
    entry_points={
        'console_scripts': [
            'fungwas-stage1=fungwas_stage1.cli:main',
        ],
    },
)

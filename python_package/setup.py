import sys
from glob import glob
from setuptools import setup, Extension
from pybind11.setup_helpers import Pybind11Extension
from setuptools.command.build_ext import build_ext
import sys, re
import setuptools
import pybind11


# Available at setup time due to pyproject.toml
from pybind11 import get_cmake_dir
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

__version__ = "0.8"



# The main interface is through Pybind11Extension.
# * You can add cxx_std=11/14/17, and then build_ext can be removed.
# * You can set include_pybind11=false to add the include directory yourself,
#   say from a submodule.
#
# Note:
#   Sort input source files if you glob sources to ensure bit-for-bit
#   reproducible builds (https://github.com/pybind/python_example/pull/53)

# cxx_std=14

ext_modules = [
        Pybind11Extension("fasthare",
        sorted(glob("src/*.cpp")),
        include_dirs=["src/",pybind11.get_include(False),pybind11.get_include(True )],
        extra_compile_args=["-O3","-std=c++14"],  # [string]
        # Example: passing in the version to the compiled code
        define_macros = [('VERSION_INFO', __version__)],
        ),
]

setup(
    name="fasthare",
    version=__version__,
    author="Thang Dinh/Phuc Thai",
    author_email="thangdn@gmail.com",
    url="https://github.com/thang-dinh/FastHare",
    description="Reduce QUBO/Hamiltonian (exact reduction). Source code for the paper 'FastHare: Fast Hamiltonian Reduction for Large-scale Quantum Annealing, IEEE Int. Conf. on Quantum Computing and Engineering (QCE) 2022' https://arxiv.org/abs/2205.05004",
    long_description="",
    ext_modules=ext_modules,
    extras_require={"test": "pytest"},
    # Currently, build_ext only provides an optional "highest supported C++
    # level" feature, but in the future it may provide more features.
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.7",
)

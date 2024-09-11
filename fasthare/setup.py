import sys
from glob import glob
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext
import pybind11

__version__ = "1.0.4"

# Read the long description from README.md
with open("README.md", "r") as fh:
    long_description = fh.read()

ext_modules = [
    Pybind11Extension(
        "fasthare",
        sorted(glob("src/*.cpp")),
        include_dirs=["src/", pybind11.get_include()],
        extra_compile_args=["-O3"],  # Optimization flags
        cxx_std=14,  # Specify C++14 standard
        define_macros=[('VERSION_INFO', __version__)],
    ),
]

setup(
    name="fasthare",
    version=__version__,
    author="Thang Dinh/Phuc Thai",
    author_email="thangdn@gmail.com",
    url="https://github.com/thang-dinh/FastHare",
    description="Reduce QUBO/Hamiltonian (exact reduction). Source code for the paper 'FastHare: Fast Hamiltonian Reduction for Large-scale Quantum Annealing, IEEE Int. Conf. on Quantum Computing and Engineering (QCE) 2022' https://arxiv.org/abs/2205.05004",
    long_description=long_description,
    long_description_content_type="text/markdown",
    ext_modules=ext_modules,
    extras_require={"test": "pytest"},
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.10",
)

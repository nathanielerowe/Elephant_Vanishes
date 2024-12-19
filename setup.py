from setuptools import setup, Extension
from pybind11.setup_helpers import Pybind11Extension, build_ext
import numpy

def readme():
    with open('README.md') as f:
        return f.read()

ext_modules = [
    Pybind11Extension(
        '_profit',
        ['pylib/profit.cxx'],
        include_dirs=[
           # 'inc', 
           'install/include',
           'install/include/eigen3',
           '/opt/homebrew/include',
           '/opt/homebrew/Cellar/root/6.32.06_1/include/root',
           'build/_deps/sbnanaobj-src'
        ],
        library_dirs=["/opt/homebrew/opt/root/lib/root/", "build/src/", "install/lib"],
        libraries=["Core", "PROfitLib", "tinyxml2", "TreePlayer"],
        language='c++',
        cxx_std=17
    ),
]

setup(
    name='profit',
    ext_modules=ext_modules,
    include_dirs=[numpy.get_include(), "inc/"],
    version='0.1',
    description='Python library for PROfit fitting package.',
    url='https://github.com/gputnam/Elephant_Vanishes/',
    author='Gray Putnam',
    author_email='gputnam@fnal.gov',
    license='MIT',
    packages=['profit'],
    cmdclass={'build_ext': build_ext},
    install_requires=[
        'numpy'
    ],
    zip_safe=False)

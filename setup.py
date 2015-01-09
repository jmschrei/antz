from distutils.core import setup
import numpy as np

setup(
    name='antz',
    version='0.0.1',
    author='Jacob Schreiber',
    author_email='jmschreiber91@gmail.com',
    packages=['antz'],
    url='http://pypi.python.org/pypi/antz',
    license='LICENSE.txt',
    description='antz is a bioinformatics package written for python.',
    install_requires=[
        "numpy >= 1.8.0",
        "matplotlib >= 1.3.1"
    ],
)

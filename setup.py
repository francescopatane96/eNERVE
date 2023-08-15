from setuptools import setup, find_packages

setup(
    name='enerve',
    version='1.0',
    packages=find_packages(),
    install_requires=[
    'numpy == 1.23.0',
    'scikit-learn == 1.2.1',
    'biopython == 1.78',
    'epitopepredict == 0.5.0',
    'seaborn == 0.12.2',
    'pandas == 2.0.1',
    'requests == 2.31.0',
    'tensorflow == 2.13.0',
    'joblib == 1.2.0',
    'matplotlib == 3.5.0',
    'ncbi-blast+',
    'git+https://github.com/francescopatane96/tmhmm.py.git' 
    ],
)

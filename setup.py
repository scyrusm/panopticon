#from distutils.core import setup
from setuptools import setup, find_packages
import os
import io

def _read(*parts, **kwargs):
    filepath = os.path.join(os.path.dirname(__file__), *parts)
    encoding = kwargs.pop('encoding', 'utf-8')
    with io.open(filepath, encoding=encoding) as fh:
        text = fh.read()
    return text

def get_requirements(path):
    content = _read(path)
    return [
        req
        for req in content.split("\n")
        if req != '' and not req.startswith('#')
    ]

install_requires = get_requirements('requirements.txt')

setup(
    name='panopticon',
    version='0.1dev',
    author='Samuel Markson',
    author_email='SamuelC_Markson@dfci.harvard.edu',
#    packages=['frogment' ],
    install_requires=install_requires, 
    packages=find_packages(),
#    scripts=['scripts/WindowedMeanExpressionClustering', 'scripts/RnaDnaCorrespondence'],
    scripts=['scripts/panopticon'],
    license='MIT',
    description='for combined DNA/RNA analysis in a longitudinal context',
    include_package_data = True,
    long_description=open('README.md').read(),
)


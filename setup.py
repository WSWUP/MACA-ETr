import io, re
from setuptools import setup

with io.open("README.rst", "rt", encoding="utf8") as f:
    readme = f.read()

with io.open("macaetr/__init__.py", "rt", encoding="utf8") as f:
    version = re.search(r"__version__ = \'(.*?)\'", f.read()).group(1)

requires = [
    'netCDF4',
    'numpy',
    'pandas',
    'refet>=0.3.7',
    'scipy>=1.1.0',
    'xarray'
]

tests_require = ['pytest']

classifiers = [
    'License :: OSI Approved :: Apache Software License',
    'Programming Language :: Python :: 3.7',
    'Environment :: Console',
    'Development Status :: 4 - Beta',
    'Topic :: Scientific/Engineering',
    'Intended Audience :: Science/Research'
]

setup(
    name='macaetr',
    version=version,
    description='Download MACA downscaled data, calculate ASCE ETrz',
    long_description=readme,
    author='John Volk',
    author_email='john.volk@dri.edu',
    license='Apache',
    url='https://github.com/WSWUP/MACA-ETr',
    platforms=['Windows','Linux','Mac OS X'],
    classifiers=classifiers,
    packages=['macaetr'],
    install_requires=requires,
    tests_require=tests_require,
    package_data={'macaetr': ['metadata/*']},
    include_package_data=True
)

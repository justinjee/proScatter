import os
from setuptools import find_packages, setup


def extract_version():
    """
    Extracts version values from the main matplotlib __init__.py and
    returns them as a dictionary.
    """
    with open('proscatter/__init__.py') as fd:
        for line in fd.readlines():
            if (line.startswith('__version__')):
                exec(line.strip())
    return locals()["__version__"]


def get_package_data():
    baseline_images = [
        'tests/baseline_images/%s/*' % x
        for x in os.listdir('tests/baseline_images')]

    return {
        'proscatter':
        baseline_images +
        [
            "examples/*.html",
            "examples/*.txt"
        ]}

setup(
    name="proScatter",
    version=extract_version(),
    author="Justin Jee",
    author_email="justin.jee@gmail.com",
    url="https://github.com/justinjee/proScatter",
    license="MIT",
    packages=find_packages(),
    package_dir={"proscatter": "proscatter"},
    package_data=get_package_data(),
    description="for visualizing pLink data from one or more experiments",
    # run pandoc --from=markdown --to=rst --output=README.rst README.md
    long_description=open("README.rst").read(),
    # numpy is here to make installing easier... Needs to be at the
    # last position, as that's the first installed with
    # "python setup.py install"
    install_requires=["bokeh",
                      "pandas >= 0.16.0",
                      "numpy"],
    entry_points={
        'console_scripts': [
            'proscatter=proScatter:main',
            ],
        },
    classifiers=['Intended Audience :: Science/Research',
                 'Programming Language :: Python',
                 'Topic :: Scientific/Engineering :: Bio-Informatics',
                 'Topic :: Scientific/Engineering :: Visualization',
                 'Operating System :: Microsoft :: Windows',
                 'Operating System :: POSIX',
                 'Operating System :: Unix',
                 'Operating System :: MacOS',
                 'Programming Language :: Python :: 2',
                 'Programming Language :: Python :: 2.7',
                 'Programming Language :: Python :: 3',
                 'Programming Language :: Python :: 3.3'],
    zip_safe=False)

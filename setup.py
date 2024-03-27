import pathlib
from setuptools import setup, find_packages


# The directory containing this file
HERE = pathlib.Path(__file__).parent
README = (HERE / "README.md").read_text()

setup(
    name="cas-tools",
    version="0.0.1.dev36",
    description="Cell Annotation Schema tools.",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/cellannotation/cas-tools",
    author="",
    license="Apache-2.0 license",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    package_dir={'': 'src'},
    packages=find_packages(where='src', exclude=("test*",)),
    include_package_data=True,
    install_requires=[
        "anndata==0.10.3",
        "cellxgene-census",  # No version specified, latest will be used
        "openpyxl",          # No version specified, latest will be used
        "dataclasses_json",  # No version specified, latest will be used
        "pandas",            # No version specified, latest will be used
        "ruamel.yaml",       # No version specified, latest will be used
        "jsonschema==4.4.0",
        "ordered-set==4.1.0",
        "deepmerge==1.1.0",
        "cell-annotation-schema",  # No version specified, latest will be used
        "h5py>=3.10.0",
        "numpy>=1.23.5",
        "setuptools>=65.5.1"
    ],
    entry_points={
        "console_scripts": [
            "cas=cas.__main__:main",
        ]
    },
)

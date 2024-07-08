import pathlib
from setuptools import setup, find_packages


# The directory containing this file
HERE = pathlib.Path(__file__).parent
README = (HERE / "README.md").read_text()

setup(
    name="cas-tools",
    version="1.0.1",
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
        "cellxgene-census==1.10.2",
        "openpyxl==3.1.2",
        "dataclasses-json==0.6.4",
        "pandas==2.2.1",
        "ruamel.yaml==0.18.6",
        "jsonschema==4.4.0",
        "ordered-set==4.1.0",
        "deepmerge==1.1.0",
        "cell-annotation-schema==0.2b0",
        "h5py==3.10.0",
        "numpy==1.26.4",
        "marshmallow==3.21.1",
        "python-dateutil==2.9.0",
        "schema-automator==0.5.1",
        "linkml==1.7.8",
        "oaklib==0.6.4",
        "defusedxml==0.7.1",
        "rdflib==7.0.0",
        "requests==2.32.0"
    ],
    entry_points={
        "console_scripts": [
            "cas=cas.__main__:main",
        ]
    },
)

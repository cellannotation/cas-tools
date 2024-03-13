import pathlib
from setuptools import setup, find_packages


# The directory containing this file
HERE = pathlib.Path(__file__).parent
README = (HERE / "README.md").read_text()

setup(
    name="cas-tools",
    version="0.0.1.dev31",
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
    install_requires=["anndata==0.10.3", "dataclasses_json", "pandas",
                      "ruamel.yaml", "jsonschema"],
    entry_points={
        "console_scripts": [
            "cas=cas.__main__:main",
        ]
    },
)

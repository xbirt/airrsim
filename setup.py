from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="airrsim",
    version="0.1.0",
    author="Cristian Dumitrescu",
    author_email="cristid@gmail.com",
    description="A tool for simulating Adaptive Immune Receptor Repertoires and sequence read data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/xbirt/airrsim",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    python_requires=">=3.7",
    install_requires=[
        "biopython>=1.79",
        "numpy>=1.21.0",
        "requests>=2.26.0",
        "beautifulsoup4>=4.9.3",
    ],
    entry_points={
        "console_scripts": [
            "airrsim=src.main:main",
        ],
    },
)
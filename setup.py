from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="vorpy3",
    version="3.0.1",
    author="John Ericson",
    author_email="jackericson98@gmail.com",
    description="A Python package for Voronoi analysis of molecular structures",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jackericson98/vorpy",
    packages=find_packages(),
    package_data={
        'vorpy': ['data/*.pdb', 'data/*.gro', 'data/*.txt'],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib",
        "pandas",
        "numba",
        "shapely",
        "plotly",
        "sympy"
    ],
    extras_require={
        "dev": [
            "pytest",
            "pytest-cov",
            "black",
            "flake8",
            "mypy",
        ],
    },
    include_package_data=True,
    license="MIT",
    entry_points={
        'console_scripts': [
            'vorpy=vorpy.__main__:run',
        ],
    },
) 
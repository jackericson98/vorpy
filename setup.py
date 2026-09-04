from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="vorpy3",
    version="3.5.2",
    author="John Ericson",
    author_email="jackericson98@gmail.com",
    description="A Python package for Voronoi analysis of molecular structures",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jackericson98/vorpy",
    packages=find_packages(),
    package_data={
        'vorpy': ['data/*.pdb', 'data/*.gro', 'data/*.txt', 'src/GUI/Images/*.png', 'src/GUI/Images/*.ico', 'workbench/assets/*.png'],
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
    ],
    python_requires=">=3.9",
    install_requires=[
        "numpy",
        "scipy", 
        "matplotlib",
        "pandas",
        "numba==0.62.1; sys_platform == 'darwin' and platform_machine == 'x86_64'",
        "numba; sys_platform != 'darwin' or platform_machine != 'x86_64'",
        "shapely",
        "Pillow",
    ],
    extras_require={
        "gui": [
            "PySide6>=6.7",
            "pyvista>=0.44",
            "pyvistaqt>=0.11",
        ],
        "dev": [
            "pytest",
            "pytest-cov",
            "black",
            "flake8",
            "mypy",
            "hypothesis",
            "ruff>=0.6",
        ],
    },
    include_package_data=True,
    license="MIT",
    entry_points={
        'console_scripts': [
            'vorpy=vorpy.__main__:main',
        ],
    },
)

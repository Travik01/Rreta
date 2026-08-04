"""
Setup script для установки пакета.
"""

from setuptools import setup, find_packages
from pathlib import Path

# Чтение README
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

# Чтение requirements
requirements = (this_directory / "requirements.txt").read_text().splitlines()

setup(
    name="protein-ligand-dynamics",
    version="1.0.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="Production-ready pipeline for protein-ligand molecular dynamics",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/your-repo/protein-ligand-md",
    packages=find_packages(where="src") + ["scripts"],
    package_dir={"": "src", "scripts": "scripts"},
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.10",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "protein-md=scripts.full_pipeline:main",
            "protein-analyze=scripts.analyze_results:main",
            "protein-viz=scripts.visualize:main",
        ],
    },
    include_package_data=True,
    package_data={
        "": ["*.yaml", "*.json"],
    },
)

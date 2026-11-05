from pathlib import Path

from setuptools import find_packages, setup

root = Path(__file__).parent

setup(
    name="protein-ligand-dynamics",
    version="0.1.0",
    description="Pipeline for protein-ligand molecular dynamics",
    long_description=(root / "README.md").read_text(encoding="utf-8"),
    long_description_content_type="text/markdown",
    url="https://github.com/Travik01/Rreta",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.10",
    install_requires=[
        line.strip()
        for line in (root / "requirements.txt").read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    ],
)

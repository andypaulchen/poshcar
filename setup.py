from setuptools import setup, find_packages

setup(
    name="poshcar",
    version="2.1.0",
    description="Crystal structure engine for VASP POSCAR files",
    packages=find_packages(),
    install_requires=[
        "ase",
        "pyvis",
        "rdkit",
        "scipy",
        "pymatgen",
        "chgnet",
        "nglview",
        "pandas",
        "networkx"
    ],
)

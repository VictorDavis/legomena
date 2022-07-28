from setuptools import setup

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name="legomena",
    version="1.2.0",
    description="Tool for exploring types, tokens, and n-legomena relationships",
    license="MIT",
    long_description=long_description,
    long_description_content_type="text/markdown",
    py_modules=["legomena"],
    author="Victor Davis",
    author_email="vadsql@gmail.com",
    url="https://github.com/VictorDavis/legomena",
    install_requires=["mpmath", "nltk", "numpy", "pandas", "scipy"],
    python_requires=">=3.5",
)

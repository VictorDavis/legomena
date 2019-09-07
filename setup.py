from setuptools import setup

with open("README", "r") as f:
    long_description = f.read()

setup(
    name="legomena",
    version="1.0.0",
    description="Tool for exploring types, tokens, and n-legomena relationships in text.",
    license="MIT",
    long_description=long_description,
    author="Victor Davis",
    author_email="vadsql@gmail.com",
    url="http://legomena.herokuapp.com/",
    packages=["legomena"],
    install_requires=["nltk", "numpy", "pandas", "scipy"],
    scripts=["scripts/legomena"],
)

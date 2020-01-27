import setuptools

with open("README.md", "r") as fh : 
    long_description = fh.read

setuptools.setup(
    name='noisydiff',
    version= "0.0.1",
    author="Meriem Bensouda",
    author_email="meriem.bensoudakoraichi@gmail.com",
    description="A package to extract clone expansion information precisely",
    long_description= long_description,
    long_description_content_type = "text/markdown",
    url="",
    keywords = "package, T-cells, RNAseq, gDNAseq, expansion",
    packages= setuptools.find_packages(),
    classifiers= [ 
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent"

     ],
)
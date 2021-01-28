import setuptools 

with open("README.md", "r") as fh :
   long_description = fh.read()

setuptools.setup(
   name="noisets",
   version= "0.0.1",
   author= "Meriem Bensouda Koraichi ",
   author_email = "meriem.bensoudakoraichi@gmail.com",
   description = "A package assessing experimental sampling Noise and follows TCR dynamics with time",
   long_description= long_description,
   long_description_content_type = "text/markdown",
   url=" ",
   keywords= 'package numbers calculations ',
   packages= setuptools.find_packages(), 
   classifiers = [ 
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent"

     ],

  entry_points = {'console_scripts': [
            'noiset-noise=noisets.noise_learning:main',
            'noiset-detection=noisets.detection_learning:main',
            'noiset-nullgenerator=noisets.generator:main',
            ], },
 )
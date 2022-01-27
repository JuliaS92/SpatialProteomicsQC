from setuptools import setup

with open("requirements.txt") as requirements_file:
    requirements = []
    for line in requirements_file:
        requirements.append(line)

setup(name='domaps',
      version='0.1',
      description='Python library for dynamic organellar maps',
      url='https://github.com/valbrecht/SpatialProteomicsQC',
      author='Julia Schessner, Vincent Albrecht',
      author_email='schessner@biochem.mpg.de',
      license='MIT',
      packages=['domaps'],
      install_requires = requirements,
      include_package_data=True,
      zip_safe=False,
      python_requires='>=3.6',)

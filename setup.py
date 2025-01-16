from setuptools import setup

with open("requirements.txt") as requirements_file:
    requirements = []
    for line in requirements_file:
        requirements.append(line)

setup(name='domaps',
      version='1.0.1',
      description='Python library for dynamic organellar maps',
      url='https://github.com/JuliaS92/SpatialProteomicsQC',
      author='Julia Schessner, Vincent Albrecht',
      author_email='schessner@biochem.mpg.de',
      license='MIT',
      packages=['domaps', 'domaps.gui', 'domaps.network'],
      install_requires = requirements,
      package_data={'domaps': ["annotations/*", "referencedata/*", "img/*"]},
      include_package_data=True,
      zip_safe=False,
      python_requires='>=3.7',)

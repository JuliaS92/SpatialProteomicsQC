#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from setuptools import setup

setup(name='domaps',
      version='0.1',
      description='Python library for dynamic organellar maps',
      url='https://github.com/valbrecht/SpatialProteomicsQC',
      author='Julia Schessner, Vincent Albrecht',
      author_email='schessner@biochem.mpg.de',
      license='MIT',
      packages=['domaps'],
      include_package_data=True,
      zip_safe=False)


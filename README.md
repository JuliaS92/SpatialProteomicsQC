[![Python package](https://github.com/valbrecht/SpatialProteomicsQC/actions/workflows/python-package.yml/badge.svg)](https://github.com/valbrecht/SpatialProteomicsQC/actions/workflows/python-package.yml)
[![DOI](https://zenodo.org/badge/313335423.svg)](https://zenodo.org/badge/latestdoi/313335423)



# SpatialProteomicsQC
This is a quality control tool for spatial proteomics data that has been developed in the Borner lab. Primary results illustrating the functionality of the tool have been published [here](https://doi.org/10.1038/s41467-023-41000-7) and presented [here](https://www.youtube.com/watch?v=dUrOxYHJihc). The graphical user interface is available at https://domqc.bornerlab.org. Alternatively you can clone this repository and run the user interface locally or use the underlying python library in custom code. If you encounter any issues with the website or the code please let us know through the issues on this repository and if you have any unanswered questions regarding the science behind it don't hesitate to contact us directly.

Instructions for uploading data as well as sample data are included in the interface.

## Running locally

### System requirements
There are no special hardware requirements for running the QCtool. The code was tested on Linux and Windows and should also run on Mac.

### Setting up the environment
To run the app locally create a new python environment (required version >=3.8, <= 3.10). This can either be done using [anaconda](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands) or a [native virtual environment](https://docs.python.org/3/library/venv.html). Illustration for unix systems:

```
python3 -m venv path/to/env
source path/to/env/bin/activate
```

Next, clone this repository (requires git) or download it here and install the package including its dependencies (see requirements file):
```
git clone https://github.com/valbrecht/SpatialProteomicsQC

pip install ".[gui]"
```

This step takes roughly 20 minutes, depending on download speed and hardware configuration, limited by the plotly library. If any pip installations fail due to external dependencies missing (e.g. zlib requirement for Pillow installation from source on python 3.6, or Microsoft Visual C++ 14 on windows) try installing that library using a binary file, thereby avoiding the dependency:
```
pip install --only-binary :all: Pillow
```

### Running from jupyter notebook
To run the tool on your computer and to be able to interact with the data and the interface through pyhton code use jupyter. This is not included in the requirements so you need to add it to your environment if you haven't used anaconda. You can then run the jupyter notebook in the webapp folder to run a locally hosted interface. The only time you *need* to do this is when working with files larger than 80 MB, since they can't be loaded through the interface directly. Instead adjust one of the two code cells at the end of the notebook, to load the file directly from your file system.

### Running a (local) server
To access the GUI without directly interacting with the code, or to host your own instance of the tool, you can run a [panel serve](https://panel.holoviz.org/user_guide/Deploy_and_Export.html#launching-a-server-on-the-commandline) command from the command line:
```
cd domaps/webapp
panel serve QCtool.py
```


### Code Contributions
Generally code contributions are welcome - if you have best practices or benchmarks that you want to see represented here get in touch, or open a PR.

For releases we follow the alphaX ecosystem https://github.com/MannLabs/alphashared/blob/main/.github/workflows/README.md.

Requirements are split between what is required for data analysis and what is required to run the panel based GUI. Please don't add code for data processing that requires any of the GUI libraries. All newly added code should be tested by doctests or pytest and documented with docstrings. For major additions please include a jupyter notebook.

This is a collection of routines useful for analysis of MRI data.
=================================================================

Code is organized in a standard python package layout as `discussed here <http://guide.python-distribute.org/creation.html>`_.

    The root directory contains the SARlabpy package directory, the doc directory (which largely is populated automatically from docstrings), the README.rst file and some other minor files.

        The SARlabpy directory contains:

          * __init__.py which is needed correct reading of the package content
          * **io** -- input output routines to read BRUKER data and other things
          * **test** a directory containing unit tests for code in all the other directories
          * **fmoosvi** Firas' sandbox
          * **DCE** DCE fitting and AIF simulations - mostly Tammo's work

The main objectives for this code is

  * Reading, and reconstructing BRUKER data
  * Analyzing DCEMRI data
  * Data presentation
  * Integration and analysis of immunohistochemistry data


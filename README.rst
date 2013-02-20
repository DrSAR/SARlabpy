=====================================================================
SARlabpy -- A collection of routines used for analysis of MRI data.
=====================================================================
The code is organized in a standard python package layout as `discussed here <http://guide.python-distribute.org/creation.html>`_. More handy tips on organizing this can be had from `infinite monkeys <http://infinitemonkeycorps.net/docs/pph/>`_.

    The root directory contains 
    
    * the doc directory (which contains files used for the automated documenation creation), 
    * the html directory: this directory contains the orphan branch hg-pages which is pushed to the
      documentation that's hosted on `code.sarlab.ca <http://code.sarlab.ca>`_
    * the SARlabpy directory contains:

          * __init__.py which is used to initialize the packageo (see `here for more information <http://docs.python.org/2/tutorial/modules.html#packages>`_)
          * **test** a directory containing unit tests for code in all the other directories
          * **io** -- input output routines to read BRUKER data and other things
          * **fmoosvi** Firas' sandbox
          * **DCE** DCE fitting and AIF simulations - mostly Tammo's work

The main objectives for this code is

  * Reading, and reconstructing BRUKER data
  * Analyzing DCEMRI data
  * Data presentation
  * Integration and analysis of immunohistochemistry data


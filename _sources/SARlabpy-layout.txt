Package Organization for SARlabpy
-----------------------------------
The code is organized in a standard python package layout as `discussed here <http://guide.python-distribute.org/creation.html>`_. More handy tips on organizing this can be had from `infinite monkeys <http://infinitemonkeycorps.net/docs/pph/>`_.

The root directory contains 
    
* *doc*: directory of files used for the automated documenation creation, 
* *html*: directory that contains the orphan branch hg-pages which is pushed to the
      documentation that's hosted on `code.sarlab.ca <http://code.sarlab.ca>`_
* *SARlabpy* contains:

   * __init__.py which is used to initialize the packageo (see `here for more information <http://docs.python.org/2/tutorial/modules.html#packages>`_)
   * **test** a directory containing unit tests for code in all the other directories
   * **io** -- input output routines to read BRUKER data and other things
   * **fmoosvi** Firas' sandbox
   * **DCE** DCE fitting and AIF simulations - mostly Tammo's work



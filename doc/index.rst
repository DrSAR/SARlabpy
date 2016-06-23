.. sarpy documentation master file, created by
   sphinx-quickstart on Fri Feb 15 16:52:01 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SARpy's documentation!
==========================

The code is hosted in the `private repository by github <https://github.com/DrSAR/SARlabpy.git>`_ and `by gitlab <https://gitlab.com/DrSAR/SARlabpy>`_.
If you are stuck with branches and commits, try this quickstart: (assuming you have no uncommitted local changes or local commits not pushed to the central repo that you care about):

.. code-block:: bash

   rm -rf SARpy # remove any traces of the local copy of the repository
   git clone git@pfeifer.phas.ubc.ca:SARlabpy # make a new fresh clone
   cd SARlabpy # you will now be on the master branch
   git checkout --track origin/thatbranchthatIusedbefore # if you are continuing on an old existing branch
   git checkout -b newfeaturebranch # if you are starting on something new

.. include:: ../README.rst

Howto, Guides, Examples
----------------------------------

.. toctree::
   :maxdepth: 2

   sarpy-layout
   Workflow
   Documenting-Testing
   BRUKERIOexamples
   SARloggerexamples

Doc's harvested from source code
----------------------------------

.. toctree::
   :maxdepth: 2

   modules

Indices and tables
-------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


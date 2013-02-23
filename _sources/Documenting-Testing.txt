Documenting and Testing SARlabpy
================================

*Wherein I say unto you: verily document and test thy code*

To make this code usable, it needs to be documented. The easiest and at the same time most helpful way to write documentation are a few examples (copied verbatim from your interactive session) into the docstring. The amazing thing is that this docstring can be parsed by the doctest module in order to **automatically** run those examples and compare them to the expected output.

We should do more of that!

Inner Workings
---------------
The documentation in this package makes use of Sphinx, the automated documentation project. Essentially it harvests the docstrings of all modules and puts it together in a nice and appetizing way. The instruction for this automated compilation of HTML documentation is in the **doc** folder (Makefile, config.py and some static pre-written, global document pages). Simply going to that directory and issuing *make html* will recreate the documentation and store it in the **html** folder. The **html** folder is bizarrely the working copy of a different, special (orphan) branch named gh-pages. When we push the html pages to github they are processed and will be hosted on the `github pages http://DrSAR.github.com/SARlabpy <http://DrSAR.github.com/SARlabpy>`_. That domain is also pointed to from our own domain `code.sarlab.ca <http://code.sarlab.ca>`_ 

Working with Sphinx -- If you really must
-----------------------------------------
If you want to work with the documentation (outside of the in-file documentation of your code) do the following in your local repository root:

.. code-block:: bash

  mkdir html #
  cd html
  git clone git@pfeifer.phas.ubc.ca:SARlabpy .
  git checkout --track origin/gh-pages
  rm -rf * # yes really, delete everything in html/

You have now created a directory that contains a repository (in a repository so to say). The orphan branch (i.e. not based off of anything) is checkout out. Note how you change branches simply by entereing and leaving that directory. (Not sure how this works on a gui client...). Now recreate the documentation:

.. code-block:: bash

   # modify the ReSTructured text files in doc/
   # modify the documentation in the source code
   cd ../doc # this is the directory where the doc setup script is
   make html # this will recompile the html files, store them in html/ 
   cd ../html # go back where the html files are
   git add newfiles # if any
   git commit -a -m 'documentation updated'

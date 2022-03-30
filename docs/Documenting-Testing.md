Documenting and Testing sarpy
================================

To make this code usable, it needs to be documented. The easiest and at the same time most helpful way to write documentation are a few examples (copied verbatim from your interactive session) into the docstring. The amazing thing is that this docstring can be parsed by the doctest module in order to **automatically** run those examples and compare them to the expected output.

We should do more of that!

Inner Workings
---------------
The documentation in this package makes use of Sphinx, the automated documentation project. Essentially it harvests the
docstrings of all modules and puts it together in a nice and appetizing way. The instruction for this automated
compilation of HTML documentation is in the `doc` folder (`Makefile`, `config.py` and some static pre-written, global
document pages). Simply going to that directory and issuing `make html` will recreate the documentation and store it in
the `html` folder. The `html` folder is bizarrely the working copy of a different, special (orphan) branch named
gh-pages. When we push the html pages to github they are processed and will be hosted on (github pages
https://DrSAR.github.com/SARlabpy)[https://DrSAR.github.com/SARlabpy]. That domain is also pointed to from our own
domain https://code.sarlab.ca. 

`make html` will write a lot of stuff into _build/html. 
One way to publish this by hand is 
 * to checkout gh-pages,
 * to wipe the entire tree (without getting rid of the `.git` directory and `CNAME`),
 * to copy the content of html (temporarily stashed somewhere?) into the project root directory, and, finally, 
 * to (force-?)push to gh-pages upstream.


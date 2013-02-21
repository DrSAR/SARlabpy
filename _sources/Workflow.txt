Suggested Workflow for the Code Repository
===========================================

*Wherein I lay down the law on how to git it done.*

There are a great many good suggested workflows, e.g. the `Successful git branching model <http://nvie.com/posts/a-successful-git-branching-model/>`_. A slightly opposing view is discussed by Randy Fay where he favours `a Rebase workflow <http://randyfay.com/content/rebase-workflow-git>`_. One might want to look at both article and try to understand them.

I suggest we do the following:

* Branch master is a branch that should in most situations have (almost) fully functional code. This branch can have tags attached so that you can go back and reproduce your work.
* Any given feature, bugfix, enhancement should be worked on in a (local) branch. Feel free to push them to the central repo (as that branch) but nobody has any right to expect to find workable code in there.

  .. code-block:: bash

        git checkout master
        git branch mynewfeature
        # do work ...
        git commit # maybe even a few times
        git push origin mynewfeature  # if you want a backup

  * If time comes to make the feature, bugfix, whatever available: 

  .. code-block:: bash

        git checkout master
        git merge --no-ff -m 'feature solves #33' mynewfeature
        git push

* This is more or less following the merge workflow. According to Randy Fay, this might get messy. If it becomes a problem we might switch to the rebase workflow.

* As a convention for the use of branches I suggest:

  * *master* the slighlty more sacrosanct branch
  * *hg-pages* do not touch! this is a special orphan branch that is used for auto-documentation
  * *firas-xxx* replace firas with your name and xxx if you find that there is a feature that you are working on at the moment
  * *xxx* branch dedicated to feature xxx

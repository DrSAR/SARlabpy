Suggested Workflow for the Code Repository
===========================================

*Wherein I lay down the law on how to git it done.*

There are a great many good suggested workflows, e.g. the `Successful git branching model <http://nvie.com/posts/a-successful-git-branching-model/>`_. A slightly opposing view is discussed by Randy Fay where he favours `a Rebase workflow <http://randyfay.com/content/rebase-workflow-git>`_. One might want to look at both article and try to understand them.

I have come to believe in the rebase workflow and hence I suggest we do the following:

* Branch master is a branch that should in most situations have (almost) fully functional code. This branch can have tags attached so that you can go back and reproduce your work. -- Currently only the central integrator (SAR) can push to it.
* Any given feature, bugfix, enhancement should be worked on in a (local) branch. Feel free to push them to the central repo (as that branch) but nobody has any right to expect to find workable code in there.

  .. code-block:: bash

        git checkout master
        git pull                 # this brings your local master in line with what's on origin
        git checkout -b mynewfeature # same as: git branch mynewfeature ; git checkout mynewfeature
        # do work ...
        git commit # maybe even a few times
        git push origin mynewfeature  # if you want a backup

* The objective here is that your code will at somepoint work with the rest of the code (which, of course, is also evolving). Hence to effect that merge at some future date, you need to stay on top of things. You do this by 'rebasing' all your as yet unmerged changes on the most recent version of master.

  .. code-block:: bash

        git fetch
        git rebase origin/master # you do this while still checked out on mynewfeature

* If time comes to make the feature, bugfix, whatever available, **tell me**. Make sure your branch is (a) rebased on top of recent master (b) pushed to the central repo so I can get to it. If you have pushed before then a simple `git push origin mynewfeature` will fail ("Updates were rejected because the tip of your current branch is behind its remote counterpart") You absolutely have to resist to follow git's usually helpful advice in merging the remote changes, pulling or whatever. What is required is a **force push** (a feature which is absent from a few GUI frontends to git! Even though it is quite safe if used on branches that are your own alone):
        
  .. code-block:: bash

        git push -f origin mynewfeature
 
And this is what I will do:

  .. code-block:: bash

        git checkout master
        git merge --no-ff mynewfeature
        git push


* As a convention for the use of branches I suggest:

  * *master* the slighlty more sacrosanct branch
  * *hg-pages* do not touch! this is a special orphan branch that is used for auto-documentation
  * *firas/xxx* replace firas with your name or initials and xxx if you find that there is a feature that you are working on at the moment. A branch that is labelled with a name is totally off-limits for anyone else to force push to! We have log files. And we will find you.

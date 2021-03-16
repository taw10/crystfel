Contributing to CrystFEL
========================

Citing CrystFEL
---------------

The easiest way to ensure the future of CrystFEL is to use it and cite it in your publications!  If your paper doesn't already appear on the [citations page](https://www.desy.de/~twhite/crystfel/citations.html), or if the information there is incomplete, be sure to get in touch about it.


Reporting issues
----------------

If you notice a problem, don't suffer in silence!  There are three ways to report a problem:

* Send an email to me (taw@physics.org).
* Open an issue on the [DESY GitLab](https://gitlab.desy.de/thomas.white/crystfel), if you have a DESY or [Helmholtz AAI](https://login.helmholtz.de/home/) account.
* Open an issue on [Github](https://github.com/taw10/crystfel).  You will need a (free) Github account for this.

Opening an issue on one of the project pages is the best way if your problem is well defined and reproducible, for example "I clicked on a button and the program crashed".

When the problem is less well defined, for example "I upgraded to the latest version and now my data doesn't come out as well as before", it's better to get in touch directly because I'll probably need to look at your data to find the problem.  In the case of a change in final data quality, I'd rather hear about it directly than read about it an a journal article, which is rarely useful for making further improvements (even if the change is positive!).

The CrystFEL project on the [DESY Jira](https://agira.desy.de/projects/CRYS/issues) server is now regarded as deprecated.  Please do not open any new issues there.


Contributing patches
--------------------

Patches to fix issues are always welcome.  Use one of the following methods, in rough order of preference:

* Create a merge request on the [DESY GitLab](https://gitlab.desy.de/thomas.white/crystfel), if you have a DESY or [Helmholtz AAI](https://login.helmholtz.de/home/) account.
* Open a pull request on [Github](https://github.com/taw10/crystfel).  You will need a (free) Github account for this.
* Send a patch to me (taw@physics.org).

The CrystFEL repository on the [DESY Bitbucket](https://stash.desy.de/projects/CRYS/repos/crystfel/) server is now regarded as deprecated.  Please do not create any new pull requests there.

For larger changes, such as adding a new feature or extending the current functionality to a new type of data, it's best to get in touch before starting more than a proof of concept.  It's very important that new features are compatible with all the existing methods, and it can be difficult to achieve full compatability without a deep knowledge of the internals of CrystFEL. Hundreds of scientific results depend on us getting this right.  Bigger changes to CrystFEL will take much longer to review and merge, and we won't rush it.  Not even if you have a paper deadline coming up.  Feel free to fork the project in the meantime.

All contributions must be under the [GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html).  The copyright holder of your contributed code depends on your situation, and we do not require any re-assignment of copyright when contributing.  In Germany, the copyright holder will most likely be your employer.  Please update the copyright and authorship details at the top of each file you make changes to.  Including your email address in the source code is optional.

Finally, please make sure that your branch or patch is "clean".  That means:

* Follow the formatting style of the existing code rigorously.
* Remove *all* trailing whitespace.
* Filter out irrelevant changes (*always* use `git add -p`).
* Rebase and squash commits to (as far as possible) form a series of self-contained changes.


Translations
------------

Are you interested in translating CrystFEL into another language?  Changes to CrystFEL to enable easy translation are planned.  It would be helpful to know how much demand there is for this, as well as which languages might be added.  Please get in touch if you're interested.

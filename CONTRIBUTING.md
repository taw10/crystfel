Contributing to CrystFEL
========================

Citing CrystFEL
---------------

The easiest way to ensure the future of CrystFEL is to use it and cite it in your publications.  If your paper doesn't already appear on the [citations page](https://www.desy.de/~twhite/crystfel/citations.html), or if the information there is incomplete, be sure to get in touch about it.  Citing any of the following papers should be enough to trigger your paper's inclusion in the list:

* T. A. White. "Processing serial crystallography data with CrystFEL: a step-by-step guide". Acta Cryst. D75 (2019).  [doi:10.1107/S205979831801238X](https://doi.org/10.1107/S205979831801238X)
* T. A. White, V. Mariani, W. Brehm, O. Yefanov, A. Barty, K. R. Beyerlein, F. Chervinskii, L. Galli, C. Gati, T. Nakane, A. Tolstikova, K. Yamashita, C. H. Yoon, K. Diederichs and H. N. Chapman. "Recent developments in CrystFEL". J. Applied Crystallography 49 (2016) p680-689.  [doi:10.1107/S1600576716004751](http://dx.doi.org/10.1107/S1600576716004751)
* T. A. White, A. Barty, F. Stellato, J. M. Holton, R. A. Kirian, N. A. Zatsepin and H. N. Chapman. "Crystallographic data processing for free-electron laser sources". Acta Cryst. D69 (2013), p1231–1240.  [doi:10.1107/S0907444913013620](http://dx.doi.org/10.1107/S0907444913013620)
* T. A. White, R. A. Kirian, A. V. Martin, A. Aquila, K. Nass, A. Barty and H. N. Chapman. "CrystFEL: a software suite for snapshot serial crystallography". J. Appl. Cryst. 45 (2012), p335–341.  [doi:10.1107/S0021889812002312](http://dx.doi.org/10.1107/S0021889812002312)

You can also cite CrystFEL directly using one of the following DOIs:

* All versions: [10.5281/zenodo.8183384](https://zenodo.org/record/8183384).
* Version 0.10.2: [10.5281/zenodo.8183385](https://zenodo.org/record/8183385).


Getting started
---------------

Want to get involved with the CrystFEL project, but don't know where to start?  Look for [issues tagged as "good first issue"](https://gitlab.desy.de/thomas.white/crystfel/-/issues?label_name%5B%5D=good+first+issue) or ["help wanted"](https://gitlab.desy.de/thomas.white/crystfel/-/issues?label_name%5B%5D=help+wanted).  These issues are not necessarily easy, but they don't involve invasive changes to CrystFEL's guts.  In fact, some of them don't involve programming at all.

Perhaps you're interested in translating CrystFEL into another language?  Changes to CrystFEL to enable easy translation are planned.  It would be helpful to know how much demand there is for this, as well as which languages might be added.  Please get in touch if you're interested.


Reporting issues
----------------

If you notice a problem, don't suffer in silence!  There are three ways to report a problem:

* Send an email to me (taw@physics.org).
* Open an issue on the [DESY GitLab](https://gitlab.desy.de/thomas.white/crystfel), if you have a DESY or [Helmholtz AAI](https://login.helmholtz.de/home/) account.
* Open an issue on [Github](https://github.com/taw10/crystfel).  You will need a (free) Github account for this.

Opening an issue on one of the project pages is the best way if your problem is well defined and reproducible, for example "I clicked on a button and the program crashed".

When the problem is less well defined, for example "I upgraded to the latest version and now my data doesn't come out as well as before", it's better to get in touch directly because I'll probably need to look at your data to find the problem.  In the case of a change in final data quality, I'd rather hear about it directly than read about it an a journal article, which is rarely useful for making further improvements (even if the change is positive!).

The CrystFEL project on the [DESY Jira](https://agira.desy.de/projects/CRYS/issues) server is now regarded as deprecated.  Please do not open any new issues there.


"Getting your new paper into CrystFEL"
--------------------------------------

You're working on an exciting new analysis method for serial crystallography, and you want to implement your new method in CrystFEL?  Great!  Be aware, however, that new methods often require deeper changes within CrystFEL.  It's very important that new methods compose well with, or at least are compatible with, all the existing features.  For example, your new method will need generalise to all crystal symmetries and all types of detector (especially with multiple panels or tilted panels), and use the same data formats.  Hundreds of scientific results depend on us getting this right, so we won't rush it - not even if you have a paper deadline coming up.

The best way to proceed is to base your paper on a proof of concept, probably as a completely separate program.  For preference, use CrystFEL's data formats rather than inventing your own (if you're programming in C/C++, you can even directly use `libcrystfel`).  This makes the publication of your paper independent of how long it takes to implement in mainline CrystFEL.  Afterwards, we can look together at how it should be implemented in CrystFEL.

Feel free to get in touch for advice and support at all stages!


Contributing patches
--------------------

Use one of the following methods, in rough order of preference:

* Create a merge request on the [DESY GitLab](https://gitlab.desy.de/thomas.white/crystfel), if you have a DESY or [Helmholtz AAI](https://login.helmholtz.de/home/) account.
* Open a pull request on [Github](https://github.com/taw10/crystfel).  You will need a (free) Github account for this.
* Send a patch to me (taw@physics.org).

The CrystFEL repository on the [DESY Bitbucket](https://stash.desy.de/projects/CRYS/repos/crystfel/) server is now regarded as deprecated.  Please do not create any new pull requests there.

All contributions must be under the [GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html), with option to re-distribute under a later GPL version.  The copyright holder of your contributed code depends on your situation, and we do not require any re-assignment of copyright when contributing.  In Germany, the copyright holder will probably be your employer.  Please update the copyright and authorship details at the top of each file you make changes to.  Including your email address in the source code is optional.

Finally, please make sure that your branch or patch is "clean".  That means:

* Follow the formatting style of the existing code rigorously.
* Remove *all* trailing whitespace.
* Filter out irrelevant changes (*always* use `git add -p`).
* Rebase and squash commits to (as far as possible) form a series of self-contained changes.

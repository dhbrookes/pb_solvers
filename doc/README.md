# Sphinx tool for generating github pages

This is the directory for generating gh pages
through python.

## Installing necessary packages

The steps that I had to go through to get it working are:

1. pip install sphinx
2. pip install -e git+git://github.com/michaeljones/sphinx-to-github.git#egg=sphinx-to-github
3. pip install sphinx_rtd_theme

## Setup of gh-pages branch for running with sphinx

1. To set up, follow the instructions 
[here](http://lucasbardella.com/blog/2010/02/hosting-your-sphinx-docs-in-github)
2. Change in the Makefile on the master branch directory `docs`, the keyword
BUILDDIR to the location of the gh-pages branch of the repo from step 1.

## Editing the website

1. Make changes in the master branch of the rst files, 
located in `doc/source/`
2. `cd ../`
3. make html - this builds in the gh-pages repo, located at
BUILDDIR.
4. Go to the gh-pages repo directory.
5. `git add * && git commit -m "..."`
6. If `git push origin gh-pages` fails, then
  5. `git fetch origin gh-pages`  
  6. `git merge --strategy=ours origin/gh-pages`
  7. `git push origin gh-pages`

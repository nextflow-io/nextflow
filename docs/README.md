# Nextflow Documentation 

Nextflow documentation is written using [Sphinx](http://www.sphinx-doc.org/) which 
uses the [reStructuredText](https://en.wikipedia.org/wiki/ReStructuredText) file format.

A quick intro to reStructuredText is available at [this link](http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html).

To edit and contribute to the documentation you only need a text editor to change the
appropriate `.rst` files in this directory.

Once you have edited the documentation files verify the your changes are correctly applied
using the command below to generate the HTML files:

```
make html
```


Then open the `_build/html/index.html` file with your browser and navigate the documentation
you have modified.


### Dependencies

Sphinx can be installed either with

```
pip install -U Sphinx
```

or

```
conda install sphinx
```

### Theme 

Docs uses the [sphinx_rtd_theme](https://github.com/readthedocs/sphinx_rtd_theme) theme. 

To update it, clone the bove repo, then copy `sphinx_rtd_theme/sphinx_rtd_theme` directory 
into `docs/_themes`.  

### License

Nextflow documentation is distributed under 
[Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0) license](https://creativecommons.org/licenses/by-sa/4.0/).

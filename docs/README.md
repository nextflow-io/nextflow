# Nextflow Documentation

Nextflow documentation is written using [Sphinx](http://www.sphinx-doc.org/) which 
uses the [reStructuredText](https://en.wikipedia.org/wiki/ReStructuredText) file format.

A quick intro to reStructuredText is available at [this link](http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html).

To edit and contribute to the documentation you only need a text editor to change the
appropriate `.rst` files in this directory.

Once you have edited the documentation files verify the your changes are correctly applied
using the command below to generate the HTML files:

```
make clean html
```


Then open the `_build/html/index.html` file with your browser and navigate the documentation
you have modified.


### Dependencies

Sphinx can be installed either with

```
pip install -U sphinx==3.5.4
```

or

```
conda install sphinx==3.5.4
```

### Theme

Nextflow documentation uses the [Read The Docs theme for Sphinx](https://github.com/readthedocs/sphinx_rtd_theme) theme.

In order to build the document install the Read The Docs theme using the command:

```
pip install sphinx_rtd_theme
```

### License

Nextflow documentation is distributed under
[Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0) license](https://creativecommons.org/licenses/by-sa/4.0/).

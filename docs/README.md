# Nextflow Documentation 

Nextflow documentation is written using [Sphinx](http://www.sphinx-doc.org/) which 
uses the [reStructuredText](https://en.wikipedia.org/wiki/ReStructuredText) file format.

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

# Nextflow Documentation

Nextflow documentation is written using [Sphinx](http://www.sphinx-doc.org/), [MyST](https://myst-parser.readthedocs.io/en/latest/) which is an extended version of Markdown for Sphinx, and the [Read The Docs theme for Sphinx](https://github.com/readthedocs/sphinx_rtd_theme).


## Dependencies

The most convenient approach is to create a Conda environment with Python 3.7 (other versions may work but haven't been tested).

The build dependencies can be installed with `pip`:

```bash
pip install sphinx==5.3.0 myst-parser==0.18.1 sphinx-rtd-theme==1.1.1
```


## Contributing

To edit and contribute to the documentation, you only need a text editor to change the appropriate `.md` files in this directory.

Once you have made your changes, run the following command to build the HTML files:

```bash
make clean html
```

Then open the `_build/html/index.html` file in your browser and verify your changes.


## License

Nextflow documentation is distributed under
[Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0) license](https://creativecommons.org/licenses/by-sa/4.0/).

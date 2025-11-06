# Nextflow Documentation

Nextflow documentation is written using [Sphinx](http://www.sphinx-doc.org/), [MyST](https://myst-parser.readthedocs.io/en/latest/) which is an extended version of Markdown for Sphinx, and the [Read The Docs theme for Sphinx](https://github.com/readthedocs/sphinx_rtd_theme).


## Dependencies

The most convenient approach is to create a Conda environment with Python 3.7 (other versions may work but haven't been tested).

The build dependencies can be installed with `pip`:

```bash
cd docs
pip install -r requirements.txt
```

Alternatively, you can use the Dockerfile to build the docs in a container (see below).


## Contributing

To edit and contribute to the documentation, you only need a text editor to change the appropriate `.md` files in this directory.

Once you have made your changes, run the following command to build the HTML files:

```bash
make clean html
```

Alternatively, you can use the Dockerfile to build the docs in a container:

```bash
docker build -t nextflow/sphinx:5.3.0 .
docker run -v $(pwd):/tmp nextflow/sphinx:5.3.0 -- make html
```

Then start up a local http server and open `localhost:8080` in your browser to verify the changes:

```bash
python -m http.server 8080 --directory _build/html/
```


## License

Nextflow documentation is distributed under
[Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0) license](https://creativecommons.org/licenses/by-sa/4.0/).

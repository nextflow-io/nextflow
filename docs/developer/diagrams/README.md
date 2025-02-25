# Class Diagrams

This directory contains class diagrams of the Nextflow source code, abridged and annotated for relevance and ease of use.

Each node is a class. Fields are selectively documented in order to show only core data structures and the classes that "own" them. Methods are not explicitly documented, but they are mentioned in certain links where appropriate.

Links between classes denote one of the following relationships:

- Inheritance (`A <|-- B`): `B` is a subclass of `A`
- Composition (`A --* B`): `A` contains `B`
- Instantiation (`A --> B : f`): `A` creates instance(s) of `B` at runtime via `A::f()`

Some links are commented out or not included at all, in order to focus on the most important classes and relationships. You can view these "hidden" links by simply uncommenting them, but I have found that their significance is sufficiently clear from the description files.

A separate diagram description is provided for each package. These files are interoperable, which means that you can combine any subset of files into a larger diagram description. The `merge-diagrams.sh` can create a merged file for you automatically, and it includes a sensible default set of packages.

You can use the [Mermaid Live Editor](https://mermaid-js.github.io/mermaid-live-editor/edit) or the [Mermaid CLI](https://github.com/mermaid-js/mermaid-cli) to render the diagram in a variety of image formats.

# AGENTS.md

Nextflow is a scientific workflow management system written in Groovy and Java, with core modules in `modules/` and plugins in `plugins/`.

## Build and test

- `make compile` builds and writes `.launch.classpath`; `./launch.sh run script.nf` runs that development build.
- Before committing, run the tests of each module you changed: `make test module=nf-lang`. For a single spec, `./gradlew test --tests "SomeSpec"`.

## Commits

- New source files carry the Apache 2.0 license header.
- Sign off every commit (`git commit -s`). Nextflow PRs require DCO sign-off.
- `Assisted-by: <tool and version>` may be used for AI tooling. `Co-authored-by:` should only be used for human co-authors.

## Further reading

- [Testing](.agents/docs/testing.md): S3 integration tests, integration workflows, cloud validation.
- [Release](.agents/docs/release.md): cutting a Nextflow or plugin release.
- `docs/adr/`: architecture decision records; check the relevant one before changing a design it covers.

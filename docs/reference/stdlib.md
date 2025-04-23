(stdlib-page)=

# Standard library

This section describes the Nextflow standard library. The standard library consists of built-in constants, functions, and types.

```{toctree}
:maxdepth: 1

stdlib-constants
stdlib-functions
stdlib-types
```

## Groovy and Java classes

Any Groovy or Java class that is available to Nextflow at runtime can be used in a Nextflow script. The following classes are imported by default:

- `groovy.lang.*`
- `groovy.util.*`
- `java.io.*`
- `java.lang.*`
- `java.math.BigDecimal`
- `java.math.BigInteger`
- `java.net.*`
- `java.util.*`

All other classes must be referenced by their fully-qualified name:

```nextflow
def vals = [1, 2, 3]
println groovy.json.JsonOutput.toJson(vals)
```

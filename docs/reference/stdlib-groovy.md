(stdlib-groovy)=

# Groovy and Java classes

This page describes...

## Default classes

`groovy.lang.*`
: Description

`groovy.util.*`
: Description

`java.io.*`
: Description

`java.lang.*`
: Description

`java.math.BigDecimal`
: Description

`java.math.BigInteger`
: Description

`java.net.*`
: Description

`java.util.*`
: Description

## Other classes

All other classes must be referenced by their fully-qualified name. For example:

```nextflow
def vals = [1, 2, 3]
println groovy.json.JsonOutput.toJson(vals)
```

:::{note}
The set of classes in Nextflow's runtime classpath can change between different Nextflow versions. As a best practice, any code that uses classes outside the Nextflow standard library should either (1) be refactored to only use the Nextflow standard library or (2) be refactored as a {ref}`plugin <plugins-dev-page>` with explicit dependencies.
:::
(stdlib-groovy)=

# Groovy and Java classes

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

:::{note}
The set of classes in Nextflow's runtime classpath can change between different Nextflow versions. As a best practice, any code that uses classes outside the Nextflow standard library should either be refactored to only use the Nextflow standard library or be refactored as a {ref}`plugin <dev-plugins-page>` with explicit dependencies.
:::


# `nextflow.ast`

The `nextflow.ast` package implements the Nextflow language extensions as AST transforms.

## Class Diagram

```{mermaid} diagrams/nextflow.ast.mmd
```

```{note}
Some classes may be excluded from the above diagram for brevity.
```

## Notes

The Nextflow scripting language is essentially Groovy with some extensions, implemented as transformations to the abstract syntax tree (AST). Every Nextflow script is syntactically (but not semantically) valid Groovy.

You can see the effect of Nextflow's AST transforms by using the Nextflow console:

1. Run `nextflow console` to open the console
2. Enter a Nextflow script
3. Execute the script
4. Go to **Script** > **Inspect AST**

Here is the example from {ref}`getstarted-first`:

```groovy
params.str = 'Hello world!'

process splitLetters {
  output:
    path 'chunk_*'

  """
  printf '${params.str}' | split -b 6 - chunk_
  """
}

process convertToUpper {
  input:
    path x
  output:
    stdout

  """
  cat $x | tr '[a-z]' '[A-Z]'
  """
}

workflow {
  splitLetters | flatten | convertToUpper | view { it.trim() }
}
```

Here it is after being parsed and de-sugared by Groovy:

```groovy
params.str = 'Hello world!'

process( splitLetters( {
  output:
    path('chunk_*')

  """
  printf '${params.str}' | split -b 6 - chunk_
  """
} ))

process( convertToUpper( {
  input:
    path(x)
  output:
    stdout

  """
  cat $x | tr '[a-z]' '[A-Z]'
  """
} ))

workflow({
  splitLetters | flatten | convertToUpper | view { it.trim() }
})
```

Here it is after being transformed by Nextflow (whitespace edited for readability):

```groovy
import static nextflow.Nextflow.*

import org.apache.commons.lang.StringUtils as StringUtils
import groovy.transform.Field as Field
import java.nio.file.Path as Path
import nextflow.Channel as Channel
import nextflow.util.Duration as Duration
import nextflow.util.MemoryUnit as MemoryUnit
import nextflow.io.ValueObject as ValueObject
import nextflow.Channel as channel

@groovy.transform.BaseScript
public class script1677225313239 extends nextflow.script.BaseScript { 

    public script1677225313239() {
        nextflow.script.ScriptMeta.get(this).setDsl1ProcessNames(['splitLetters', 'convertToUpper'])
    }

    public script1677225313239(final groovy.lang.Binding context) {
        super.setBinding(context)
        nextflow.script.ScriptMeta.get(this).setDsl1ProcessNames(['splitLetters', 'convertToUpper'])
    }

    public static void main(final java.lang.String[] args) {
        org.codehaus.groovy.runtime.InvokerHelper.runScript(script1677225313239, args)
    }

    @groovy.transform.Generated
    protected java.lang.Object runScript() {
        params.str = 'Hello world!'

        this.process('splitLetters', { 
            this._out_path('chunk_*')
            new nextflow.script.BodyDef(
                {
                    "printf '$params.str' | split -b 6 - chunk_"
                },
                '"""\n  printf \'${params.str}\' | split -b 6 - chunk_\n  """\n',
                'script',
                [
                    new nextflow.script.TokenValRef('params.str', 8, 13)
                ]
            )
        })

        this.process('convertToUpper', { 
            this._in_path(new nextflow.script.TokenVar('x'))
            this._out_stdout()
            new nextflow.script.BodyDef(
                { 
                    "cat $x | tr '[a-z]' '[A-Z]'"
                },
                '"""\n  cat $x | tr \'[a-z]\' \'[A-Z]\'\n  """\n',
                'script',
                [
                    new nextflow.script.TokenValRef('x', 19, 8)
                ]
            )
        })

        this.workflow({ 
            new nextflow.script.BodyDef(
                {
                    splitLetters | flatten | convertToUpper | this.view({ 
                        it.trim()
                    })
                },
                '  splitLetters | flatten | convertToUpper | view { it.trim() }\n',
                'workflow',
                [
                    new nextflow.script.TokenValRef('flatten', 24, 18),
                    new nextflow.script.TokenValRef('splitLetters', 24, 3),
                    new nextflow.script.TokenValRef('convertToUpper', 24, 28)
                ]
            )
        })
    }

}
```

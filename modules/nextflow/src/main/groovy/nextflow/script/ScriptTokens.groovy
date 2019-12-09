/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package nextflow.script

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import groovy.transform.TupleConstructor
/**
 * Presents a variable definition in the script context.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ToString
@EqualsAndHashCode
@TupleConstructor
class TokenVar {

    /** The variable name */
    String name

}

/**
 *  A token used by the DSL to identify a 'file' declaration in a 'set' parameter, for example:
 *      <pre>
 *      input:
 *      set( file('name'), ... )
 *      </pre>
 *
 */
class TokenFileCall {
    final target
    TokenFileCall(target)  { this.target = target }
}

/**
 *  A token used by the DSL to identify a 'path' declaration in a 'set' parameter, for example:
 *      <pre>
 *      input:
 *      set( path('name'), ... )
 *      </pre>
 *
 */
class TokenPathCall {
    final target
    final Map opts

    TokenPathCall(target) {
        this.target = target
        this.opts = Collections.emptyMap()
    }

    TokenPathCall(Map opts, target) {
        this.target = target
        this.opts = opts
    }
}

/**
 * An object of this class replace the {@code stdin} token in input map declaration. For example:
 * <pre>
 * input:
 *   map( stdin, .. ) from x
 * </pre>
 *
 * @see nextflow.ast.NextflowDSLImpl
 * @see nextflow.script.params.TupleInParam#bind(java.lang.Object[])
 */
class TokenStdinCall { }

/**
 * An object of this class replace the {@code stdout} token in input map declaration. For example:
 * <pre>
 * input:
 *   map( stdout, .. ) into x
 * </pre>
 *
 * @see nextflow.ast.NextflowDSLImpl
 * @see nextflow.script.params.TupleOutParam#bind(java.lang.Object[])
 */
class TokenStdoutCall { }

/**
 * Token used by the DSL to identify a environment variable declaration, like this
 *     <pre>
 *     input:
 *     set( env(X), ... )
 *     <pre>
 */
@ToString
@EqualsAndHashCode
@TupleConstructor
class TokenEnvCall {
    Object val
}


/**
 * This class is used to identify a 'val' when used like in this example:
 * <pre>
 *  input:
 *  set ( val(x), ...  )
 *
 *  output:
 *  set( val(y), ...  )
 *
 * </pre>
 *
 */
@ToString
@EqualsAndHashCode
@TupleConstructor
class TokenValCall {
    Object val
}


@ToString
@EqualsAndHashCode
class TokenValRef {
    String name
    int lineNum
    int colNum

    TokenValRef( String name, int lineNum = -1, int colNum = -1 ) {
        this.name = name
        this.lineNum = lineNum
        this.colNum = colNum
    }
}


@ToString
@EqualsAndHashCode
@TupleConstructor
@CompileStatic
class TokenBranchDef {
    Closure<TokenBranchChoice> closure
    List<String> branches
}

@ToString
@EqualsAndHashCode
@TupleConstructor
@CompileStatic
class TokenBranchChoice {
    Object value
    String choice
}

@ToString
@EqualsAndHashCode
@TupleConstructor
@CompileStatic
class TokenMultiMapDef {
    Closure<Map> closure
    List<String> names
}

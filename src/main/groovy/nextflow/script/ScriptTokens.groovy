/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.script
import groovy.transform.Canonical
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
/**
 * Presents a variable definition in the script context.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Canonical
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

    final Object target

    TokenFileCall( value ) {
        this.target = value
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
 * @see SetInParam#bind(java.lang.Object[])
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
 * @see SetOutParam#bind(java.lang.Object[])
 */
class TokenStdoutCall { }

/**
 * Token used by the DSL to identify a environment variable declaration, like this
 *     <pre>
 *     input:
 *     set( env(X), ... )
 *     <pre>
 */
@Canonical
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
@Canonical
class TokenValCall {
    Object val
}



/**
 * Holds process script meta-data
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskBody {

    /** The actual task code */
    Closure closure

    /** The task code as source, eventually used to report diagnostic message in case of error */
    String source

    /** The script type, either a system text script or groovy native code*/
    ScriptType type

    /** When true the script is interpreted a shell script in which $xxx variable are not interpreted by groovy */
    boolean isShell

    /**
     * The set of token variable and properties referenced in the task script
     */
    Set<TokenValRef> valRefs

    /**
     * Hold the use provided task body
     *
     * @param closure The task body closure
     * @param source The task source string
     * @param section The section in which it has been, declared. One of
     *      {@code exec}, {@code script} or {@code shell} literals
     */
    TaskBody( Closure closure, String source, String section = 'script' ) {
        this.closure = closure
        this.source = source
        setType(section)
    }

    TaskBody( Closure closure, String source, String section, TokenValRef... values ) {
        this(closure, source, section)
        this.valRefs = values != null ? values as Set : new HashSet<>()
    }

    /**
     * Given a code section literal sets the task body {@link #type} and {@link #isShell} flag
     *
     * @param section
     * @return
     */
    protected TaskBody setType( String section )  {

        switch( section ) {
            case 'exec':
                type = ScriptType.GROOVY
                isShell = false
                break

            case 'script':
                type = ScriptType.SCRIPTLET
                isShell = false
                break

            case 'shell':
                this.type = ScriptType.SCRIPTLET
                isShell = true
                break

            default:
                throw new IllegalArgumentException("Not a valid process body definition: $section")
        }

        return this
    }

    String toString() {
        "TaskBody[closure: $closure; types: $type; shell: $isShell; variables: $valRefs; source:\n$source \n]"
    }

    /**
     * @return
     *      The set of variable and properties referenced in the user script.
     *      NOTE: it includes properties in the form {@code object.propertyName}
     */
    Set<String> getValNames() {
        valRefs *. name
    }
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

/*
 * Copyright (c) 2012, the authors.
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

/**
 * Presents a variable definition in the script context.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Canonical
class ScriptVar {

    /** The variable name */
    String name

}

/**
 *
 *
 */
class ScriptFileWrap {

    final String name

    final String filePattern

    ScriptFileWrap( value ) {

        if( value instanceof ScriptVar ) {
            name = value.name
            filePattern = '*'
        }
        else if( value instanceof Map ) {
            def entry = value.entrySet().first()
            name = entry.key
            filePattern = entry.value?.toString()
        }
        else if( value instanceof String ) {
            name = value
            filePattern = value
        }
        else {
            throw new IllegalArgumentException()
        }

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
class ScriptStdinWrap { }

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
class ScriptStdoutWrap { }

@Canonical
class ScriptEnvWrap {
    String name
}

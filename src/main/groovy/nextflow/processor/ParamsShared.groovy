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

package nextflow.processor
import java.nio.file.Path

import groovy.transform.ToString
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.ast.ProcessVarRef
import nextflow.script.BaseScript
/**
 * 'SharedParam' Marker interface
 */
interface SharedParam extends InParam, OutParam {

}

/**
 * Model a process shared parameter that behaves both as input and output parameter
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
abstract class BaseSharedParam extends BaseParam implements SharedParam {

    protected Object target

    protected String name

    protected Object into

    private inChannel

    private outChannel


    protected BaseSharedParam( Script script ) {
        super(script)
    }

    BaseSharedParam setup() {

        /*
         * initialize *input* channel
         */
        if( target instanceof Closure ) {
            inChannel = sharedValToChannel( target.call() )
        }
        else if( target != null ) {
            inChannel = sharedValToChannel( target )
        }
        else {
            inChannel = sharedValToChannel( defValue() )
        }

        /*
         * initialize *output* channel
         */
        if( into != null ) {
            log.trace "shared output using > channel ref: $into"
            outChannel = outputValToChannel( script, into, DataflowVariable )
        }

        return this
    }

    /**
     * Implements the {@code as} keyword for the shared param declaration
     * NOTE: since {@code as} is a keyword for the groovy programming language
     * the method as to be named {@code _as}.
     *
     * A special pre-process will replace the "as" from the user script to the "_as"
     *
     * @see nextflow.ast.SourceModifierParserPlugin
     *
     * @param value
     * @return
     */
    SharedParam _as( Object value ) {
        alias(value)
    }

    /*
     * Just a synonym for the "as" method
     *
     * @param value
     * @return
     */
    protected BaseSharedParam alias( Object value ) {
        this.name = value
        return this
    }

    BaseSharedParam to( Object value ) {
        this.into = value
        return this
    }

    String getName() { name }

    DataflowReadChannel getInChannel() {
        lazyInit()
        inChannel
    }

    DataflowWriteChannel getOutChannel() {
        lazyInit()
        outChannel
    }

}



class ValueSharedParam extends BaseSharedParam {

    ValueSharedParam ( Script script, def val ) {
        super(script)

        if( val instanceof ProcessVarRef ) {
            this.name = val.name
            this.target = getScriptVar(val.name)
        }
        else {
            this.target = val
        }
    }

    @Override
    Object defValue() { null }

}





@ToString(includePackage=false, includeSuper = true)
class FileSharedParam extends BaseSharedParam {

    private String fileName

    @Lazy
    private Path tempFile = { ((BaseScript)script).tempFile(fileName) } ()

    FileSharedParam( Script script, Object val ) {
        super(script)

        if( val instanceof ProcessVarRef ) {
            this.name = val.name
            this.target = getScriptVar(val.name)
        }
        else if( val instanceof String ) {
            fileName = val
        }
    }

    @Override
    protected FileSharedParam alias( Object str ) {
        this.fileName = str
        return this
    }

    Path defValue() {
        tempFile
    }

}


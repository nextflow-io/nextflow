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
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import nextflow.script.BaseScript
/**
 * Model a process shared parameter that behaves both as input and output parameter
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@ToString(includePackage=false, includeNames = true)
abstract class SharedParam {

    final protected Script script

    protected String name

    protected Object using

    protected Object into

    @Lazy
    DataflowReadChannel input = { inputRef() }()

    /**
     * The channel over which entries are sent
     */
    @Lazy
    DataflowWriteChannel output = outputRef()


    SharedParam( Script script, String name ) {
        this.script = script
        this.name = name
    }

    SharedParam using( Object value ) {
        this.using = value
        return this
    }

    SharedParam into( Object value ) {
        this.into = value
        return this
    }


    def String getName() { name }


    abstract Object defValue()


    protected DataflowReadChannel inputRef() {

        def value
        if( using != null ) {
            value = using instanceof Closure ? using.call() : using
        }
        else if( script && script.getBinding().hasVariable(name) ) {
            value = script.getBinding().getVariable(name)
        }
        else {
            value = defValue()
        }

        return asChannel(value)
    }

    protected DataflowWriteChannel outputRef() {

        if( into != null ) {
            log.trace "shared output using > channel ref: $into"
            return OutParam.channelRef( script, into, { new DataflowVariable() } )
        }

        return null
    }


    DataflowReadChannel asChannel( def value ) {

        if( value instanceof DataflowExpression ) {
            return value
        }
        else if( value instanceof DataflowReadChannel ) {
            throw new IllegalArgumentException()
        }

        def result = new DataflowVariable()
        result.bind(value)
        result
    }

}



/**
 *  Shared value parameter.
 *  <p>
 *      {@code shared: val <name> }
 *      <br>
 *      Declares a thread safe variable in the execution context with the specified 'name'.
 *      When exists a variable with the same name in the script context, tried to use it
 *      as the initial value to which it binds on the input channel
 *      <br>
 *      No output channel is defined
 * <p>
 *      {@code shared: val <name> using obj }
 *      Declares a thread safe variable in the execution context with the specified 'name'.
 *      The initial value is taken from the obj declared by the 'using' keyword
 *
 * <p>
 *      {@code shared: val <name> using obj into channelName }
 *      Declares a thread safe variable in the execution context with the specified 'name'.
 *      The initial value is taken from the obj declared by the 'using' keyword.
 *      When processor completes, bind the final value to the channel specified by the 'into' keyword
 *
 */
@ToString(includePackage=false, includeSuper = true)
class ValueSharedParam extends SharedParam {

    ValueSharedParam ( Script script, String name ) {
        super(script,name)
    }

    @Override
    Object defValue() { null }

}




/**
 *  Model a process *file* input parameter
 */
@ToString(includePackage=false, includeSuper = true)
class FileSharedParam extends SharedParam  {

    String filePattern

    FileSharedParam( Script script, String name ) {
        super(script,name)
        this.filePattern = name
    }

    Object defValue() {
        ((BaseScript)script).tempFile(name)
    }

    @Override
    protected DataflowReadChannel inputRef() {

        if( !into ) {
            this.filePattern = '*'
        }
        super.inputRef()

    }


}



/**
 * Container to hold all process outputs
 */
class SharedList implements List<SharedParam> {

    @Delegate
    List<SharedParam> target = new LinkedList<>()

    List<DataflowReadChannel> getInChannels() { target *.input }

    List<DataflowWriteChannel> getOutChannels() { target *.output }

    List<String> getNames() { target *. name }


}


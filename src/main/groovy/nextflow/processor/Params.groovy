package nextflow.processor

import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Nextflow

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
abstract class InParam {

    DataflowReadChannel channel

    String name

    /**
     * Parameter factory method, given the map of attributes create the respective input parameter object
     * <p>
     *     For example {@code file: 'input.fasta', from: my_channel} will create a
     *     {@code InFileParam} instance
     *
     *
     * @param args
     * @return
     */
    static InParam parse( Map args ) {
        assert args

        if( !args.from )  {
            throw new IllegalArgumentException('Missing \'from\' definition in input declaration')
        }


        if( args.file ) {
            return new InFileParam( name: args.file, channel: asChannel(args.from) )
        }

        if( args.env ) {
            return new InEnvParam( name: args.env, channel: asChannel(args.from) )
        }

        if( args.val ) {
            return new InValueParam( name:args.val, channel: asChannel(args.from) )
        }

        throw new IllegalArgumentException("Unknown 'in' definition: $args")
    }


    static DataflowReadChannel asChannel( def value ) {

        if ( value instanceof DataflowBroadcast )  {
            return value.createReadChannel()
        }

        if( value instanceof DataflowReadChannel ) {
            return value
        }

        // wrap any collections with a DataflowQueue
        if( value instanceof Collection ) {
            return Nextflow.channel(value)
        }

        // wrap any array with a DataflowQueue
        if ( value && value.class.isArray() ) {
            return Nextflow.channel(value as List)
        }

        // wrap a single value with a DataflowVariable
        return Nextflow.val(value)

    }

}


class InFileParam extends InParam  { }

class InEnvParam extends InParam { }

class InValueParam extends InParam { }


class InputList implements List<InParam> {

    @Delegate
    List<InParam> target = new LinkedList<>()

    List<DataflowReadChannel> getChannels() { target *.channel }

    List<String> getNames() { target *. name }

    List<InParam> ofType( Class... classes ) { target.findAll { it.class in classes } }

    void eachParam (Closure closure) {
        target.each { InParam param -> closure.call(param.name, param.channel) }
    }

}


abstract class OutParam {

    String name

    DataflowWriteChannel channel

}

class OutFileParam extends OutParam {

    boolean grouped

}

class OutputList implements List<OutParam> {

    @Delegate
    List<OutParam> target = new LinkedList<>()

    List<DataflowWriteChannel> getChannels() { target *.channel }

    List<String> getNames() { target *. name }

    void eachParam (Closure closure) {
        target.each { OutParam param -> closure.call(param.name, param.channel) }
    }

}
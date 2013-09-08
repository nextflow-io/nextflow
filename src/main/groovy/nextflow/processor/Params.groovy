package nextflow.processor

import groovy.transform.ToString
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.Nextflow

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@ToString(includePackage=false, includeNames = true)
abstract class InParam {

    String name

    DataflowReadChannel channel


    /**
     * Parameter factory method, given the map of attributes create the respective input parameter object
     * <p>
     *     For example {@code file: 'input.fasta', from: my_channel} will create a
     *     {@code FileInParam} instance
     *
     *
     * @param args
     * @return
     */
    static InParam create( Map args ) {
        log.debug "Input ${args}"
        assert args

        if( !args.from )  {
            throw new IllegalArgumentException('Missing \'from\' definition in input declaration')
        }


        if( args.file ) {
            def nm = args.file
            def ch = asChannel(args.from)
            def obj = nm == '-' ? new StdInParam() : new FileInParam()

            obj.name = nm
            obj.channel = ch

            return obj
        }

        if( args.env ) {
            return new EnvInParam( name: args.env, channel: asChannel(args.from) )
        }

        if( args.val ) {
            return new ValueInParam( name:args.val, channel: asChannel(args.from) )
        }

        if( args.each ) {
            if ( !(args.from instanceof Collection) ) {
                throw new IllegalArgumentException("Illegal form: ${args.from} for 'each' input iterator")
            }

            // the collection is wrapper like a scalar value
            def wrap = Nextflow.val( args.from )

            return new EachInParam( name: args.each, channel: wrap )
        }

        throw new IllegalArgumentException("Illegal'input' definition: $args")
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

@ToString(includePackage=false)
class FileInParam extends InParam  { }

@ToString(includePackage=false)
class EnvInParam extends InParam { }

@ToString(includePackage=false)
class ValueInParam extends InParam { }

@ToString(includePackage=false)
class StdInParam extends InParam { StdInParam() { name='-' } }

@ToString(includePackage=false)
class EachInParam extends InParam {}

@Slf4j
@ToString(includePackage=false, includeNames = true)
abstract class OutParam {

    /** The out parameter name */
    String name

    /** The channel over which entries are sent */
    DataflowWriteChannel channel

    /** Whenever the channel has to closed on task termination */
    Boolean autoClose = Boolean.TRUE


    /**
     * Output parameter factory method. Given the map of named attributes for an output item
     * create the right instance of {@code OutParam} object
     *
     * @param attributes
     * @param script
     * @return
     */
    static OutParam create( Map attributes, Script script ) {
        log.debug "Output: ${attributes}"
        assert attributes

        def OutParam result
        def name = attributes.file as String
        if( name == '-' ) {
            result = new StdOutParam(name:name)
        }
        else if( name ) {
            result = new FileOutParam(name:name)
        }
        else {
            name = attributes.val
            if( !name ) throw new IllegalStateException("Missing output parameter name")
            result = new ValueOutParam(name:name)
        }


        def channel = attributes.into

        // the receiving channel may not be defined explicitly
        // in that case the specified file name will used to
        // reference it in the script context
        if( channel == null ) {
            log.trace "output > channel not defined"
            channel = ( name != '-' ) ? name.replaceAll(/\./, '_') :  Nextflow.channel()
        }

        if( channel instanceof String ) {
            // the channel is specified by name
            def local = channel
            log.trace "output > channel name: $local"

            // look for that name in the 'script' context
            channel = getScriptVariable(script, local)
            if( channel instanceof DataflowWriteChannel ) {
                // that's OK -- nothing to do
            }
            else {
                if( channel == null ) {
                    log.debug "output > channel unknown: $local -- creating a new instance"
                }
                else {
                    log.warn "Duplicate output channel name: '$channel' in the script context -- it's worth to rename it to avoid possible conflicts"
                }

                // instantiate the new channel
                channel = Nextflow.channel()
                // bind it to the script on-fly
                if( local != '-' && script) {
                    // bind the outputs to the script scope
                    script.setProperty(local, channel)
                }
            }
        }

        // the special file name '-' has to be mapped to special parameter ot type StdOutParam
        result.channel = channel


        /*
         * set any other attributes on the target object other than 'file' and 'into'
         */
        attributes.each { String key, value ->
            if( key in ['file','into', 'val'] || value == null ) { return }

            def prop = result.getMetaClass().getMetaProperty(key)
            if( !prop ) {
                log.warn "Invalid output attribute name: '$key'"
                return
            }

            if( prop.getType().isInstance(value) ) {
                prop.setProperty(result, value)
            }
            else {
                log.warn "Output type mismatch -- cannot assign value: '$value' to output attribute: '$key'"
            }
        }


        return result
    }



    /*
     * Try to access to the script property defined by the
     * specified {@code name}
     * <p>
     * If it's not defined, it will return {@code null}
     *
     */
    static private getScriptVariable( Script script, String name ) {
        assert script
        try {

            script.getBinding().getVariable(name)
        }
        catch( MissingPropertyException e ) {
            return null
        }
    }



}

@ToString(includePackage=false, includeSuper = true, includeNames = true)
class FileOutParam extends OutParam {

    /**
     * Whenever multiple files matching the same name pattern have to be output as grouped collection
     * or bound as single entries
     */
    Boolean joint = Boolean.FALSE

    /**
     * The character used to separate multiple names (pattern) in the output specification
     */
    String separatorChar = ':'

}

@ToString(includePackage=false, includeSuper = true, includeNames = true)
class StdOutParam extends OutParam { StdOutParam() { name='-' } }


@ToString(includePackage=false, includeSuper = true, includeNames = true)
class ValueOutParam extends OutParam { }



class InputsList implements List<InParam> {

    @Delegate
    List<InParam> target = new LinkedList<>()

    List<DataflowReadChannel> getChannels() { target *.channel }

    List<String> getNames() { target *. name }

    def <T extends InParam> List<T> ofType( Class<T> clazz ) { (List<T>) target.findAll { it.class == clazz } }

    void eachParam (Closure closure) {
        target.each { InParam param -> closure.call(param.name, param.channel) }
    }

}

class OutputsList implements List<OutParam> {

    @Delegate
    List<OutParam> target = new LinkedList<>()

    List<DataflowWriteChannel> getChannels() { target *.channel }

    List<String> getNames() { target *. name }

    def <T extends OutParam> List<T> ofType( Class<T> clazz ) { (List<T>) target.findAll { it.class == clazz } }

    void eachParam (Closure closure) {
        target.each { OutParam param -> closure.call(param.name, param.channel) }
    }

}

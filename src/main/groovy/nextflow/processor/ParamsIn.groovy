package nextflow.processor

import groovy.transform.InheritConstructors
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowReadChannel
import nextflow.Nextflow
/**
 * Model a process generic input parameter
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@ToString(includePackage=false, includeNames = true)
abstract class InParam {

    final protected Script script

    protected String name

    protected Object using

    @Lazy
    DataflowReadChannel channel = { getLazyChannel() } ()


    protected DataflowReadChannel getLazyChannel() {

        if( using ) {
            return asChannel(using)
        }
        else if( script && script.getBinding().hasVariable(name) ) {
            return asChannel(script.getBinding().getVariable(name))
        }

        throw new IllegalStateException("Missing channel for input parameter: $name")

    }

    InParam( Script script, String name ) {
        this.script = script
        this.name = name
    }

    InParam using( Object value ) {
        this.using = value
        return this
    }

    def String getName() { name }



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

/**
 *  Model a process *file* input parameter
 */
@ToString(includePackage=false, includeSuper = true)
class FileInParam extends InParam  {

    String filePattern

    FileInParam( Script script, String name ) {
        super(script,name)
        this.filePattern = name
    }

    @Override
    protected DataflowReadChannel getLazyChannel() {

        if( !using ) {
            this.filePattern = '*'
        }
        super.getLazyChannel()

    }


}

/**
 *  Model a process *environment* input parameter
 */
@InheritConstructors
@ToString(includePackage=false, includeSuper = true)
class EnvInParam extends InParam { }

/**
 *  Model a process *value* input parameter
 */
@InheritConstructors
@ToString(includePackage=false, includeSuper = true)
class ValueInParam extends InParam { }

/**
 *  Model a process *stdin* input parameter
 */
@ToString(includePackage=false, includeSuper = true)
class StdInParam extends InParam { StdInParam(Script script) { super(script,'-') } }

/**
 *  Model a process input *iterator* parameter
 */
@InheritConstructors
@ToString(includePackage=false, includeSuper = true)
class EachInParam extends InParam {

    def EachInParam using( Object value ) {

        // everything is mapped to a collection
        def list = Nextflow.list(value)
        // the collection is wrapped to a "scalar" dataflow variable
        super.using( Nextflow.val(list) )
        return this

    }

}



/**
 * Container to hold all process outputs
 */
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


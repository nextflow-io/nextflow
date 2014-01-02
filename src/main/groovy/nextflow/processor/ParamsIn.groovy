package nextflow.processor
import groovy.transform.InheritConstructors
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.expression.DataflowExpression
import nextflow.Nextflow
import nextflow.ast.ProcessVarRef

@Slf4j
abstract class BaseParam {

    final protected Script script

    protected BaseParam ( Script script ) {
        this.script = script
    }

    /**
     * Initialize the parameter before use it
     */
    abstract BaseParam setup()

    private boolean initialized

    final protected void lazyInit() {
        if( initialized ) return
        setup()
        initialized = true
    }

    /**
     * Get the value of variable {@code name} in the script context
     *
     * @param name The variable name
     * @param strict If {@code true} raises a {@code MissingPropertyException} when the specified variable does not exist
     * @return The variable object
     */
    protected getScriptVar(String name, boolean strict = false) {
        if( script.getBinding().hasVariable(name) ) {
            return script.getBinding().getVariable(name)
        }

        if( strict )
            throw new MissingPropertyException(name,this.class)

        return null
    }

    final protected DataflowReadChannel inputValToChannel( def value ) {

        if ( value instanceof DataflowBroadcast )  {
            return value.createReadChannel()
        }

        if( value instanceof DataflowReadChannel ) {
            return value
        }

        // wrap any collections with a DataflowQueue
        if( value instanceof Collection ) {
            return Nextflow.channel(value as List)
        }

        // wrap any array with a DataflowQueue
        if ( value && value.class.isArray() ) {
            return Nextflow.channel(value as List)
        }

        // wrap a single value with a DataflowVariable
        return Nextflow.val(value)

    }

    final protected DataflowWriteChannel outputValToChannel( Script script, Object channel, Class<DataflowWriteChannel> factory ) {

        if( channel instanceof String ) {
            // the channel is specified by name
            def local = channel

            def binding = script.getBinding()

            // look for that name in the 'script' context
            channel = binding.hasVariable(local) ? binding.getVariable(local) : null
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
                channel = factory.newInstance()

                // bind it to the script on-fly
                if( local != '-' && script) {
                    // bind the outputs to the script scope
                    binding.setVariable(local, channel)
                }
            }
        }

        if( channel instanceof DataflowWriteChannel ) {
            return channel
        }

        throw new IllegalArgumentException("Invalid output channel reference")
    }

    final DataflowReadChannel sharedValToChannel( def value ) {

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

interface InParam {

    String getName()

    DataflowReadChannel getInChannel()

    InParam _as( Object value )

    Object defValue()

}

/**
 * Model a process generic input parameter
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@Slf4j
@ToString(includePackage=false, includeNames = true)
abstract class BaseInParam extends BaseParam implements InParam {

    /**
     * The name used to bind this parameter in the process execution context (Map)
     */
    protected String name

    /**
     * the target object which hold the value of the parameter to bind in the process execution context
     */
    protected Object target

    Object defValue() {
        throw new IllegalStateException()
    }

    private inChannel

    DataflowReadChannel getInChannel() {
        lazyInit()
        return inChannel
    }

    /**
     * @param script The global script object
     * @param obj
     */
    protected BaseInParam( Script script, obj ) {
        super(script)

        if( obj instanceof ProcessVarRef ) {
            // when the value is a variable reference
            // - use that name for the parameter itself
            // - get the variable value in the script binding
            this.name = obj.name
            this.target = getScriptVar(obj.name, true)  // <-- true: raise an MissingPropertyException when it does not exist
        }
        else {
            // just set the parameter to the specified value
            // it must follow an 'as' keyword in the parameter definition
            this.target = obj
        }
    }

    @Override
    BaseInParam setup() {

        // initialize the *inChannel* object based on the 'target' attribute
        def result
        if( target instanceof Closure ) {
            result = target.call()
        }
        else if( target != null ) {
            result = target
        }
        else {
            // fallback on the default value
            result = defValue()
        }

        inChannel = inputValToChannel(result)

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
    InParam _as( Object value ) {
        alias(value)
    }

    /*
     * Just a synonym for the "as" method
     *
     * @param value
     * @return
     */
    protected InParam alias( Object value ) {
        this.name = value
        return this
    }

    def String getName() { name }

}

/**
 *  Model a process *file* input parameter
 */
@ToString(includePackage=false, includeSuper = true)
class FileInParam extends BaseInParam  {

    FileInParam( Script script, obj ) {
        super(script, obj)

        if( obj instanceof String ) {
            name = obj
        }

        if( !name ) {
            name = '*'
        }
    }

}

/**
 *  Model a process *environment* input parameter
 */
@InheritConstructors
@ToString(includePackage=false, includeSuper = true)
class EnvInParam extends BaseInParam { }

/**
 *  Model a process *value* input parameter
 */
@ToString(includePackage=false, includeSuper = true)
class ValueInParam extends BaseInParam {

    ValueInParam( Script script, Object val ) {
        super(script, val)
    }

}

/**
 *  Model a process *stdin* input parameter
 */
@ToString(includePackage=false, includeSuper = true)
class StdInParam extends BaseInParam {

    StdInParam( Script script, def obj ) {
        super(script,obj)
        this.name = '-'
    }

    def StdInParam alias( Object obj ) {
        throw new IllegalAccessException("keyword 'as' not supported for 'stdin' definition")
    }

}

/**
 *  Model a process input *iterator* parameter
 */
@ToString(includePackage=false, includeSuper = true)
class EachInParam extends BaseInParam {

    def EachInParam( Script script, Object value ) {
        super(script,value)

        // everything is mapped to a collection
        // the collection is wrapped to a "scalar" dataflow variable
        def list = Nextflow.list(target)
        target = Nextflow.val(list)
    }

}



/**
 * Container to hold all process outputs
 */
class InputsList implements List<InParam> {

    @Delegate
    private List<InParam> target = new LinkedList<>()

    List<DataflowReadChannel> getChannels() {
        target.collect { InParam it -> it.getInChannel() }
    }

    List<String> getNames() { target *. name }


    def <T extends InParam> List<T> ofType( Class<T> clazz ) {
        (List<T>) target.findAll { it.class == clazz }
    }


}


package nextflow.extension

import java.lang.reflect.Method
import java.lang.reflect.ParameterizedType
import java.util.concurrent.Callable

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowBroadcast
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import nextflow.NF
import nextflow.dag.NodeMarker
import nextflow.exception.ScriptRuntimeException
import nextflow.script.ChannelOut
import org.codehaus.groovy.runtime.InvokerHelper
/**
 * Represents an nextflow operation invocation
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class OpCall implements Callable {

    final static private List<String> SPECIAL_NAMES = ["choice","merge","separate"]

    final static private String SET_OP_hack = 'set'

    static ThreadLocal<OpCall> current = new ThreadLocal<>()

    private OperatorEx owner
    private String methodName
    private Method method
    private Object source
    private Object[] args
    private Set inputs = new HashSet(5)
    private Set outputs = new HashSet<>(5)
    boolean ignoreDagNode

    static OpCall create(String methodName, Object args) {
        new OpCall(methodName, InvokerHelper.asArray(args))
    }

    static OpCall create(String methodName) {
        new OpCall(methodName, InvokerHelper.asArray(null))
    }

    OpCall(OperatorEx owner, Object source, String method, Object[] args ) {
        assert owner
        assert method
        this.owner = owner
        this.methodName = method
        this.args = ChannelOut.spread(args).toArray()
        this.setSource(source)
    }

    OpCall(String method, Object[] args ) {
        assert method
        this.owner = OperatorEx.instance
        this.methodName = method
        this.args = ChannelOut.spread(args).toArray()
    }

    OpCall setSource(ChannelOut left) {

        if( methodName == SET_OP_hack ) {
            source = left
            return this
        }

        if( left.size()== 0 ) {
            throw new ScriptRuntimeException("Operator '${methodName}' cannot be applied to an undefined output")
        }
        if( left.size()==1 ) {
            this.source = left[0] as DataflowWriteChannel
            return this
        }

        if( args.size() )
            throw new ScriptRuntimeException("Multi-channel output cannot be applied to operator ${methodName} for which argument is already provided")

        source = left[0] as DataflowWriteChannel
        args = left[1..-1] as Object[]
        return this
    }

    OpCall setSource( obj ) {
        if( obj instanceof ChannelOut )
            return setSource(obj)
        else
            source = obj
        return this
    }

    OpCall setSource(DataflowWriteChannel channel) {
        this.source = channel
        return this
    }

    OpCall setArgs(Object[] args) {
        this.args = args
        return this
    }

    @Override
    Object call() throws Exception {
        if( source==null )
            throw new IllegalStateException("Missing operator source channel")
        current.set(this)
        try {
            return invoke()
        }
        finally {
            current.remove()
        }
    }

    Set getInputs() { inputs }

    Set getOutputs() { outputs }

    String getMethodName() { methodName }

    Object[] getArgs() { args }

    private <T> T read0(source){
        if( source instanceof DataflowBroadcast )
            return (T)CH.getReadChannel(source)

        if( source instanceof DataflowQueue )
            return (T)CH.getReadChannel(source)

        else
            return (T)source
    }

    private Object[] read1(Object[] args) {
        if( methodName != 'separate' && methodName != 'choice' ) {
            Object[] params = new Object[args.length]
            for( int i=0; i<args.length; i++ )
                params[i] = read0(args[i])

            return params
        }
        return args
    }


    protected Object invoke() {
        if( methodName==SET_OP_hack ) {
            // well this is ugly, the problem is that `set` is not a real operator
            // but it's exposed as such. let's live whit this for now
            return invoke1('set', [source, args[0]] as Object[])
        }

        final DataflowReadChannel source = read0(source)
        final result = invoke0(source, read1(args))
        if( !ignoreDagNode )
            addGraphNode(result)
        return result
    }

    protected void addGraphNode(Object result) {
        // infer inputs
        inputs.add( source )
        inputs.addAll( getInputChannels() )
        outputs.addAll( getOutputChannels(result))

        NodeMarker.addOperatorNode(methodName, inputs, outputs)
    }

    protected List getInputChannels() {
        def result = new ArrayList(5)

        if( declaresParamType(DataflowReadChannel, method) ) {
            for( int i=0; i<args.length; i++ ) {
                fetchChannels(args[i], result)
            }
        }

        return result
    }

    protected List getOutputChannels(Object value) {
        def result = new ArrayList(5)

        if( declaresParamType(DataflowWriteChannel, method) ) {
            for( int i=0; i<args.length; i++ ) {
                fetchChannels(args[i], result)
            }
        }

        if( declaresReturnType(DataflowWriteChannel, method) ) {
            fetchChannels(value, result)
        }
        else if( declaresReturnType(ChannelOut, method) ) {
            fetchChannels(value, result)
        }

        return result
    }

    protected void fetchChannels(Object item, List result) {
        if( item instanceof DataflowWriteChannel )
            result.add(item)
        else if( item instanceof List ) {
            final arr = (item as List)
            for( int i=0; i<arr.size(); i++ ) {
                if( arr[i] instanceof DataflowWriteChannel )
                    result.add(arr[i])
            }
        }
        else if( item instanceof Object[] ) {
            final arr = (item as Object[])
            for( int i=0; i<arr.size(); i++ ) {
                if( arr[i] instanceof DataflowWriteChannel )
                    result.add(arr[i])
            }
        }
    }

    protected boolean declaresReturnType(Class type, Method method) {
        def class1 = method.getReturnType()
        if( type.isAssignableFrom(class1) )
            return true

        if( class1.isArray() && type.isAssignableFrom(class1.getComponentType()) )
            return true

        if( List.isAssignableFrom(class1) ) {
            if( method.getGenericReturnType() instanceof ParameterizedType ) {
                final param = (ParameterizedType)method.getGenericReturnType()
                final class2 = (Class)param.getActualTypeArguments()[0]
                return type.isAssignableFrom(class2)
            }
        }

        return false
    }

    protected boolean declaresParamType(Class type, Method method) {
        def allParams = method.getParameterTypes()
        // note: starts from 1 because the first parameter is expected
        // to be a DataflowReadChannel by definition
        for( int i=1; i<allParams.length; i++ ) {
            if( declaresParamType(type, method, i))
                return true
        }
        return false
    }

    protected boolean declaresParamType(Class type, Method method, int index) {
        final paramTypes = method.getParameterTypes()
        final class1 = paramTypes[index]
        if( type.isAssignableFrom(class1) )
            return true

        if( class1.isArray() && type.isAssignableFrom(class1.getComponentType()) )
            return true

        if( List.isAssignableFrom(class1) ) {
            final class2 = method.getGenericParameterTypes()[index]
            if( class2 instanceof ParameterizedType ) {
                final param = (ParameterizedType)class2
                final class3 = (Class)param.getActualTypeArguments()[0]
                return type.isAssignableFrom(class3)
            }
        }
        return false
    }


    protected Object invoke0(DataflowReadChannel channel, Object[] args) {

        if (checkOpenArrayDataflowMethod(SPECIAL_NAMES, methodName, args)) {
            // make it possible to invoke some dataflow operator specifying an open array and ending with a closure argument
            // For example:
            //  DataflowReadChannel# separate(final List<DataflowWriteChannel<?>> outputs, final Closure<List<Object>> code)
            // can be invoked as:
            //  queue.separate( x, y, z ) { ... }

            Object[] params = new Object[3]
            params[0] = channel
            params[1] = toListOfChannel(args)
            params[2] = args[args.length - 1]

            return invoke1(methodName, params)
        }
        else {
            // create the invocation parameters array
            Object[] params = new Object[args.length + 1]
            params[0] = channel
            for (int i = 0; i < args.length; i++) {
                params[i + 1] = args[i];
            }

            return invoke1(methodName, params);
        }
    }

    protected Method getMethod0(String methodName, Object[] args) {
        def meta = owner.metaClass.getMetaMethod(methodName, args)
        if( meta == null )
            throw new MissingMethodException(methodName, owner.getClass())
        method = owner.getClass().getMethod(methodName, meta.getNativeParameterTypes())
    }

    protected Object invoke1(String methodName, Object[] args) {
        method = getMethod0(methodName, args)
        checkDeprecation(method)
        return owner.metaClass.invokeMethod(owner, methodName, args)
    }

    protected void checkDeprecation(Method method) {
        if( method.getAnnotation(Deprecated) ) {
            log.warn "Operator `$methodName` is deprecated -- it will be removed in a future release"
        }
        else if( method.getAnnotation(DeprecatedDsl2) && NF.isDsl2() ) {
            def annot = method.getAnnotation(DeprecatedDsl2)
            def messg = annot.message() ?: "Operator `$methodName` is deprecated -- it will be removed in a future release".toString()
            log.warn messg
        }
    }

    @CompileStatic
    private static boolean checkOpenArrayDataflowMethod(List<String> validNames, String methodName, Object[] args) {
        if( !validNames.contains(methodName) )
            return false
        if( args == null || args.length<2 )
            return false
        if( !(args[args.length-1] instanceof Closure) )
            return false
        for( int i=0; i<args.length-1; i++ ) {
            if( !(args[i] instanceof DataflowWriteChannel ))
                return false
        }

        return true
    }


    /**
     * Given an array of arguments copy the first n-1 to a list of {@link groovyx.gpars.dataflow.DataflowWriteChannel}
     */
    @CompileStatic
    private static List<DataflowWriteChannel> toListOfChannel( Object[] args )  {
        List<DataflowWriteChannel> result = new ArrayList<>(args.length-1);
        for( int i=0; i<args.length-1; i++ ) {
            result.add(i, (DataflowWriteChannel)args[i]);
        }
        return result;
    }

}

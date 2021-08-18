package groovy.runtime.metaclass;

/**
 * Interface to define channel factory extensions e.g.
 *
 *   Channel.foo.fromSomething(params)
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public interface ChannelFactory {

    Object invokeExtensionMethod(String method, Object[] args);

}

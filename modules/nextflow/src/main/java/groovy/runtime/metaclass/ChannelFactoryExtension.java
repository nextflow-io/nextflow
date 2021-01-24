package groovy.runtime.metaclass;

/**
 * Interface to create Nextflow channel factory extensions e.g
 *
 *   Channel.foo.fromSomething(params)
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public interface ChannelFactoryExtension {

    Object invokeExtensionMethod(String method, Object[] args);

}

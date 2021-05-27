package nextflow.util

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

/**
 * Implements helper methods for {@link Process} objects
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ProcessHelper {

    static long pid(Process process) {
        try {
            // with Java 9 and later pid can be access with `.toHandle().getPid()` method
            final toHandle = Process.class.getMethod("toHandle")
            final handle = toHandle.invoke(process)
            return (Long) Class.forName("java.lang.ProcessHandle").getMethod("pid").invoke(handle)
        }
        catch (Exception e) {
            // ignore it and fallback on next
        }

        // Fallback for Java6+ on unix. This is known not to work on Windows.
        try {
            final pidField = process.class.getDeclaredField("pid")
            pidField.setAccessible(true)
            return pidField.getLong(process)
        }
        catch(Exception e) {
            throw new UnsupportedOperationException("Unable to access process pid for class: ${process.getClass().getName()}",e )
        }
    }

}

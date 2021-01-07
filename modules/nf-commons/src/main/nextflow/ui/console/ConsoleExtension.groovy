package nextflow.ui.console


import org.pf4j.ExtensionPoint
/**
 * Define the extension interface for console app
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
interface ConsoleExtension extends ExtensionPoint {

    void run(String...args)

}

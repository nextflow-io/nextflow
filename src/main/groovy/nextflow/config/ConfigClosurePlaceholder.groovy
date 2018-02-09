package nextflow.config

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.PackageScope

/**
 * Placeholder class that replacing closure definitions in the nextflow configuration
 * file in order to print the closure content itself
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@EqualsAndHashCode
@CompileStatic
@PackageScope
class ConfigClosurePlaceholder {

    private String str

    ConfigClosurePlaceholder(String str) {
        this.str = str
    }

    @Override String toString() { str }
}

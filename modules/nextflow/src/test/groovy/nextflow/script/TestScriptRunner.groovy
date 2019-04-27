package nextflow.script

import groovy.transform.InheritConstructors
import test.TestHelper

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@InheritConstructors
class TestScriptRunner extends ScriptRunner {

    TestScriptRunner setScript( String str ) {
        def file = TestHelper.createInMemTempFile('main.nf', str)
        super.setScript(file)
        return this
    }

}

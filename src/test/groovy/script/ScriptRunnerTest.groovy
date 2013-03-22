package script

import nextflow.script.CliRunner
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ScriptRunnerTest extends Specification {



    def 'test version' () {

        when:
        def opt = CliRunner.parseMainArgs('-v')

        then:
        assert opt.version

    }

    def 'test help' () {

        when:
        def opt = CliRunner.parseMainArgs('-h')

        then:
        assert opt.help

    }



    def 'test usage' () {

        when:
        CliRunner.parseMainArgs([] as String)
        CliRunner.jcommander.usage()

        then:
        noExceptionThrown()

    }

    def 'validateWorkDir' () {

        /*
         * the user specified a path that does not exist
         * it will created
         */
        when:
        def path = CliRunner.validateWorkDirectory( 'x1/testValidateFolder', 'scriptName' )

        then:
        path == new File('x1/testValidateFolder').canonicalFile
        path.exists()

        cleanup:
        path?.delete()

    }

    def 'validateWorkDir 2' () {

        setup:
        def folder = new File('x2/testValidateFolder')
        folder.mkdirs()
        def file = File.createTempFile('test','test',folder)

        /*
         * The folder specified by the user contains a file
         * -> it is not a valid folder, an exception is thrown
         */
        when:
        def path = CliRunner.validateWorkDirectory( folder.toString(), 'scriptName' )

        then:
        thrown(CliRunner.CliArgumentException.class)

        cleanup:
        folder?.delete()
        file?.delete()

    }

    def 'validateWorkDir 3' () {

        /*
         * No folder is specified, create a temp folder in the current dir, using the script name
         */
        when:
        def path = CliRunner.validateWorkDirectory( null, 'myScript' )

        then:
        path == new File('run-myScript').canonicalFile
        path.exists()

        cleanup:
        path?.delete()

    }


}

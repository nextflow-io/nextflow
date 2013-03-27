package script

import nextflow.script.CliRunner
import org.apache.commons.io.FileUtils
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CliRunnerTest extends Specification {


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
        def path = CliRunner.validateWorkDirectory( 'x1_testValidateFolder', 'scriptName' )

        then:
        path == new File('x1_testValidateFolder').canonicalFile
        path.exists()

        cleanup:
        if( path ) FileUtils.deleteDirectory(path)

    }

    def 'validateWorkDir 2' () {

        setup:
        def folder = new File('x2_testValidateFolder')
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
        if( folder ) FileUtils.deleteDirectory(folder)

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

    def 'buildConfigObject' () {

        setup:
        def env = [PATH:'/local/bin', HOME:'/home/my']

        when:
        def config1 = CliRunner.buildConfig0(env,false, null)
        def config2 = CliRunner.buildConfig0(env,true, null)

        then:
        ! ('PATH' in config1.env )
        ! ('HOME' in config1.env )
        ! ('XXX' in config1.env )

        ('PATH' in config2.env )
        ('HOME' in config2.env )
        !('XXX' in config2.env )

        config2.env.PATH == '/local/bin'
        config2.env.HOME == '/home/my'

    }

    def 'buildConfigObject 2 ' () {

        setup:
        def env = [HOME:'/home/my', PATH:'/local/bin', 'dot.key.name':'any text']

        def text1 = '''
        task { field1 = 1; field2 = 'hola'; }
        env { alpha = 'a1'; beta  = 'b1'; HOME="$HOME:/some/path"; }
        '''

        def text2 = '''
        task { field2 = 'Hello' }
        env { beta = 'b2'; delta = 'd2'; HOME="$HOME:/local/path"; XXX="$PATH:/xxx"; YYY = "$XXX:/yyy"   }
        '''

        when:
        def config1 = CliRunner.buildConfig0(env,true, [text1])
        def config2 = CliRunner.buildConfig0(env,false, [text1, text2])
        def config3 = CliRunner.buildConfig0(env,true, [text1, text2])

        // note: configuration object can be modified like any map
        config3.env ['ZZZ'] = '99'

        then:
        config1.task.field1 == 1
        config1.task.field2 == 'hola'
        config1.env.HOME == '/home/my:/some/path'
        config1.env.PATH == '/local/bin'
        config1.env.'dot.key.name' == 'any text'

        config2.task.field1 == 1
        config2.task.field2 == 'Hello'
        config2.env.alpha == 'a1'
        config2.env.beta == 'b2'
        config2.env.delta == 'd2'
        config2.env.HOME == '/home/my:/local/path'
        config2.env.XXX == '/local/bin:/xxx'
        !('PATH' in config2.env)

        config3.task.field1 == 1
        config3.task.field2 == 'Hello'
        config3.env.alpha == 'a1'
        config3.env.beta == 'b2'
        config3.env.delta == 'd2'
        config3.env.HOME == '/home/my:/local/path'
        config3.env.XXX == '/local/bin:/xxx'
        config3.env.PATH == '/local/bin'
        config3.env.YYY == '/local/bin:/xxx:/yyy'
        config3.env.ZZZ == '99'

    }


    def 'test validateConfigFiles '() {
        when:
        def files = CliRunner.validateConfigFiles(['file1','file2'])

        then:
        thrown(CliRunner.CliArgumentException)

        when:
        def f1 = File.createTempFile('file1','x').absoluteFile
        def f2 = File.createTempFile('file1','x').absoluteFile
        files = CliRunner.validateConfigFiles([f1.toString(), f2.toString()])

        then:
        files == [f1, f2]

        cleanup:
        f1?.delete()
        f2?.delete()

    }


    def 'test config to map ' () {

        setup:
        def config = new ConfigSlurper().parse( 'task {field1=1; field2="two"}; env { x = 99 } ' )

        when:
        def map = CliRunner.configToMap(config)
        map.env.PATH = '/local/bin'

        then:
        map.task.field1 == 1
        map.task.field2 == "two"
        map.env.x == 99
        map.env.y == null
        map.env.PATH == '/local/bin'

    }




}

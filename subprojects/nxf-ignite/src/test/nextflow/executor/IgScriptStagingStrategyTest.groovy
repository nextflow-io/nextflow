package nextflow.executor
import java.nio.file.Paths

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class IgScriptStagingStrategyTest extends Specification {

    def 'should check delegate methods' () {
        given:
        def PATH = Paths.get('/any/path')
        def delegate = new SimpleFileCopyStrategy()
        def strategy = new IgScriptStagingStrategy(delegate: delegate)

        expect:
        strategy.getBeforeStartScript() == null
        strategy.getStageInputFilesScript(['foo.txt': PATH]) == null
        strategy.getUnstageOutputFilesScript(['foo.txt'], PATH) == null
        strategy.resolveForeignFiles(['foo.txt': PATH]) == ['foo.txt': PATH]

        strategy.touchFile(PATH) == delegate.touchFile(PATH)
        strategy.fileStr(PATH) == delegate.fileStr(PATH)
        strategy.copyFile('foo.txt', PATH) == delegate.copyFile('foo.txt', PATH)
        strategy.exitFile(PATH) == delegate.exitFile(PATH)
        strategy.pipeInputFile(PATH) == delegate.pipeInputFile(PATH)
        strategy.getEnvScript([FOO:1],null) == delegate.getEnvScript([FOO:1],null)
    }
}

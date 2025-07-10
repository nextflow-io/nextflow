package nextflow.script

import java.nio.file.Path

import nextflow.Session
import nextflow.file.FileHelper
import nextflow.exception.ScriptRuntimeException
import spock.lang.Specification
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ParamsDslTest extends Specification {

    def 'should declare workflow params with CLI overrides'() {
        given:
        def cliParams = [input: './data', chunk_size: '3']
        def configParams = [outdir: 'results']
        def session = new Session([params: configParams + cliParams])
        session.init(null, null, cliParams, configParams)

        when:
        def dsl = new ParamsDsl()
        dsl.declare('input', Path)
        dsl.declare('chunk_size', Integer, 1)
        dsl.declare('save_intermeds', Boolean, false)
        dsl.apply(session)
        then:
        session.binding.getParams() == [input: FileHelper.asPath('./data'), chunk_size: 3, save_intermeds: false, outdir: 'results']
    }

    def 'should report error for missing required param'() {
        given:
        def cliParams = [:]
        def configParams = [outdir: 'results']
        def session = new Session()
        session.init(null, null, cliParams, configParams)

        when:
        def dsl = new ParamsDsl()
        dsl.declare('input', Path)
        dsl.declare('save_intermeds', Boolean, false)
        dsl.apply(session)
        then:
        def e = thrown(ScriptRuntimeException)
        e.message == 'Parameter `input` is required but was not specified on the command line, params file, or config'
    }

    def 'should report error for invalid param'() {
        given:
        def cliParams = [inputs: './data']
        def configParams = [outdir: 'results']
        def session = new Session()
        session.init(null, null, cliParams, configParams)

        when:
        def dsl = new ParamsDsl()
        dsl.declare('input', Path)
        dsl.declare('save_intermeds', Boolean, false)
        dsl.apply(session)
        then:
        def e = thrown(ScriptRuntimeException)
        e.message == 'Parameter `inputs` was specified on the command line or params file but is not declared in the script or config'
    }

    def 'should report error for invalid type'() {
        given:
        def cliParams = [input: './data', save_intermeds: 42]
        def configParams = [:]
        def session = new Session()
        session.init(null, null, cliParams, configParams)

        when:
        def dsl = new ParamsDsl()
        dsl.declare('input', Path)
        dsl.declare('save_intermeds', Boolean, false)
        dsl.apply(session)
        then:
        def e = thrown(ScriptRuntimeException)
        e.message == 'Parameter `save_intermeds` with type Boolean cannot be assigned to 42 [Integer]'
    }

}

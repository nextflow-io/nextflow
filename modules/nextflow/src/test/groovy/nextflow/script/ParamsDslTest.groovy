package nextflow.script

import nextflow.Session
import nextflow.exception.ScriptRuntimeException
import spock.lang.Specification
/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ParamsDslTest extends Specification {

    def 'should declare workflow params with CLI overrides'() {
        given:
        def cliParams = [input: './data']
        def config = [
            params: [outdir: 'results']
        ]
        def session = new Session(config)
        session.init(null, null, cliParams)

        when:
        def dsl = new ParamsDsl()
        dsl.declare('input')
        dsl.declare('save_intermeds', false)
        dsl.apply(session)
        then:
        session.binding.getParams() == [input: './data', save_intermeds: false, outdir: 'results']
    }

    def 'should report error for missing required param'() {
        given:
        def cliParams = [:]
        def config = [
            params: [outdir: 'results']
        ]
        def session = new Session(config)
        session.init(null, null, cliParams)

        when:
        def dsl = new ParamsDsl()
        dsl.declare('input')
        dsl.declare('save_intermeds', false)
        dsl.apply(session)
        then:
        def e = thrown(ScriptRuntimeException)
        e.message == 'Parameter `input` is required but was not specified on the command line, params file, or config'
    }

    def 'should report error for invalid param'() {
        given:
        def cliParams = [inputs: './data']
        def config = [
            params: [outdir: 'results'] + cliParams
        ]
        def session = new Session(config)
        session.init(null, null, cliParams)

        when:
        def dsl = new ParamsDsl()
        dsl.declare('input')
        dsl.declare('save_intermeds', false)
        dsl.apply(session)
        then:
        def e = thrown(ScriptRuntimeException)
        e.message == 'Parameter `inputs` was specified on the command line or params file but is not declared in the script or config'
    }

}

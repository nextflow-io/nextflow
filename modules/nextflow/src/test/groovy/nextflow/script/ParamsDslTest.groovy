package nextflow.script

import nextflow.script.types.Bag

import java.nio.file.Files
import java.nio.file.Path

import nextflow.Session
import nextflow.file.FileHelper
import nextflow.exception.ScriptRuntimeException
import spock.lang.Specification
import spock.lang.Unroll
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
        dsl.declare('input', Path, false)
        dsl.declare('chunk_size', Integer, false, 1)
        dsl.declare('save_intermeds', Boolean, false, false)
        dsl.apply(session)
        then:
        session.binding.getParams() == [input: FileHelper.asPath('./data'), chunk_size: 3, save_intermeds: false, outdir: 'results']
    }

    def 'should allow optional param'() {
        given:
        def session = new Session()
        session.init(null)

        when:
        def dsl = new ParamsDsl()
        dsl.declare('input', Path, true)
        dsl.apply(session)
        then:
        noExceptionThrown()
    }

    def 'should report error for missing required param'() {
        given:
        def cliParams = [:]
        def configParams = [outdir: 'results']
        def session = new Session()
        session.init(null, null, cliParams, configParams)

        when:
        def dsl = new ParamsDsl()
        dsl.declare('input', Path, false)
        dsl.declare('save_intermeds', Boolean, false, false)
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
        dsl.declare('input', Path, false)
        dsl.declare('save_intermeds', Boolean, false, false)
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
        dsl.declare('input', Path, false)
        dsl.declare('save_intermeds', Boolean, false, false)
        dsl.apply(session)
        then:
        def e = thrown(ScriptRuntimeException)
        e.message == 'Parameter `save_intermeds` with type Boolean cannot be assigned to 42 [Integer]'
    }

    @Unroll
    def 'should validate float param with default value'() {
        given:
        def session = new Session()
        session.init(null)

        when:
        def dsl = new ParamsDsl()
        dsl.declare('factor', Float, false, DEF_VALUE)
        dsl.apply(session)
        then:
        noExceptionThrown()

        where:
        DEF_VALUE << [ 0.1f, 0.1d, 0.1g ]
    }

    @Unroll
    def 'should validate integer param with default value'() {
        given:
        def session = new Session()
        session.init(null)

        when:
        def dsl = new ParamsDsl()
        dsl.declare('factor', Integer, false, DEF_VALUE)
        dsl.apply(session)
        then:
        noExceptionThrown()

        where:
        DEF_VALUE << [ 100i, 100l, 100g ]
    }

    def 'should load collection param from CSV file'() {
        given:
        def csvFile = Files.createTempFile('test', '.csv')
        csvFile.text = '''\
            id,name,value
            1,sample1,100
            2,sample2,200
            3,sample3,300
            '''.stripIndent()
        def cliParams = [samples: csvFile.toString()]
        def session = new Session()
        session.init(null, null, cliParams, [:])

        when:
        def dsl = new ParamsDsl()
        dsl.declare('samples', List, false)
        dsl.apply(session)

        then:
        def samples = session.binding.getParams().samples
        samples instanceof List
        samples.size() == 3
        samples[0].id == '1'
        samples[0].name == 'sample1'
        samples[0].value == '100'
        samples[1].id == '2'
        samples[2].id == '3'

        cleanup:
        csvFile?.delete()
    }

    def 'should load collection param from JSON file'() {
        given:
        def jsonFile = Files.createTempFile('test', '.json')
        jsonFile.text = '''\
            [
              {"id": 1, "name": "sample1", "value": 100},
              {"id": 2, "name": "sample2", "value": 200},
              {"id": 3, "name": "sample3", "value": 300}
            ]
            '''.stripIndent()
        def jsonFile2 = Files.createTempFile('test2', '.json')
        jsonFile2.text = '''\
            {"id": 1, "name": "sample1", "value": 100}
            '''.stripIndent()
        def cliParams = [samplesList: jsonFile.toString(),
                         samplesIter: jsonFile.toString(),
                         samplesBag: jsonFile.toString(),
                         samplesSet: jsonFile.toString(),
                         inputMap: jsonFile2.toString()]
        def session = new Session()
        session.init(null, null, cliParams, [:])

        when:
        def dsl = new ParamsDsl()
        dsl.declare('samplesList', List, false)
        dsl.declare('samplesIter', Iterable, false)
        dsl.declare('samplesBag', Bag, false)
        dsl.declare('samplesSet', Set, false)
        dsl.declare('inputMap', Map, false)
        dsl.apply(session)

        then:
        def samplesList = session.binding.getParams().samplesList
        samplesList instanceof List
        samplesList.size() == 3
        samplesList[0].id == 1
        samplesList[0].name == 'sample1'
        samplesList[0].value == 100
        samplesList[1].id == 2
        samplesList[2].id == 3

        def samplesIter = session.binding.getParams().samplesIter
        samplesIter instanceof Iterable
        samplesIter.size() == 3

        def samplesBag = session.binding.getParams().samplesBag
        samplesBag instanceof Bag
        samplesBag.size() == 3

        def samplesSet = session.binding.getParams().samplesSet
        samplesSet instanceof Set
        samplesSet.size() == 3

        def inputMap = session.binding.getParams().inputMap
        inputMap instanceof Map
        inputMap.id == 1
        inputMap.name == 'sample1'
        inputMap.value == 100


        cleanup:
        jsonFile?.delete()
        jsonFile2?.delete()
    }

    def 'should load collection param from YAML file'() {
        given:
        def yamlFile = Files.createTempFile('test', '.yml')
        yamlFile.text = '''\
            - id: 1
              name: sample1
              value: 100
            - id: 2
              name: sample2
              value: 200
            - id: 3
              name: sample3
              value: 300
            '''.stripIndent()
        def cliParams = [samples: yamlFile.toString()]
        def session = new Session()
        session.init(null, null, cliParams, [:])

        when:
        def dsl = new ParamsDsl()
        dsl.declare('samples', List, false)
        dsl.apply(session)

        then:
        def samples = session.binding.getParams().samples
        samples instanceof List
        samples.size() == 3
        samples[0].id == 1
        samples[0].name == 'sample1'
        samples[0].value == 100
        samples[1].id == 2
        samples[2].id == 3

        cleanup:
        yamlFile?.delete()
    }

    def 'should load collection param from file specified in config'() {
        given:
        def jsonFile = Files.createTempFile('test', '.json')
        jsonFile.text = '[{"x": 1}, {"x": 2}]'
        def configParams = [items: jsonFile.toString()]
        def session = new Session()
        session.init(null, null, [:], configParams)

        when:
        def dsl = new ParamsDsl()
        dsl.declare('items', List, false)
        dsl.apply(session)

        then:
        def items = session.binding.getParams().items
        items instanceof List
        items.size() == 2
        items[0].x == 1
        items[1].x == 2

        cleanup:
        jsonFile?.delete()
    }

    def 'should report error for unrecognized file format'() {
        given:
        def txtFile = Files.createTempFile('test', '.txt')
        txtFile.text = 'some text'
        def cliParams = [items: txtFile.toString()]
        def session = new Session()
        session.init(null, null, cliParams, [:])

        when:
        def dsl = new ParamsDsl()
        dsl.declare('items', List, false)
        dsl.apply(session)

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains("Unrecognized file format 'txt'")
        e.message.contains("supplied for parameter `items` -- should be CSV, JSON, or YAML")

        cleanup:
        txtFile?.delete()
    }

    def 'should report error for invalid file content type'() {
        given:
        def jsonFile = Files.createTempFile('test', '.json')
        jsonFile.text = '{"not": "a list"}'
        def cliParams = [items: jsonFile.toString()]
        def session = new Session()
        session.init(null, null, cliParams, [:])

        when:
        def dsl = new ParamsDsl()
        dsl.declare('items', List, false)
        dsl.apply(session)

        then:
        def e = thrown(ScriptRuntimeException)
        e.message.contains('Parameter `items` with type List cannot be assigned to contents of')

        cleanup:
        jsonFile?.delete()
    }

}

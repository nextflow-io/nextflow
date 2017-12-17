package nextflow.trace

import groovy.json.JsonSlurper
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ReportObserverTest extends Specification {

    def setupSpec() {
        TraceRecord.TIMEZONE = TimeZone.getTimeZone('UTC') // note: set the timezone to be sure the time string does not change on CI test servers
    }

    def 'should render json data' () {

        given:
        def now = 1429821425141
        def r1 = new TraceRecord()
        r1.task_id = '1'
        r1.name = 'foo'
        r1.process = 'alpha'

        def r2 = new TraceRecord()
        r2.task_id = '2'
        r2.name = 'bar'
        r2.submit = now
        r2.start = now + 100
        r2.complete = now + 500
        r2.realtime = 400
        r2.duration = 500
        r2.process = 'alpha'

        def r3 = new TraceRecord()
        r3.task_id = '3'
        r3.name = 'baz'
        r3.submit = now
        r3.start = now + 200
        r3.complete = now + 700
        r3.realtime = 500
        r3.duration = 700
        r3.process = 'beta'

        def observer = [:] as ReportObserver
        def data = [r1, r2, r3]

        when:
        def str = observer.renderJsonData(data)
        def json = new JsonSlurper().parseText(str)

        then:
        json.trace[0].task_id == '1'
        json.trace[0].name == 'foo'
        json.trace[0].process == 'alpha'

        json.trace[1].task_id == '2'
        json.trace[1].name == 'bar'
        json.trace[1].submit == '1429821425141'
        json.trace[1].start == '1429821425241'
        json.trace[1].complete == '1429821425641'
        json.trace[1].realtime == '400'
        json.trace[1].duration == '500'
        json.trace[1].process == 'alpha'

        json.trace[2].task_id == '3'
        json.trace[2].name == 'baz'
        json.trace[2].submit == '1429821425141'
        json.trace[2].start == '1429821425341'
        json.trace[2].complete == '1429821425841'
        json.trace[2].realtime == '500'
        json.trace[2].duration == '700'
        json.trace[2].process == 'beta'

    }

}

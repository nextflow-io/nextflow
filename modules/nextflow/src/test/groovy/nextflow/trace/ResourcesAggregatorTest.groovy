package nextflow.trace

import spock.lang.Specification

import java.util.concurrent.Executors

import groovy.json.JsonSlurper
import nextflow.Session

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ResourcesAggregatorTest extends Specification {

    private r1 = new TraceRecord([process: 'bwa-mem', name: 'bwa-mem-1',          '%cpu':1_000, peak_rss:2_000,                   realtime:3_000,  time: 5_000,  rchar:4_000, wchar:5_000 ])
    private r2 = new TraceRecord([process: 'bwa-mem', name: 'bwa-mem-2',          '%cpu':6_000, peak_rss:7_000,   memory: 10_000, realtime:8_000,  time: 10_000, rchar:9_000, wchar:10_000 ])
    private r3 = new TraceRecord([process: 'bwa-mem', name: 'bwa-mem-3', cpus: 2, '%cpu':10_000, peak_rss:12_000, memory: 10_000, realtime:13_000, time: 10_000, rchar:14_000, wchar:15_000 ])
    private r4 = new TraceRecord([process: 'multiqc', name: 'multiqc-1',          '%cpu':16_000, peak_rss:17_000, memory: 20_000, realtime:18_000, time: 20_000, rchar:19_000, wchar:20_000 ])
    private r5 = new TraceRecord([process: 'multiqc', name: 'multiqc-2', cpus: 2, '%cpu':21_000, peak_rss:22_000, memory: 20_000, realtime:23_000, time: 20_000, rchar:24_000, wchar:25_000 ])

    def 'should render summary json' () {
        given:
        def executor = Executors.newCachedThreadPool()
        def session = Mock(Session) {
            getExecService() >> executor
        }

        def observer = new ResourcesAggregator(session)
        observer.aggregate(r1)
        observer.aggregate(r2)
        observer.aggregate(r3)
        observer.aggregate(r4)
        observer.aggregate(r5)

        when:
        def map = observer.computeSummaryMap()
        then:
        map.size() == 2
        map.'bwa-mem'.cpu."min" == 1000
        map.'bwa-mem'.cpu."max" == 10000
        map.'multiqc'.mem."min" == 17000.0
        map.'multiqc'.mem."max" == 22000.0

        when:
        def json = observer.renderSummaryJson()
        def result = new JsonSlurper().parseText(json)
        then:
        result instanceof List
        result.size() == 2

        result[0].process == 'bwa-mem'
        result[0].cpu instanceof Map
        result[0].mem instanceof Map
        result[0].time instanceof Map
        result[0].reads instanceof Map
        result[0].writes instanceof Map

        result[0].time.min == 3000
        result[0].time.max == 13000
        result[0].time.q1 == 5500
        result[0].time.q2 == 8000
        result[0].time.q3 == 10500
        result[0].time.mean == 8000

        result[0].timeUsage.min == 60
        result[0].timeUsage.max == 130
        result[0].timeUsage.q1 == 70
        result[0].timeUsage.q2 == 80
        result[0].timeUsage.q3 == 105
        result[0].timeUsage.mean == 90

        result[0].cpu."min" == 1000
        result[0].cpu."minLabel" == 'bwa-mem-1'
        result[0].cpu."max" == 10000
        result[0].cpu."maxLabel" == "bwa-mem-3"
        result[0].cpu."q1" == 3500
        result[0].cpu."q2" == 6000
        result[0].cpu."q3" == 8000
        result[0].cpu."mean" == 5666.67

        result[0].cpuUsage."min" == 1000
        result[0].cpuUsage."minLabel" == 'bwa-mem-1'
        result[0].cpuUsage."max" == 6000
        result[0].cpuUsage."maxLabel" == "bwa-mem-2"
        result[0].cpuUsage."q1" == 3000
        result[0].cpuUsage."q2" == 5000
        result[0].cpuUsage."q3" == 5500
        result[0].cpuUsage."mean" == 4000

        result[1].mem."min" == 17000.0
        result[1].mem."minLabel" == "multiqc-1"
        result[1].mem."max" == 22000.0
        result[1].mem."maxLabel" == "multiqc-2"
        result[1].mem."q1" == 18250.0
        result[1].mem."q2" == 19500.0
        result[1].mem."q3" == 20750.0
        result[1].mem."mean" == 19500.0

        result[1].memUsage."min" == 85
        result[1].memUsage."minLabel" == "multiqc-1"
        result[1].memUsage."max" == 110
        result[1].memUsage."maxLabel" == "multiqc-2"
        result[1].memUsage."q1" == 91.25
        result[1].memUsage."q2" == 97.50
        result[1].memUsage."q3" == 103.75
        result[1].memUsage."mean" == 97.50

        result[1].time.min == 18000
        result[1].time.max == 23000
        result[1].time.q1 == 19250
        result[1].time.q2 == 20500
        result[1].time.q3 == 21750
        result[1].time.mean == 20500

        result[1].timeUsage.min == 90
        result[1].timeUsage.max == 115
        result[1].timeUsage.q1 == 96.25
        result[1].timeUsage.q2 == 102.50
        result[1].timeUsage.q3 == 108.75
        result[1].timeUsage.mean == 102.50

        cleanup:
        observer?.executor?.shutdown()
    }


    def 'should compute summary list' () {

        given:
        def executor = Executors.newCachedThreadPool()
        def session = Mock(Session) {
            getExecService() >> executor
        }

        def observer = new ResourcesAggregator(session)
        observer.aggregate(r1)
        observer.aggregate(r2)
        observer.aggregate(r3)
        observer.aggregate(r4)
        observer.aggregate(r5)

        when:
        def result = observer.computeSummaryList()
        then:
        result.size() == 2

        result[0].process == 'bwa-mem'
        result[0].cpu."min" == 1000
        result[0].cpu."max" == 10000
        result[0].time.min == 3000
        result[0].time.max == 13000

        result[1].process == 'multiqc'
        result[1].mem."min" == 17000.0
        result[1].mem."max" == 22000.0
        result[1].time.min == 18000
        result[1].time.max == 23000

        cleanup:
        observer?.executor?.shutdown()
    }


    def 'should maintain insertion order' () {
        given:
        def executor = Executors.newCachedThreadPool()
        def session = Mock(Session) {
            getExecService() >> executor
        }

        def observer = new ResourcesAggregator(session)
        observer.aggregate(new TraceRecord([process: 'gamma', name: 'gamma-1']))
        observer.aggregate(new TraceRecord([process: 'delta', name: 'delta-1']))
        observer.aggregate(new TraceRecord([process: 'delta', name: 'delta-2']))
        observer.aggregate(new TraceRecord([process: 'omega', name: 'omega-1']))
        observer.aggregate(new TraceRecord([process: 'gamma', name: 'gamma-2']))

        when:
        def list = observer.computeSummaryList()
        then:
        list.size() == 3
        and: // insertion order is maintained
        list[0].process == 'gamma'
        list[1].process == 'delta'
        list[2].process == 'omega'
    }

}

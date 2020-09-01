package nextflow

import spock.lang.Specification

import java.text.SimpleDateFormat

import static nextflow.extension.Bolts.DATETIME_FORMAT

import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class NextflowMetaTest extends Specification {

    private String dateToString(Date date) {
        def tz = TimeZone.getTimeZone('UTC')
        def fmt = new SimpleDateFormat(DATETIME_FORMAT)
        fmt.setTimeZone(tz)
        fmt.format(date) + ' ' + tz.getDisplayName( true, TimeZone.SHORT )
    }

    def 'should convert to map' () {
        given:
        def meta = new NextflowMeta('10.12.0', 123, Const.APP_TIMESTAMP_UTC)
        meta.enableDsl2()

        when:
        def map = meta.toJsonMap()
        then:
        map.version == '10.12.0'
        map.build == 123
        map.enable.dsl == 2
        dateToString((Date)map.timestamp) == Const.APP_TIMESTAMP_UTC

    }

    @Unroll
    def 'should find dsl2 declaration: #SCRIPT' () {
        given:
        def meta = new NextflowMeta()

        when:
        meta.checkDsl2Mode(SCRIPT)
        then:
        meta.isDsl2() == DSL2
        meta.isDsl2Final() == FINAL

        where:
        DSL2    | FINAL | SCRIPT
        false   | false | 'hello'
        false   | false | 'nextflow.preview.dsl=1'
        and:
        true    | false | 'nextflow.preview.dsl=2'
        true    | false | 'nextflow.preview.dsl = 2'
        true    | false | 'nextflow.preview.dsl =  2;'
        true    | false | '#!/bin/env nextflow\nnextflow.preview.dsl=2\nprintln foo'
        and:
        true    | true | 'nextflow.enable.dsl=2'
        true    | true | 'nextflow.enable.dsl = 2'
        true    | true | 'nextflow.enable.dsl =  2;'
        true    | true | '#!/bin/env nextflow\nnextflow.enable.dsl=2\nprintln foo'
    }
}

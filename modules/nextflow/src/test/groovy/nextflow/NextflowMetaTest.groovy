package nextflow

import spock.lang.Specification

import java.text.SimpleDateFormat

import static nextflow.extension.Bolts.DATETIME_FORMAT

import test.MockScriptRunner

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

    def dsl_eval(String str) {
        new MockScriptRunner().setScript(str).execute()
    }

    def 'should set nf meta' () {

        when:
        def ret = dsl_eval '''
        nextflow.preview.dsl=2
        
        return nextflow
        '''

        then:
        with(ret as NextflowMeta) {
            isDsl2()
            !isStrictModeEnabled()
            !isDsl2Final()
        }

    }
}

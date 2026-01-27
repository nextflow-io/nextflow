package nextflow

import static nextflow.extension.Bolts.*

import java.text.SimpleDateFormat

import spock.lang.Specification
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
        fmt.format(date) + ' ' + tz.getDisplayName(true, TimeZone.SHORT)
    }

    def 'should convert to map'() {
        given:
        def meta = new NextflowMeta('10.12.0', 123, BuildInfo.timestampUTC)

        when:
        def map = meta.toJsonMap()
        then:
        map.version == '10.12.0'
        map.build == 123
        dateToString((Date) map.timestamp) == BuildInfo.timestampUTC

    }

}

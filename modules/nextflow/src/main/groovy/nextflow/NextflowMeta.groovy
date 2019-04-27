package nextflow

import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import nextflow.util.VersionNumber

/**
 * Models nextflow script properties and metadata
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Singleton(strict = false)
@ToString(includeNames = true)
@EqualsAndHashCode
class NextflowMeta {

    static public final Closure JSON_CONVERTER = { meta, key -> key=='timestamp' ? new Date(Const.APP_TIMESTAMP) : meta.getProperty(key) }

    @Slf4j
    static class Preview {
        float dsl

        void setDsl( float num ) {
            if( num != 2 && num != 1 )
                throw new IllegalArgumentException("Not a valid DSL version number: $num")
            if( num == 2 )
                log.warn "DSL 2 IS AN EXPERIMENTAL FEATURE UNDER DEVELOPMENT -- SYNTAX MAY CHANGE IN FUTURE RELEASE"
            dsl = num
        }
    }

    final VersionNumber version
    final int build
    final String timestamp
    final Preview preview = new Preview()

    private NextflowMeta() {
        version = new VersionNumber(Const.APP_VER)
        build = Const.APP_BUILDNUM
        timestamp = Const.APP_TIMESTAMP_UTC
    }

    protected NextflowMeta(String ver, int build, String timestamp ) {
        this.version = new VersionNumber(ver)
        this.build = build
        this.timestamp = timestamp
    }

    boolean isDsl2() {
        preview.dsl == 2
    }

    void enableDsl2() {
        preview.dsl = 2
    }

    void disableDsl2() {
        preview.dsl = 1
    }
}

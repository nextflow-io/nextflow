package nextflow

import java.text.SimpleDateFormat
import java.util.regex.Pattern

import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import nextflow.util.VersionNumber
import static nextflow.extension.Bolts.DATETIME_FORMAT

/**
 * Models nextflow script properties and metadata
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Singleton(strict = false)
@ToString(includeNames = true)
@EqualsAndHashCode
class NextflowMeta {

    private static final Pattern DSL2_DECLARATION = ~/(?m)^\s*(nextflow\.(preview|enable)\.dsl\s*=\s*2)\s*;?\s*$/

    private static boolean ignoreWarnDsl2 = System.getenv('NXF_IGNORE_WARN_DSL2')=='true'

    static trait Flags {
        abstract float dsl
        abstract boolean strict
    }

    @Slf4j
    static class Preview implements Flags {
        volatile float dsl
        boolean strict

        void setDsl( float num ) {
            if( num != 2 && num != 1 )
                throw new IllegalArgumentException("Not a valid DSL version number: $num")
            if( num == 2 && !ignoreWarnDsl2 )
                log.warn1 "DSL 2 IS AN EXPERIMENTAL FEATURE UNDER DEVELOPMENT -- SYNTAX MAY CHANGE IN FUTURE RELEASE"
            dsl = num
        }
    }

    static class Features implements Flags {
        volatile float dsl
        boolean strict
    }

    final VersionNumber version
    final int build

    /*
     * Timestamp as dd-MM-yyyy HH:mm UTC formatted string
     */
    final String timestamp

    final Preview preview = new Preview()

    final Features enable = new Features()

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

    Map featuresMap() {
        final result = new LinkedHashMap()
        if( isDsl2() )
            result.dsl = 2i
        if( isStrictModeEnabled() )
            result.strict = true
        return result
    }

    Map toJsonMap() {
        final result = new LinkedHashMap<>(5)
        result.version = version.toString()
        result.build = build
        result.timestamp = parseDateStr(timestamp)
        if( isDsl2Final() ) {
            result.enable = featuresMap()
        }
        else if( isDsl2() ) {
            result.preview = featuresMap()
        }
        return result
    }

    private Date parseDateStr(String str) {
        def fmt = new SimpleDateFormat(DATETIME_FORMAT + ' Z')
        fmt.parse(str)
    }

    boolean isDsl2() {
        preview.dsl == 2 || enable.dsl == 2
    }

    boolean isDsl2Final() {
        enable.dsl == 2
    }

    void enableDsl2(boolean preview=false) {
        if( preview )
            this.preview.dsl = 2
        else
            this.enable.dsl = 2
    }

    void disableDsl2() {
        enable.dsl = 1
        preview.dsl = 1
    }

    boolean isStrictModeEnabled() {
        preview.strict || enable.strict
    }

    void checkDsl2Mode(String script) {
        final matcher = DSL2_DECLARATION.matcher(script)
        final mode = matcher.find() ? matcher.group(2) : null
        if( !mode )
            return
        if( mode == 'enable' )
            enableDsl2()
        else if( mode == 'preview' )
            enableDsl2(true)
        else
            throw new IllegalArgumentException("Unknown nextflow mode=${matcher.group(1)}")
    }
}

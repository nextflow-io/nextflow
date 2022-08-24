package nextflow.plugin

import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.util.VersionNumber
import org.pf4j.DefaultVersionManager
/**
 * Extends default version manager adding the ability to support
 * Calendar versioning
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class CustomVersionManager extends DefaultVersionManager {

    public static Pattern CAL_VER = ~/2\d\.\d\d\.\S+/

    @Override
    boolean checkVersionConstraint(String version, String constraint) {
        if( !version || !constraint || constraint=='*' || version==constraint )
            return true

        try {
            return safeCheck0(version, constraint)
        }
        catch (Throwable e) {
            log.debug "Failed to check version constraint - version: $version; constraint: $constraint"
            return false
        }
    }

    private boolean safeCheck0(String version, String constraint) {
        if( version =~ CAL_VER ) {
            if( constraint.startsWith('nextflow@') ) {
                constraint = constraint.substring('nextflow@'.size())
            }
            if( version.endsWith('-SNAPSHOT') ) {
                version = stripSuffix(version)
                constraint = stripSuffix(constraint)
            }
            return new VersionNumber(version).matches(constraint)
        }
        else {
            return super.checkVersionConstraint(version, constraint)
        }
    }

    private String stripSuffix(String str) {
        int p = str.lastIndexOf('-')
        return p==-1 ? str : str.substring(0,p)
    }

    @Override
    int compareVersions(String v1, String v2) {
        return v1=~CAL_VER || v2=~CAL_VER
                ? new VersionNumber(v1) <=> new VersionNumber(v2)
                : super.compareVersions(v1, v2)
    }

}

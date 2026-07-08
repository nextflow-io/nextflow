/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package nextflow.plugin

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.attribute.PosixFilePermissions

import nextflow.SysEnv
import nextflow.exception.AbortOperationException
import spock.lang.IgnoreIf
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Claude <noreply@anthropic.com>
 */
@IgnoreIf({ os.windows })
class PluginSecurityTest extends Specification {

    private Path tempDirWithPerms(String perms) {
        final dir = Files.createTempDirectory('test-plugin-sec')
        Files.setPosixFilePermissions(dir, PosixFilePermissions.fromString(perms))
        return dir
    }

    @Unroll
    def 'should classify trust by uid and perms - owner=#OWNER current=#CUR perms=#PERMS trusted=#TRUSTED untrusted=#UNTRUSTED' () {
        expect:
        PluginSecurity.isUntrusted(PosixFilePermissions.fromString(PERMS), OWNER, CUR, TRUSTED) == UNTRUSTED

        where:
        OWNER  | CUR    | PERMS       | TRUSTED           | UNTRUSTED
        // self-owned: private, read-only, or umask-002 group-writable -> trusted
        1000L  | 1000L  | 'rwx------' | ([] as Set)       | false
        1000L  | 1000L  | 'rwxr-xr-x' | ([] as Set)       | false
        1000L  | 1000L  | 'rwxrwxr-x' | ([] as Set)       | false   // umask 002 default - the key regression
        1000L  | 1000L  | 'rwxrwx---' | ([] as Set)       | false
        // self-owned but world-writable -> untrusted
        1000L  | 1000L  | 'rwxrwxrwx' | ([] as Set)       | true
        // root-owned admin cache -> trusted; root-owned world-writable -> untrusted
        0L     | 1000L  | 'rwxr-xr-x' | ([] as Set)       | false
        0L     | 1000L  | 'rwxrwxrwx' | ([] as Set)       | true
        // foreign owner -> untrusted even with safe perms (the PoC case)
        1201L  | 1202L  | 'rwx------' | ([] as Set)       | true
        1201L  | 1202L  | 'rwxr-xr-x' | ([] as Set)       | true
        // foreign owner listed as trusted (service-account cache) -> trusted
        1201L  | 1202L  | 'rwx------' | ([1201L] as Set)  | false
        // running as root: trust any owner, but still flag world-writable
        1201L  | 0L     | 'rwx------' | ([] as Set)       | false
        1201L  | 0L     | 'rwxrwxrwx' | ([] as Set)       | true
    }

    @Unroll
    def 'should accept a self-owned dir that is not world-writable - mode=#MODE perms=#PERMS' () {
        given:
        def dir = tempDirWithPerms(PERMS)

        when:
        PluginSecurity.checkTrustedDir(dir, MODE)
        then:
        noExceptionThrown()

        cleanup:
        dir?.deleteDir()

        where:
        MODE     | PERMS
        'strict' | 'rwx------'
        'warn'   | 'rwx------'
        'strict' | 'rwxr-xr-x'
        'strict' | 'rwxrwxr-x'   // umask-002 group-writable, owned by me
        'warn'   | 'rwxrwx---'
    }

    @Unroll
    def 'should reject a world-writable dir in strict mode - perms=#PERMS' () {
        given:
        def dir = tempDirWithPerms(PERMS)

        when:
        PluginSecurity.checkTrustedDir(dir, 'strict')
        then:
        def e = thrown(AbortOperationException)
        e.message.contains('is not secure')

        cleanup:
        dir?.deleteDir()

        where:
        PERMS << ['rwxrwxrwx', 'rwx----w-']
    }

    def 'should only warn (not throw) for a world-writable dir in warn mode' () {
        given:
        def dir = tempDirWithPerms('rwxrwxrwx')

        when:
        PluginSecurity.checkTrustedDir(dir, 'warn')
        then:
        noExceptionThrown()

        cleanup:
        dir?.deleteDir()
    }

    def 'should skip the check entirely when mode is off' () {
        given:
        def dir = tempDirWithPerms('rwxrwxrwx')

        when:
        PluginSecurity.checkTrustedDir(dir, 'off')
        then:
        noExceptionThrown()

        cleanup:
        dir?.deleteDir()
    }

    def 'should be a no-op for a null dir' () {
        when:
        PluginSecurity.checkTrustedDir(null, 'strict')
        then:
        noExceptionThrown()
    }

    def 'should default to warn mode from the environment' () {
        given:
        SysEnv.push([:])
        expect:
        PluginSecurity.getMode() == 'warn'

        cleanup:
        SysEnv.pop()
    }

    @Unroll
    def 'should resolve mode #VALUE from the environment' () {
        given:
        SysEnv.push([NXF_PLUGINS_STRICT_MODE: VALUE])
        expect:
        PluginSecurity.getMode() == VALUE

        cleanup:
        SysEnv.pop()

        where:
        VALUE << ['warn', 'strict', 'off']
    }

    @Unroll
    def 'should normalise the env value #VALUE to #EXPECTED' () {
        given:
        SysEnv.push([NXF_PLUGINS_STRICT_MODE: VALUE])
        expect:
        PluginSecurity.getMode() == EXPECTED

        cleanup:
        SysEnv.pop()

        where:
        VALUE      | EXPECTED
        'STRICT'   | 'strict'
        'Warn'     | 'warn'
        '  off  '  | 'off'
        ' Strict ' | 'strict'
    }

    @Unroll
    def 'should fall back to warn for an unrecognised or blank env value #VALUE' () {
        given:
        SysEnv.push([NXF_PLUGINS_STRICT_MODE: VALUE])
        expect:
        PluginSecurity.getMode() == 'warn'

        cleanup:
        SysEnv.pop()

        where:
        VALUE << ['bogus', 'strictly', '', '   ']
    }

    @Unroll
    def 'should parse NXF_PLUGINS_TRUSTED_OWNERS #VALUE' () {
        given:
        SysEnv.push(VALUE != null ? [NXF_PLUGINS_TRUSTED_OWNERS: VALUE] : [:])
        expect:
        PluginSecurity.getTrustedOwners() == EXPECTED

        cleanup:
        SysEnv.pop()

        where:
        VALUE              | EXPECTED
        null               | ([] as Set)
        ''                 | ([] as Set)
        '1000'             | (['1000'] as Set)
        'nextflow'         | (['nextflow'] as Set)
        '1000, nextflow ,' | (['1000', 'nextflow'] as Set)
    }

}

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
    def 'should accept a private dir owned by current user - mode=#MODE perms=#PERMS' () {
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
        'warn'   | 'rwxr-xr-x'
    }

    @Unroll
    def 'should reject a group/world-writable dir in strict mode - perms=#PERMS' () {
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
        PERMS << ['rwxrwx---', 'rwxrwxrwx', 'rwx-w----', 'rwx----w-']
    }

    @Unroll
    def 'should only warn (not throw) for a group/world-writable dir in warn mode - perms=#PERMS' () {
        given:
        def dir = tempDirWithPerms(PERMS)

        when:
        PluginSecurity.checkTrustedDir(dir, 'warn')
        then:
        noExceptionThrown()

        cleanup:
        dir?.deleteDir()

        where:
        PERMS << ['rwxrwx---', 'rwxrwxrwx']
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

}

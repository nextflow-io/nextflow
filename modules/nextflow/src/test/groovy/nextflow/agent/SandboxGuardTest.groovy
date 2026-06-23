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
package nextflow.agent

import java.nio.file.Files
import java.nio.file.Path
import spock.lang.Specification
import spock.lang.TempDir

class SandboxGuardTest extends Specification {
    @TempDir Path tmp

    def 'should allow files inside the work dir for read and write'() {
        given:
        def work = Files.createDirectories(tmp.resolve('work'))
        def f = work.resolve('out.txt')
        expect:
        SandboxGuard.isAllowed(f, work, [] as Set, true)
        SandboxGuard.isAllowed(f, work, [] as Set, false)
    }

    def 'should reject parent-traversal escape'() {
        given:
        def work = Files.createDirectories(tmp.resolve('work'))
        def escape = work.resolve('../secret.txt')
        expect:
        !SandboxGuard.isAllowed(escape, work, [] as Set, false)
        !SandboxGuard.isAllowed(escape, work, [] as Set, true)
    }

    def 'should reject an absolute path outside all roots'() {
        given:
        def work = Files.createDirectories(tmp.resolve('work'))
        def outside = Files.createDirectories(tmp.resolve('outside')).resolve('x.txt')
        expect:
        !SandboxGuard.isAllowed(outside, work, [] as Set, false)
    }

    def 'should allow reads from a whitelisted module-output dir but not writes'() {
        given:
        def work = Files.createDirectories(tmp.resolve('work'))
        def mod = Files.createDirectories(tmp.resolve('moduleout'))
        def f = mod.resolve('result.fa')
        expect:
        SandboxGuard.isAllowed(f, work, [mod] as Set, false)
        !SandboxGuard.isAllowed(f, work, [mod] as Set, true)
    }

    def 'should reject a symlink whose target escapes the sandbox'() {
        given:
        def work = Files.createDirectories(tmp.resolve('work'))
        def outside = Files.createDirectories(tmp.resolve('outside'))
        def secret = Files.write(outside.resolve('secret.txt'), 'x'.bytes)
        def link = Files.createSymbolicLink(work.resolve('link.txt'), secret)
        expect:
        !SandboxGuard.isAllowed(link, work, [] as Set, false)
    }
}

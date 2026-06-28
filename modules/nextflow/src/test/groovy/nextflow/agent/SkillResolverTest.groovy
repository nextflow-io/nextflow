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

import nextflow.exception.ScriptRuntimeException
import spock.lang.Specification
import spock.lang.TempDir

class SkillResolverTest extends Specification {

    @TempDir Path tmp

    // -- frontmatter --

    def 'should parse valid SKILL.md frontmatter + body'() {
        when:
        def r = SkillResolver.parseFrontmatter("---\nname: greet\ndescription: a greeting skill\n---\nHello body\nline2\n")
        then:
        r.name == 'greet'
        r.description == 'a greeting skill'
        r.content == 'Hello body\nline2'
    }

    def 'should tolerate BOM, CRLF and leading blank lines'() {
        when:
        def r = SkillResolver.parseFrontmatter("﻿\r\n\r\n---\r\nname: x\r\ndescription: y\r\n---\r\nBODY\r\n")
        then:
        r.name == 'x'
        r.description == 'y'
        r.content == 'BODY'
    }

    def 'should reject malformed frontmatter'() {
        when:
        SkillResolver.parseFrontmatter(text)
        then:
        thrown(ScriptRuntimeException)
        where:
        text << [
            'no frontmatter at all',                          // missing opening ---
            "---\nname: x\ndescription: y\nBODY",             // unterminated (no closing ---)
            "---\ndescription: y\n---\nbody",                 // missing name
            "---\nname: x\n---\nbody",                        // missing description
            "---\nname: x\ndescription: y\n---\n",            // empty body
            "---\nname: x\ndescription: y\n---\n   \n",       // whitespace-only body
        ]
    }

    // -- resource loading --

    def 'should load bundled resources and skip SKILL.md'() {
        given:
        def dir = Files.createDirectories(tmp.resolve('skills/foo'))
        Files.writeString(dir.resolve('SKILL.md'), "---\nname: foo\ndescription: d\n---\nbody")
        Files.createDirectories(dir.resolve('references'))
        Files.writeString(dir.resolve('references/a.txt'), 'AAA')

        when:
        def d = SkillResolver.parseSkillDir(dir)

        then:
        d.name == 'foo'
        d.resources.size() == 1
        d.resources[0].relativePath == 'references/a.txt'
        d.resources[0].content == 'AAA'
    }

    def 'should skip .git files and symlinks when loading resources'() {
        given:
        def dir = Files.createDirectories(tmp.resolve('skills/foo'))
        Files.writeString(dir.resolve('SKILL.md'), "---\nname: foo\ndescription: d\n---\nbody")
        Files.createDirectories(dir.resolve('.git'))
        Files.writeString(dir.resolve('.git/config'), 'secret')
        def outside = Files.writeString(tmp.resolve('outside.txt'), 'OUT')
        Files.createSymbolicLink(dir.resolve('link.txt'), outside)

        when:
        def d = SkillResolver.parseSkillDir(dir)

        then:
        d.resources.every { !it.relativePath.contains('.git') && it.relativePath != 'link.txt' }
    }

    def 'should cap the number of resource files'() {
        given:
        def dir = Files.createDirectories(tmp.resolve('skills/foo'))
        Files.writeString(dir.resolve('SKILL.md'), "---\nname: foo\ndescription: d\n---\nbody")
        (1..70).each { Files.writeString(dir.resolve("r${it}.txt"), 'x') }

        when:
        def d = SkillResolver.parseSkillDir(dir)

        then:
        d.resources.size() <= 64
    }

    def 'should skip an oversized resource but keep the smaller ones (sorted so the big one is first)'() {
        given:
        def dir = Files.createDirectories(tmp.resolve('skills/foo'))
        Files.writeString(dir.resolve('SKILL.md'), "---\nname: foo\ndescription: d\n---\nbody")
        // names chosen so the oversized file sorts BEFORE the small one: proves skip (continue), not break
        Files.write(dir.resolve('a-big.bin'), new byte[300 * 1024])
        Files.writeString(dir.resolve('b-small.txt'), 'keep me')

        when:
        def d = SkillResolver.parseSkillDir(dir)

        then:
        d.resources*.relativePath == ['b-small.txt']
        d.resources[0].content == 'keep me'
    }

    // -- local resolution --

    def 'should resolve a local skill by name'() {
        given:
        def base = Files.createDirectories(tmp.resolve('proj'))
        def dir = Files.createDirectories(base.resolve('skills/greet'))
        Files.writeString(dir.resolve('SKILL.md'), "---\nname: greet\ndescription: greets\n---\ninstructions")

        when:
        def list = SkillResolver.loadLocal(base.resolve('skills'), 'greet')

        then:
        list.size() == 1
        list[0].name == 'greet'
        list[0].content == 'instructions'
    }

    def 'should error when a local skill dir is missing'() {
        given:
        def base = Files.createDirectories(tmp.resolve('proj'))
        when:
        SkillResolver.loadLocal(base.resolve('skills'), 'nope')
        then:
        thrown(ScriptRuntimeException)
    }

    def 'should error when a skill dir has no SKILL.md'() {
        given:
        def base = Files.createDirectories(tmp.resolve('proj'))
        Files.createDirectories(base.resolve('skills/empty'))
        when:
        SkillResolver.loadLocal(base.resolve('skills'), 'empty')
        then:
        thrown(ScriptRuntimeException)
    }

    // -- remote reference disambiguation --

    def 'should recognize explicit github refs as remote and bare names/registry refs as local'() {
        expect:
        SkillResolver.isRemoteRef('https://github.com/org/repo')
        SkillResolver.isRemoteRef('github.com/org/repo@v1')
        SkillResolver.isRemoteRef('git@github.com:org/repo')
        and: 'a bare name and a bare org/repo (registry-style) are NOT remote'
        !SkillResolver.isRemoteRef('greet')
        !SkillResolver.isRemoteRef('nf-core/fastqc')
        !SkillResolver.isRemoteRef(null)
    }

    def 'should parse a remote ref into clone url, repo and rev'() {
        when:
        def p = SkillResolver.parseRemoteRef('github.com/org/myrepo@abc123')
        then:
        p.url == 'https://github.com/org/myrepo.git'
        p.repo == 'myrepo'
        p.rev == 'abc123'
    }

    // -- remote fetch (offline, via a file:// clone) --

    def 'should fetch a remote skill via a file:// clone and reuse the cache'() {
        given: 'a local git repo acting as the remote, containing a skill'
        def remote = Files.createDirectories(tmp.resolve('remote'))
        def git = org.eclipse.jgit.api.Git.init().setDirectory(remote.toFile()).call()
        Files.writeString(remote.resolve('SKILL.md'), "---\nname: remoteskill\ndescription: from git\n---\nremote instructions")
        git.add().addFilepattern('.').call()
        git.commit().setMessage('init').setSign(false).setAuthor('t','t@t').setCommitter('t','t@t').call()
        git.close()
        def base = Files.createDirectories(tmp.resolve('proj'))
        def url = "file://${remote.toAbsolutePath()}/.git".toString()

        when:
        def list = SkillResolver.loadRemoteUrl(base.resolve('skills'), url, 'myrepo', null)

        then:
        list.size() == 1
        list[0].name == 'remoteskill'
        list[0].content == 'remote instructions'
        Files.isDirectory(base.resolve('skills/myrepo'))

        when: 'called again with a bogus url, the existing cache dir is reused (no re-clone)'
        def list2 = SkillResolver.loadRemoteUrl(base.resolve('skills'), 'file:///does/not/exist.git', 'myrepo', null)

        then:
        list2.size() == 1
        list2[0].name == 'remoteskill'
    }

    def 'should sanitize a slash-bearing rev into a flat cache dir name (no path traversal/nesting)'() {
        given: 'a remote repo with a branch whose name contains a slash'
        def remote = Files.createDirectories(tmp.resolve('remote'))
        def git = org.eclipse.jgit.api.Git.init().setDirectory(remote.toFile()).call()
        Files.writeString(remote.resolve('SKILL.md'), "---\nname: branchskill\ndescription: from a branch\n---\nbranch instructions")
        git.add().addFilepattern('.').call()
        git.commit().setMessage('init').setSign(false).setAuthor('t','t@t').setCommitter('t','t@t').call()
        git.branchCreate().setName('feature/x').call()
        git.close()
        def base = Files.createDirectories(tmp.resolve('proj'))
        def url = "file://${remote.toAbsolutePath()}/.git".toString()
        def skillsRoot = base.resolve('skills')

        when:
        def list = SkillResolver.loadRemoteUrl(skillsRoot, url, 'repo', 'feature/x')

        then: 'the cache dir is a single flat segment with the slash sanitized'
        list[0].name == 'branchskill'
        Files.isDirectory(skillsRoot.resolve('repo@feature_x'))
        and: 'no nested repo@feature/x path was created'
        !Files.exists(skillsRoot.resolve('repo@feature'))
    }
}

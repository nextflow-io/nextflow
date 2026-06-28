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

import java.nio.charset.StandardCharsets
import java.nio.file.Files
import java.nio.file.Path
import java.util.stream.Stream

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.ScriptRuntimeException
import org.yaml.snakeyaml.Yaml

/**
 * Resolves an agent {@code skills} directive entry into one or more portable
 * {@link SkillDescriptor}s. An entry is either a <b>local</b> skill (a
 * {@code <baseDir>/skills/<name>/} directory containing a {@code SKILL.md}) or a
 * <b>remote</b> GitHub reference (cloned and cached into that same {@code skills/}
 * directory — see {@code loadRemote}).
 *
 * <p>All filesystem and SCM work lives here in core; the {@code nf-agent} plugin only
 * receives the portable descriptors. {@code SKILL.md} is parsed with a hand-rolled
 * YAML-frontmatter split (snakeyaml has no frontmatter support): the leading
 * {@code ---}…{@code ---} block is the metadata ({@code name}, {@code description}),
 * the remainder is the skill {@code content}. Bundled files (other than {@code SKILL.md})
 * are loaded as resources, skipping {@code .git/} and symlinks, rejecting path escapes,
 * and capped for safety.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class SkillResolver {

    static final String SKILL_FILE = 'SKILL.md'
    static final String SKILLS_DIR = 'skills'

    private static final int MAX_RESOURCE_FILES = 64
    private static final long MAX_RESOURCE_BYTES = 256 * 1024

    /**
     * Parse {@code SKILL.md} text into {@code [name, description, content]}. Tolerates a
     * leading BOM, CRLF line endings and leading blank lines. Throws when the frontmatter
     * is missing/unterminated or {@code name}/{@code description} are absent.
     */
    static Map parseFrontmatter(String raw) {
        if( raw == null )
            throw new ScriptRuntimeException("Invalid SKILL.md: empty content")
        String text = raw
        if( text.startsWith('﻿') )
            text = text.substring(1)
        text = text.replace('\r\n', '\n').replace('\r', '\n')
        final List<String> lines = text.split('\n', -1) as List<String>
        int i = 0
        while( i < lines.size() && lines[i].trim().isEmpty() )
            i++
        if( i >= lines.size() || lines[i].trim() != '---' )
            throw new ScriptRuntimeException("Invalid SKILL.md: missing YAML frontmatter (expected a leading '---' line)")
        i++
        final StringBuilder fm = new StringBuilder()
        boolean closed = false
        while( i < lines.size() ) {
            if( lines[i].trim() == '---' ) { closed = true; i++; break }
            fm.append(lines[i]).append('\n')
            i++
        }
        if( !closed )
            throw new ScriptRuntimeException("Invalid SKILL.md: unterminated YAML frontmatter (missing closing '---')")
        final Object loaded = new Yaml().load(fm.toString())
        final Map data = loaded instanceof Map ? (Map) loaded : [:]
        final String name = data.get('name')?.toString()?.trim()
        final String description = data.get('description')?.toString()?.trim()
        if( !name )
            throw new ScriptRuntimeException("Invalid SKILL.md: missing 'name' in frontmatter")
        if( !description )
            throw new ScriptRuntimeException("Invalid SKILL.md: missing 'description' in frontmatter")
        final String content = (i < lines.size() ? lines.subList(i, lines.size()).join('\n') : '').trim()
        return [name: name, description: description, content: content]
    }

    /**
     * Parse a single skill directory (one containing a {@code SKILL.md}) into a descriptor.
     */
    static SkillDescriptor parseSkillDir(Path skillDir) {
        final Path md = skillDir.resolve(SKILL_FILE)
        if( !Files.exists(md) )
            throw new ScriptRuntimeException("Skill directory `${skillDir}` has no ${SKILL_FILE}")
        final Map meta = parseFrontmatter(new String(Files.readAllBytes(md), StandardCharsets.UTF_8))
        final List<SkillResource> resources = loadResources(skillDir)
        return new SkillDescriptor(meta.name as String, meta.description as String, meta.content as String, resources)
    }

    /**
     * Load the bundled resource files under a skill directory: every regular, non-symlink
     * file other than the top-level {@code SKILL.md}, excluding anything under a {@code .git}
     * directory, rejecting paths that escape the skill dir, capped at {@value #MAX_RESOURCE_FILES}
     * files / {@value #MAX_RESOURCE_BYTES} bytes total.
     */
    static List<SkillResource> loadResources(Path skillDir) {
        final Path root = skillDir.toRealPath()
        final List<SkillResource> result = new ArrayList<>()
        long totalBytes = 0
        final Stream<Path> walk = Files.walk(skillDir)
        try {
            final Iterator<Path> it = walk.iterator()
            while( it.hasNext() ) {
                final Path p = it.next()
                if( !Files.isRegularFile(p) || Files.isSymbolicLink(p) )
                    continue
                if( p.fileName?.toString() == SKILL_FILE && skillDir == p.parent )
                    continue
                if( hasGitSegment(skillDir.relativize(p)) )
                    continue
                if( !p.toRealPath().startsWith(root) )
                    continue
                if( result.size() >= MAX_RESOURCE_FILES ) {
                    log.warn("Skill `${skillDir.fileName}`: more than ${MAX_RESOURCE_FILES} resource files - ignoring the rest")
                    break
                }
                final byte[] bytes = Files.readAllBytes(p)
                if( totalBytes + bytes.length > MAX_RESOURCE_BYTES ) {
                    log.warn("Skill `${skillDir.fileName}`: resource files exceed ${MAX_RESOURCE_BYTES} bytes - ignoring `${skillDir.relativize(p)}` and the rest")
                    break
                }
                totalBytes += bytes.length
                result.add(new SkillResource(skillDir.relativize(p).toString(), new String(bytes, StandardCharsets.UTF_8)))
            }
        }
        finally {
            walk.close()
        }
        return result
    }

    private static boolean hasGitSegment(Path relative) {
        for( final Path seg : relative ) {
            if( seg.toString() == '.git' )
                return true
        }
        return false
    }

    /**
     * Resolve a local skill by name: {@code <baseDir>/skills/<name>/}. The directory itself
     * may be a single skill (has a {@code SKILL.md}) or hold multiple skills in subdirectories.
     */
    static List<SkillDescriptor> loadLocal(Path baseDir, String name) {
        final Path dir = baseDir.resolve(SKILLS_DIR).resolve(name)
        if( !Files.isDirectory(dir) )
            throw new ScriptRuntimeException("Agent skill `${name}` not found: no directory `${dir}`")
        return scanSkillRoot(dir, name)
    }

    /**
     * Treat {@code dir} as a single skill when it has a {@code SKILL.md}, otherwise scan its
     * immediate subdirectories for skills. Throws when no {@code SKILL.md} is found anywhere.
     */
    static List<SkillDescriptor> scanSkillRoot(Path dir, String label) {
        if( Files.exists(dir.resolve(SKILL_FILE)) )
            return [ parseSkillDir(dir) ]
        final List<SkillDescriptor> result = new ArrayList<>()
        final Stream<Path> kids = Files.list(dir)
        try {
            final Iterator<Path> it = kids.iterator()
            while( it.hasNext() ) {
                final Path sub = it.next()
                if( Files.isDirectory(sub) && Files.exists(sub.resolve(SKILL_FILE)) )
                    result.add(parseSkillDir(sub))
            }
        }
        finally {
            kids.close()
        }
        if( result.isEmpty() )
            throw new ScriptRuntimeException("Agent skill `${label}` has no ${SKILL_FILE} (looked in `${dir}` and its subdirectories)")
        return result
    }
}

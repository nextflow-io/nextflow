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
import java.nio.file.DirectoryNotEmptyException
import java.nio.file.FileAlreadyExistsException
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.StandardCopyOption
import java.util.regex.Matcher
import java.util.stream.Stream

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.ScriptRuntimeException
import org.eclipse.jgit.api.Git
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

    /** The {@code ---} delimiter line that fences the YAML frontmatter block. */
    private static final String FENCE = '---'

    private static final int MAX_RESOURCE_FILES = 64
    private static final long MAX_RESOURCE_BYTES = 256 * 1024

    /**
     * Parse {@code SKILL.md} text into {@code [name, description, content]}: normalize the text,
     * locate the {@code ---}…{@code ---} frontmatter fences, then extract and validate the fields.
     * Tolerates a leading BOM, CRLF line endings and leading blank lines. Throws when the frontmatter
     * is missing/unterminated or {@code name}/{@code description}/body are absent.
     */
    static Map parseFrontmatter(String raw) {
        final List<String> lines = normalizedLines(raw)
        final int open = skipBlankLines(lines, 0)
        if( open >= lines.size() || lines[open].trim() != FENCE )
            throw new ScriptRuntimeException("Invalid SKILL.md: missing YAML frontmatter (expected a leading '${FENCE}' line)")
        final int close = indexOfFence(lines, open + 1)
        if( close < 0 )
            throw new ScriptRuntimeException("Invalid SKILL.md: unterminated YAML frontmatter (missing closing '${FENCE}')")
        final String yaml = lines.subList(open + 1, close).join('\n')
        final String body = joinFrom(lines, close + 1).trim()
        return frontmatterFields(yaml, body)
    }

    /** Strip a leading BOM, normalize CRLF/CR to LF, and split into lines. */
    private static List<String> normalizedLines(String raw) {
        if( raw == null )
            throw new ScriptRuntimeException("Invalid SKILL.md: empty content")
        final String text = (raw.startsWith('﻿') ? raw.substring(1) : raw)
            .replace('\r\n', '\n').replace('\r', '\n')
        return text.split('\n', -1) as List<String>
    }

    /** Index of the first non-blank line at or after {@code from}. */
    private static int skipBlankLines(List<String> lines, int from) {
        int i = from
        while( i < lines.size() && lines[i].trim().isEmpty() )
            i++
        return i
    }

    /** Index of the next {@code ---} fence line at or after {@code from}, or {@code -1} if none. */
    private static int indexOfFence(List<String> lines, int from) {
        for( int i = from; i < lines.size(); i++ )
            if( lines[i].trim() == FENCE )
                return i
        return -1
    }

    /** Join lines from {@code from} to the end with LF, or empty when {@code from} is past the end. */
    private static String joinFrom(List<String> lines, int from) {
        return from < lines.size() ? lines.subList(from, lines.size()).join('\n') : ''
    }

    /**
     * Extract and validate the frontmatter fields. {@code name}/{@code description} come from the YAML
     * block; the {@code body} is the skill content. Validating the body here (core, langchain4j-free)
     * makes an empty SKILL.md fail with a clear Nextflow error naming the skill, rather than a raw
     * langchain4j {@code IllegalArgumentException} later.
     */
    private static Map frontmatterFields(String yaml, String body) {
        final Object loaded = new Yaml().load(yaml)
        final Map data = loaded instanceof Map ? (Map) loaded : [:]
        final String name = data.get('name')?.toString()?.trim()
        final String description = data.get('description')?.toString()?.trim()
        if( !name )
            throw new ScriptRuntimeException("Invalid SKILL.md: missing 'name' in frontmatter")
        if( !description )
            throw new ScriptRuntimeException("Invalid SKILL.md: missing 'description' in frontmatter")
        if( !body )
            throw new ScriptRuntimeException("Invalid SKILL.md: skill `${name}` has an empty body (no instructions)")
        return [name: name, description: description, content: body]
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
        return readUnderCaps(skillDir, eligibleResourceFiles(skillDir))
    }

    /**
     * The bundled files eligible to be skill resources, in deterministic (lexicographic-by-relative-path)
     * order so the same skill yields the same resource set on every host. Excludes the top-level
     * {@code SKILL.md}, symlinks, anything under {@code .git/}, and any path escaping the skill dir.
     */
    private static List<Path> eligibleResourceFiles(Path skillDir) {
        final Path root = skillDir.toRealPath()
        final List<Path> files = new ArrayList<>()
        final Stream<Path> walk = Files.walk(skillDir)
        try {
            final Iterator<Path> it = walk.iterator()
            while( it.hasNext() ) {
                final Path p = it.next()
                if( isEligibleResource(p, skillDir, root) )
                    files.add(p)
            }
        }
        finally {
            walk.close()
        }
        files.sort(byRelativePath(skillDir))
        return files
    }

    private static boolean isEligibleResource(Path p, Path skillDir, Path root) {
        if( !Files.isRegularFile(p) || Files.isSymbolicLink(p) )
            return false
        if( p.fileName?.toString() == SKILL_FILE && skillDir == p.parent )
            return false
        if( hasGitSegment(skillDir.relativize(p)) )
            return false
        return p.toRealPath().startsWith(root)
    }

    private static Comparator<Path> byRelativePath(Path skillDir) {
        return new Comparator<Path>() {
            int compare(Path a, Path b) {
                return skillDir.relativize(a).toString() <=> skillDir.relativize(b).toString()
            }
        }
    }

    /**
     * Read the given files into resources, enforcing the file-count and total-byte caps. Each file's
     * size is stat'd BEFORE reading so an oversized file is never loaded into memory (DoS-safe); an
     * over-budget file is skipped (not a hard stop) so smaller later files are still included.
     */
    private static List<SkillResource> readUnderCaps(Path skillDir, List<Path> files) {
        final List<SkillResource> result = new ArrayList<>()
        long totalBytes = 0
        for( final Path p : files ) {
            if( result.size() >= MAX_RESOURCE_FILES ) {
                log.warn("Skill `${skillDir.fileName}`: more than ${MAX_RESOURCE_FILES} resource files - ignoring the rest")
                break
            }
            final long size = Files.size(p)
            if( totalBytes + size > MAX_RESOURCE_BYTES ) {
                log.warn("Skill `${skillDir.fileName}`: skipping resource `${skillDir.relativize(p)}` - would exceed the ${MAX_RESOURCE_BYTES}-byte resource budget")
                continue
            }
            totalBytes += size
            result.add(new SkillResource(skillDir.relativize(p).toString(), new String(Files.readAllBytes(p), StandardCharsets.UTF_8)))
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
     * Resolve a local skill by name under the given skills-root directory ({@code <skillsRoot>/<name>/}).
     * That directory itself may be a single skill (has a {@code SKILL.md}) or hold multiple skills in
     * subdirectories.
     */
    static List<SkillDescriptor> loadLocal(Path skillsRoot, String name) {
        final Path dir = skillsRoot.resolve(name)
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

    // -- remote (GitHub) resolution --

    private static final java.util.regex.Pattern REMOTE_REF =
        ~/^(https:\/\/github\.com\/|git@github\.com:|github\.com\/)([^\/@:]+)\/([^\/@:]+?)(?:@(.+))?$/

    /**
     * Whether a {@code skills} entry is a remote GitHub reference. Only explicit GitHub
     * forms qualify ({@code https://github.com/<org>/<repo>}, {@code git@github.com:<org>/<repo>},
     * {@code github.com/<org>/<repo>}, each optionally {@code @rev}); a bare {@code <org>/<repo>}
     * is NOT remote (that shape is reserved for registry-style module refs) — anything that is
     * not remote is treated as a local skill name.
     */
    static boolean isRemoteRef(String ref) {
        if( !ref )
            return false
        return REMOTE_REF.matcher(ref).matches()
    }

    /** Parse a remote GitHub ref into {@code [url, repo, rev]} (rev may be null). */
    static Map parseRemoteRef(String ref) {
        final Matcher m = REMOTE_REF.matcher(ref)
        if( !m.matches() )
            throw new ScriptRuntimeException("Invalid remote skill reference `${ref}` - expected github.com/<org>/<repo>[@rev]")
        final String prefix = m.group(1)
        final String org = m.group(2)
        String repo = m.group(3)
        if( repo.endsWith('.git') )
            repo = repo.substring(0, repo.length() - 4)
        final String rev = m.group(4)
        final String url = prefix == 'git@github.com:'
            ? "git@github.com:${org}/${repo}.git".toString()
            : "https://github.com/${org}/${repo}.git".toString()
        return [url: url, repo: repo, rev: rev]
    }

    /**
     * Resolve a remote GitHub skill reference: clone (and cache) the repo into the given
     * skills-root directory, then load the skill(s) it contains.
     */
    static List<SkillDescriptor> loadRemote(Path skillsRoot, String ref) {
        final Map parsed = parseRemoteRef(ref)
        if( !parsed.rev )
            log.warn("Agent skill `${ref}` is not pinned to a commit - its instructions may change between runs; pin a commit SHA (e.g. `${ref}@<sha>`) for reproducibility")
        return loadRemoteUrl(skillsRoot, parsed.url as String, parsed.repo as String, parsed.rev as String)
    }

    /**
     * Low-level remote fetch: clone {@code cloneUrl} (checking out {@code rev} when given) into a
     * rev-keyed cache dir under {@code skillsRoot}, reusing it when present, then scan it for skills.
     * The cache dir is keyed by repo name and {@code @rev} so a pinned revision is never silently
     * served from a cache populated for a different revision.
     */
    static List<SkillDescriptor> loadRemoteUrl(Path skillsRoot, String cloneUrl, String repoName, String rev) {
        final Path cacheDir = cacheDirFor(skillsRoot, repoName, rev)
        if( !Files.isDirectory(cacheDir) )
            cloneInto(cloneUrl, rev, cacheDir)
        return scanSkillRoot(cacheDir, repoName)
    }

    /**
     * The rev-keyed cache dir under {@code skillsRoot}, so a pinned revision is never served from a
     * cache populated for a different one. The rev is sanitized for the on-disk NAME only (it may
     * contain {@code /} or {@code ..}, e.g. {@code feature/x}, {@code refs/tags/v1}) — keeping the
     * cache a single flat segment and preventing path traversal out of {@code skillsRoot}; the real
     * rev still drives the git checkout.
     */
    private static Path cacheDirFor(Path skillsRoot, String repoName, String rev) {
        final String safeRev = rev ? rev.replaceAll(/[^A-Za-z0-9._-]/, '_') : null
        final String cacheName = safeRev ? "${repoName}@${safeRev}".toString() : repoName
        final Path cacheDir = skillsRoot.resolve(cacheName).normalize()
        if( cacheDir.parent == null || !cacheDir.startsWith(skillsRoot.normalize()) )
            throw new ScriptRuntimeException("Invalid skills cache path `${cacheDir}` for skill repo `${repoName}`")
        return cacheDir
    }

    /**
     * Clone into a sibling temp dir then atomically rename into {@code cacheDir}, so an
     * interrupted/failed clone never leaves a half-populated cache. A shallow clone is used only
     * for a moving (default-branch) remote ref; a pinned {@code rev} or a local {@code file://}
     * source is full-cloned so the requested commit is reachable.
     */
    private static void cloneInto(String cloneUrl, String rev, Path cacheDir) {
        Files.createDirectories(cacheDir.parent)
        final Path tmp = Files.createTempDirectory(cacheDir.parent, '.skill-clone-')
        Files.delete(tmp)  // JGit creates the directory itself
        try {
            cloneAndCheckout(cloneUrl, rev, tmp)
            publishClone(tmp, cacheDir)
        }
        catch( Exception e ) {
            deleteQuietly(tmp)
            throw new ScriptRuntimeException("Unable to fetch remote skill from `${redactUrl(cloneUrl)}`${rev ? " (@${rev})" : ''}: ${e.message}", e)
        }
    }

    /**
     * Clone {@code cloneUrl} into {@code target} and check out {@code rev} when given. A shallow clone
     * is used only for a moving (default-branch) remote ref; a pinned {@code rev} or a {@code file://}
     * source is full-cloned so the requested commit is reachable.
     */
    private static void cloneAndCheckout(String cloneUrl, String rev, Path target) {
        final clone = Git.cloneRepository().setURI(cloneUrl).setDirectory(target.toFile())
        if( !rev && !cloneUrl.startsWith('file:') )
            clone.setDepth(1)
        final Git git = clone.call()
        try {
            if( rev )
                checkoutRev(git, rev)
        }
        finally {
            git.close()
        }
    }

    /**
     * Check out a rev so a SHA, a tag, OR a branch all work: a fresh clone exposes branches only as
     * remote-tracking refs ({@code origin/<branch>}), so a bare {@code checkout(branch)} would fail —
     * resolve the rev to a concrete {@code ObjectId} and check that out.
     */
    private static void checkoutRev(Git git, String rev) {
        org.eclipse.jgit.lib.ObjectId id = git.repository.resolve(rev)
        if( id == null )
            id = git.repository.resolve("origin/${rev}".toString())
        if( id == null )
            throw new ScriptRuntimeException("Revision `${rev}` not found in the remote repository")
        git.checkout().setName(id.name()).call()
    }

    /**
     * Atomically rename the temp clone into the cache dir. If a concurrent run populated the cache
     * first, discard our clone and use theirs (the atomic move guarantees no half-populated cache).
     */
    private static void publishClone(Path tmp, Path cacheDir) {
        try {
            Files.move(tmp, cacheDir, StandardCopyOption.ATOMIC_MOVE)
        }
        catch( FileAlreadyExistsException | DirectoryNotEmptyException raced ) {
            deleteQuietly(tmp)
        }
    }

    private static void deleteQuietly(Path dir) {
        try { dir.toFile().deleteDir() } catch( Exception ignore ) {}
    }

    /** Redact any {@code user:token@} userinfo from a URL before it appears in an error/log message. */
    private static String redactUrl(String url) {
        return url?.replaceAll('//[^@/]+@', '//***@')
    }
}

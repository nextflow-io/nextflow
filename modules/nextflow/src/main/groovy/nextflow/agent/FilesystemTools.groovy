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
import java.nio.file.FileSystems
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.PathMatcher
import java.util.regex.Pattern
import java.util.regex.PatternSyntaxException
import java.util.stream.Stream

import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.transform.TupleConstructor
import groovy.util.logging.Slf4j

/**
 * The {@code fs:} tool family: the six separately-named tools {@code read}, {@code write},
 * {@code edit}, {@code ls}, {@code grep} and {@code find} — both their wire-level definition
 * (the descriptors the model reads) and their driver-JVM implementation.
 *
 * <p>These are the Pi-baseline names, deliberately: a tool ref is a contract and the wire
 * name is the coordination point, so the model sees the same six names whichever runner
 * executes them — the runner's own builtins in a container, or the Groovy implementation
 * below in the driver JVM. There is no {@code exists} tool: it had no
 * counterpart on the other runner and {@code ls}/{@code read} cover it.
 *
 * <p>A tool DESCRIPTION is the only documentation the model ever reads, so each one states
 * the sandbox boundary explicitly — a model that does not know where the boundary is spends
 * turns probing it. The caps below are quoted verbatim in the descriptions for the same
 * reason: the model must be able to predict a truncated result rather than discover it.
 *
 * <p>The work dir is bound per call from the {@link DispatchContext} and is never supplied
 * by the LLM; every argument here is a path relative to it (or an absolute path that must
 * still fall inside the sandbox).
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class FilesystemTools {

    static final String READ = 'read'
    static final String WRITE = 'write'
    static final String EDIT = 'edit'
    static final String LS = 'ls'
    static final String GREP = 'grep'
    static final String FIND = 'find'

    /**
     * The six wire names, in canonical order. Sourced from the grammar inventory rather than
     * restated here so the set a directive can SELECT and the set the bridge can SERVE cannot
     * drift apart; {@link #descriptor} throws for a name it does not know, which turns any
     * future drift into a build-time failure rather than a tool the model can call but nothing
     * implements.
     */
    static final List<String> NAMES = ToolRefResolver.FS_TOOLS

    /** Default and hard cap on the number of results {@code grep}/{@code find} return. */
    static final int DEFAULT_MAX_RESULTS = 200
    static final int MAX_MAX_RESULTS = 1000

    /** Default and hard cap on the directory depth {@code grep}/{@code find} walk. */
    static final int DEFAULT_MAX_DEPTH = 20
    static final int MAX_MAX_DEPTH = 50

    /** Files larger than this are skipped by {@code grep}: a match in a bulk data file is noise. */
    static final long MAX_GREP_FILE_BYTES = 2 * 1024 * 1024L

    /**
     * Hard cap on the number of directory entries one {@code grep}/{@code find} call may VISIT.
     *
     * <p>{@code max_results} bounds the answer, not the work: a pattern that matches nothing walks
     * the whole tree, and every visited entry costs a {@code stat} plus a {@code toRealPath()}
     * sandbox check — on the agent dispatch thread, with no deadline. This bounds the walk itself,
     * and tripping it sets {@code truncated} with a {@code truncated_reason} that distinguishes
     * "you saw the first N matches" from "the search stopped early", which are opposite
     * instructions for the model: narrow the pattern versus narrow the search root.
     */
    static final int MAX_VISITED_ENTRIES = 20_000

    /** Hard cap on the number of files one {@code grep} call opens and reads; see above. */
    static final int MAX_GREP_FILES = 2_000

    /** {@code truncated_reason} when the RESULT limit ({@code max_results}) stopped the search. */
    static final String TRUNCATED_RESULTS = 'max_results'

    /** {@code truncated_reason} when the SEARCH budget stopped it: the result set is incomplete. */
    static final String TRUNCATED_SEARCH = 'search_budget'

    /** Matched lines are truncated to this many characters before being returned. */
    static final int MAX_LINE_CHARS = 300

    /**
     * The sandbox sentence shared by every description. Repeated in each tool rather than stated
     * once in the system prompt because a tool description is the only text guaranteed to travel
     * with the tool through every runner and every provider.
     */
    private static final String SANDBOX =
        'Confined to the agent sandbox: the task work dir plus any module-output files returned by an earlier tool call. A path outside it is refused.'

    /**
     * The descriptor of one {@code fs:} tool by its wire name.
     *
     * @throws IllegalArgumentException for a name outside {@link #NAMES}
     */
    static ToolDescriptor descriptor(String name) {
        switch( name ) {
            case READ:  return readDescriptor()
            case WRITE: return writeDescriptor()
            case EDIT:  return editDescriptor()
            case LS:    return lsDescriptor()
            case GREP:  return grepDescriptor()
            case FIND:  return findDescriptor()
            default:
                throw new IllegalArgumentException("Unknown filesystem tool `${name}` - known tools: ${NAMES}")
        }
    }

    /**
     * The descriptors of the given wire names, in the order given. An EMPTY collection means
     * exactly that — no fs: tool was selected — never the whole family; the only caller
     * ({@code ModuleToolBridge.filesystemDescriptors}) holds a non-null set the constructor
     * already defaulted, so there is no null case to absorb here.
     */
    static List<ToolDescriptor> descriptors(Collection<String> names) {
        return names.collect { descriptor(it) }
    }

    // -------------------------------------------------------------------------
    // the six descriptors
    // -------------------------------------------------------------------------

    private static ToolDescriptor readDescriptor() {
        final input = ToolSchema.object(
            [
                path: [type: 'string', description: 'File to read, relative to the agent work dir (or an absolute path inside the sandbox).'],
            ] as Map<String,Object>,
            ['path'] )
        final desc = 'Read the contents of a single file. ' +
            'Small text-like files (.txt .md .json .yaml .yml .csv .tsv .tab .log) are returned inline under `content`; ' +
            'every other file — binary, bulk data, or over the inline size limit — is returned as an opaque absolute `path` handle you can pass to another tool but cannot see the bytes of. ' +
            "Use `${LS}` first if you are unsure the file exists. " + SANDBOX
        return new ToolDescriptor(READ, desc, input, null)
    }

    private static ToolDescriptor writeDescriptor() {
        final input = ToolSchema.object(
            [
                path   : [type: 'string', description: 'File to write, relative to the agent work dir. Missing parent directories are created.'],
                content: [type: 'string', description: 'The full new contents of the file. The file is overwritten, not appended to.'],
            ] as Map<String,Object>,
            ['path','content'] )
        final desc = 'Create a file or replace its entire contents. Missing parent directories are created. ' +
            "To change part of an existing file use `${EDIT}` instead — this tool overwrites everything. " +
            'Writes are confined to the agent work dir only: unlike reads, a module-output path outside the work dir is NOT writable.'
        return new ToolDescriptor(WRITE, desc, input, null)
    }

    private static ToolDescriptor editDescriptor() {
        final input = ToolSchema.object(
            [
                path       : [type: 'string', description: 'File to edit, relative to the agent work dir.'],
                old_string : [type: 'string', description: 'The exact text to replace, including whitespace and indentation. Must occur exactly once in the file unless `replace_all` is true.'],
                new_string : [type: 'string', description: 'The text to put in its place. May be empty to delete the matched text.'],
                replace_all: [type: 'boolean', description: 'Replace every occurrence instead of requiring a unique one. Defaults to false.'],
            ] as Map<String,Object>,
            ['path','old_string','new_string'] )
        final desc = 'Replace an exact literal string in an existing file, leaving the rest untouched. ' +
            '`old_string` is matched literally, never as a regular expression. ' +
            'If it occurs more than once the edit is REFUSED rather than applied to the first match: either extend `old_string` with surrounding lines until it is unique, or pass `replace_all: true` when you really mean every occurrence. ' +
            'The number of replacements made is reported back. ' +
            'Edits are confined to the agent work dir only, like writes.'
        return new ToolDescriptor(EDIT, desc, input, null)
    }

    private static ToolDescriptor lsDescriptor() {
        final input = ToolSchema.object(
            [
                path: [type: 'string', description: 'Directory to list, relative to the agent work dir. Defaults to the work dir itself.'],
            ] as Map<String,Object>,
            [] )
        final desc = 'List the immediate entries of a directory (not recursive). ' +
            'Each entry reports its `name`, its `type` (`file` or `dir`) and, for a file, its `size` in bytes — so you can tell whether reading it is worthwhile. ' +
            "Use `${FIND}` to search a directory tree instead. " + SANDBOX
        return new ToolDescriptor(LS, desc, input, null)
    }

    private static ToolDescriptor grepDescriptor() {
        final input = ToolSchema.object(
            [
                pattern         : [type: 'string', description: 'Regular expression (Java/PCRE syntax) matched against each line. A plain substring is a valid pattern.'],
                path            : [type: 'string', description: 'File or directory to search, relative to the agent work dir. Defaults to the work dir itself.'],
                include         : [type: 'string', description: 'Optional glob restricting which files are searched, e.g. `*.tsv`. Matched against the file name, or against the path relative to the search root when it contains a `/`.'],
                case_insensitive: [type: 'boolean', description: 'Match case-insensitively. Defaults to false.'],
                max_results     : [type: 'integer', description: "Maximum number of matching lines to return. Defaults to ${DEFAULT_MAX_RESULTS}, capped at ${MAX_MAX_RESULTS}.".toString()],
                max_depth       : [type: 'integer', description: "Maximum directory depth to descend. Defaults to ${DEFAULT_MAX_DEPTH}, capped at ${MAX_MAX_DEPTH}.".toString()],
            ] as Map<String,Object>,
            ['pattern'] )
        final desc = 'Search file contents line by line for a regular expression, recursively. ' +
            'Returns one entry per matching line with its absolute `file`, 1-based `line` number and the matched line `text` ' +
            "(truncated to ${MAX_LINE_CHARS} characters). " +
            "Binary files and files larger than ${MAX_GREP_FILE_BYTES.intdiv(1024 * 1024)} MB are skipped. " +
            'When more lines match than the limit allows the result sets `truncated: true` with `truncated_reason: "' + TRUNCATED_RESULTS + '"` and reports the `limit` used — narrow the pattern or the `include` glob rather than assuming you saw everything. ' +
            "The search itself is also bounded (at most ${MAX_VISITED_ENTRIES} entries visited and ${MAX_GREP_FILES} files read): when THAT stops it, `truncated_reason` is `\"${TRUNCATED_SEARCH}\"` and the tree was NOT searched to the end — narrow the search root or the depth. " + SANDBOX
        return new ToolDescriptor(GREP, desc, input, null)
    }

    private static ToolDescriptor findDescriptor() {
        final input = ToolSchema.object(
            [
                pattern    : [type: 'string', description: 'Glob matched against the file NAME, e.g. `*.fastq.gz`. A pattern containing a `/` is matched against the path relative to the search root instead, e.g. `**/results/*.json`.'],
                path       : [type: 'string', description: 'Directory to search, relative to the agent work dir. Defaults to the work dir itself.'],
                type       : [type: 'string', enum: ['file','dir','any'], description: 'Restrict results to regular files or to directories. Defaults to `any`.'],
                max_results: [type: 'integer', description: "Maximum number of paths to return. Defaults to ${DEFAULT_MAX_RESULTS}, capped at ${MAX_MAX_RESULTS}.".toString()],
                max_depth  : [type: 'integer', description: "Maximum directory depth to descend. Defaults to ${DEFAULT_MAX_DEPTH}, capped at ${MAX_MAX_DEPTH}.".toString()],
            ] as Map<String,Object>,
            ['pattern'] )
        final desc = 'Find files and directories by name, recursively, returning their absolute paths. ' +
            "This searches NAMES only — use `${GREP}` to search file contents. " +
            'When more paths match than the limit allows the result sets `truncated: true` with `truncated_reason: "' + TRUNCATED_RESULTS + '"` and reports the `limit` used. ' +
            "The search itself is also bounded (at most ${MAX_VISITED_ENTRIES} entries visited): when THAT stops it, `truncated_reason` is `\"${TRUNCATED_SEARCH}\"` and the tree was NOT searched to the end — narrow the search root or the depth. " + SANDBOX
        return new ToolDescriptor(FIND, desc, input, null)
    }

    // -------------------------------------------------------------------------
    // The `fs:` family (read / write / edit / ls / grep / find).
    //
    // These are the driver-JVM implementation of the runner-native fs: tools: an in-container
    // runner serves the same six wire names with its own builtins, so what changes between
    // runners is the implementation and the sandbox mechanism, never the name or the contract.
    // Every one of them is gated by SandboxGuard and returns {"error": ...} instead of throwing,
    // because a refused path is something the model must be able to see and correct.
    // -------------------------------------------------------------------------

    /**
     * Dispatch one {@code fs:} tool call. Resolves the LLM-supplied {@code path} against the
     * dispatch context's work dir when relative, gates it through {@link SandboxGuard}, then
     * routes to the per-tool implementation.
     *
     * <p>Returns {@code {"error": "..."}} when:
     * <ul>
     *   <li>no dispatch context is active on the current thread (called outside a sandboxed invocation)</li>
     *   <li>the context work dir is not a local {@code file:} path</li>
     *   <li>the resolved path is outside the sandbox (workDir or whitelisted readable paths)</li>
     * </ul>
     *
     * <p><b>Threading note</b>: filesystem calls carry their own ThreadLocal invocation context
     * and do not need the module request gateway. The ThreadLocal itself belongs to
     * {@link ModuleToolBridge}, whose {@code setContext}/{@code clearContext} are the public API
     * the agent task brackets an invocation with; the context is passed IN here.
     *
     * @param tool           the wire name of the tool being called, one of {@link #NAMES}
     * @param args           the parsed LLM arguments
     * @param ctx            the per-invocation dispatch context, or {@code null} when called
     *                       outside a sandboxed invocation
     * @param maxInlineBytes the cap under which {@code read} inlines a file's content rather than
     *                       returning an opaque path handle
     * @return a JSON object string with the tool result or an {@code {"error": "..."}} JSON
     */
    static String call(String tool, Map args, DispatchContext ctx, long maxInlineBytes) {
        final String unavailable = requireSandbox(ctx, tool)
        if( unavailable != null )
            return unavailable
        final Path workDir = ctx.workDir

        final PathArg arg = resolvePath(tool, args, workDir)
        if( arg.error != null )
            return arg.error
        final Path resolved = arg.resolved
        final String pathStr = arg.pathStr

        final String refused = requireAllowed(tool, resolved, workDir, ctx, pathStr)
        if( refused != null )
            return refused

        switch( tool ) {
            case READ:  return fsRead(resolved, pathStr, maxInlineBytes)
            case WRITE: return fsWrite(resolved, args)
            case EDIT:  return fsEdit(resolved, pathStr, args)
            case LS:    return fsList(resolved, pathStr, ctx)
            case GREP:  return fsGrep(resolved, pathStr, args, ctx)
            case FIND:  return fsFind(resolved, pathStr, args, ctx)
            default:
                // unreachable: `call` only routes names in `filesystemTools`
                return error("unknown filesystem tool `${tool}`")
        }
    }

    /**
     * The sandbox must exist before a path can be resolved against it: a dispatch context, a work
     * dir, and a work dir that is a local {@code file:} path.
     *
     * @return the {@code {"error": …}} result to hand back verbatim, or {@code null} when usable
     */
    private static String requireSandbox(DispatchContext ctx, String tool) {
        if( ctx == null )
            return error("`${tool}` tool unavailable: no sandbox context")

        final Path workDir = ctx.workDir
        if( workDir == null )
            return error("`${tool}` tool unavailable: no work dir in sandbox context")

        // Only local file: paths are supported; cloud/remote work dirs are not supported yet
        try { workDir.toUri() } catch( UnsupportedOperationException e ) {
            log.debug("Agent `${tool}` tool: work dir `${workDir}` has no URI (non-local scheme): ${e.message}")
            return error("`${tool}` tool unavailable: work dir scheme is not a local file path")
        }
        final scheme = workDir.toUri()?.getScheme()
        if( scheme != null && scheme != 'file' )
            return error("`${tool}` tool unavailable: work dir scheme '${scheme}' is not supported (only local file: paths)")
        return null
    }

    /**
     * The {@code path} argument, resolved — or the error result that says it could not be.
     *
     * <p>Both halves travel: the normalized {@code resolved} path is what the tool acts on, while
     * the RAW {@code pathStr} is what the six-way dispatch hands to the per-tool "not found"
     * messages and what a sandbox refusal quotes. Reconstructing the latter from the former would
     * change those messages.
     */
    @TupleConstructor
    private static class PathArg {
        final String error
        final Path resolved
        final String pathStr
    }

    private static PathArg resolvePath(String tool, Map args, Path workDir) {
        // `path` addresses the single target for read/write/edit and the search root for
        // ls/grep/find, where it defaults to the work dir itself
        final boolean pathRequired = tool==READ || tool==WRITE || tool==EDIT
        final String pathStr = args.path?.toString() ?: (pathRequired ? null : '.')
        if( !pathStr )
            return new PathArg(error("`${tool}` tool: missing required argument: path"), null, null)

        // Resolve path: relative paths are resolved against the work dir; normalize for defensive correctness
        Path resolved = Path.of(pathStr)
        if( !resolved.isAbsolute() )
            resolved = workDir.resolve(pathStr).normalize()
        else
            resolved = resolved.normalize()
        return new PathArg(null, resolved, pathStr)
    }

    /**
     * The sandbox gate. Whether the call MUTATES is derived from the tool name here and nowhere
     * else, so the read/write asymmetry has exactly one definition.
     *
     * @return the {@code {"error": …}} result to hand back verbatim, or {@code null} when allowed
     */
    private static String requireAllowed(String tool, Path resolved, Path workDir, DispatchContext ctx, String pathStr) {
        // write and edit are the only mutating tools, and the guard confines them to the work dir
        // (a module-output path whitelisted for reading is NOT writable)
        final boolean isWrite = tool==WRITE || tool==EDIT
        if( !SandboxGuard.isAllowed(resolved, workDir, ctx.readablePaths, isWrite) )
            return error("path outside sandbox: ${pathStr}")
        return null
    }

    /** A tool-result error object. Every fs: failure is a RESULT the model can read and retry from. */
    private static String error(Object message) {
        return JsonOutput.toJson([error: message.toString()])
    }

    /**
     * Read one file, applying the same inline-vs-handle policy as a module tool output
     * ({@link ToolOutputReader#readOrHandle} under the {@code agent.maxToolOutputInlineSize} cap).
     *
     * <p>The result SHAPE says which of the two happened: {@code content} is only ever real file
     * content, and a file that stayed an opaque handle comes back as {@code path} + {@code note}.
     * Putting a path string under {@code content} — as the aggregate tool this replaces did —
     * reads to the model as "the file contains this path", which is a lie it cannot detect.
     */
    private static String fsRead(Path resolved, String pathStr, long maxInlineBytes) {
        if( !Files.exists(resolved) )
            return error("file not found: ${pathStr}")
        if( Files.isDirectory(resolved) )
            return error("path is a directory, not a file: ${pathStr}")
        // a non text-like format is bulk/binary data: chainable as a handle, never inlined
        final ext = ToolOutputReader.extensionOf(resolved)
        if( !ToolOutputReader.TEXT_EXTENSIONS.contains(ext) )
            return JsonOutput.toJson([
                path: resolved.toAbsolutePath().toString(),
                note: "content not inlined: `${ext ?: '(no extension)'}` is not a text-like format; pass this path to another tool".toString() ])
        final Object readResult = ToolOutputReader.readOrHandle(resolved, maxInlineBytes)
        // an oversized inline candidate arrives as an annotated [path, note] map
        if( readResult instanceof Map )
            return JsonOutput.toJson(readResult)
        // readOrHandle returns a bare String for BOTH the inlined content and its binary-content
        // safety net, so re-run the (8 KB) sniff to tell the two apart rather than labelling a
        // path handle as content
        if( ToolOutputReader.looksBinary(resolved) )
            return JsonOutput.toJson([
                path: resolved.toAbsolutePath().toString(),
                note: 'content not inlined: the file contains binary data' ])
        return JsonOutput.toJson([content: readResult])
    }

    /**
     * Replace a file's entire contents.
     *
     * <p>A missing {@code content} is an ERROR, not an empty write: {@code content} is a required
     * property of the schema, and treating its absence as {@code ''} would TRUNCATE an existing
     * file on a malformed call the model never intended. Empty is still writable — it just has to
     * be asked for explicitly, exactly as {@code edit} requires of {@code new_string}.
     */
    private static String fsWrite(Path resolved, Map args) {
        if( !args.containsKey('content') || args.content == null )
            return error('`write` tool: missing required argument: content (pass an empty string to create an empty file)')
        final byte[] bytes = args.content.toString().getBytes(StandardCharsets.UTF_8)
        final parent = resolved.getParent()
        if( parent != null )
            Files.createDirectories(parent)
        Files.write(resolved, bytes)
        return JsonOutput.toJson([path: resolved.toAbsolutePath().toString(), bytes: bytes.length])
    }

    /**
     * Literal old-string/new-string replacement.
     *
     * <p>A non-unique {@code old_string} without {@code replace_all} is an ERROR, never a
     * silent first-match edit: the model asked to change "the" occurrence, and if there are
     * three of them it does not yet know which one it changed. Reporting the count back lets it
     * either extend the match with surrounding context or opt into replacing all of them.
     */
    private static String fsEdit(Path resolved, String pathStr, Map args) {
        if( !Files.exists(resolved) )
            return error("file not found: ${pathStr}")
        if( Files.isDirectory(resolved) )
            return error("path is a directory, not a file: ${pathStr}")
        final String oldStr = args.old_string?.toString()
        if( !oldStr )
            return error('`edit` tool: missing required argument: old_string (it must be a non-empty literal string)')
        // new_string may legitimately be empty (a deletion), so it is checked for presence only
        if( !args.containsKey('new_string') || args.new_string == null )
            return error('`edit` tool: missing required argument: new_string (pass an empty string to delete the matched text)')
        final String newStr = args.new_string.toString()
        if( oldStr == newStr )
            return error('`edit` tool: old_string and new_string are identical - nothing to change')
        final boolean replaceAll = asBoolean(args.replace_all)

        final String text = new String(Files.readAllBytes(resolved), StandardCharsets.UTF_8)
        final int occurrences = countOccurrences(text, oldStr)
        if( occurrences == 0 )
            return error("`edit` tool: no match for old_string in ${pathStr} - it must match the file contents exactly, including whitespace and indentation")
        if( occurrences > 1 && !replaceAll )
            return error("`edit` tool: old_string occurs ${occurrences} times in ${pathStr} - extend it with surrounding context to make it unique, or pass replace_all: true to change all ${occurrences} occurrences")

        // literal replacement: replace() (unlike replaceAll()) treats both arguments as literals
        final String updated = replaceAll
                ? text.replace(oldStr, newStr)
                : replaceFirstLiteral(text, oldStr, newStr)
        final byte[] bytes = updated.getBytes(StandardCharsets.UTF_8)
        Files.write(resolved, bytes)
        return JsonOutput.toJson([
            path        : resolved.toAbsolutePath().toString(),
            replacements: replaceAll ? occurrences : 1,
            bytes       : bytes.length ])
    }

    /** Count the NON-OVERLAPPING occurrences of a literal string, i.e. as many as replace() makes. */
    private static int countOccurrences(String text, String literal) {
        int count = 0
        int from = 0
        while( true ) {
            final int at = text.indexOf(literal, from)
            if( at < 0 )
                return count
            count++
            from = at + literal.length()
        }
    }

    private static String replaceFirstLiteral(String text, String oldStr, String newStr) {
        final int at = text.indexOf(oldStr)
        return text.substring(0, at) + newStr + text.substring(at + oldStr.length())
    }

    /**
     * List one directory's immediate entries.
     *
     * <p>Every entry is re-checked against the sandbox for the same reason {@link #fsGrep} and
     * {@link #fsFind} re-check theirs: {@code Files.isDirectory}/{@code Files.size} FOLLOW
     * symlinks, so a link pointing outside the sandbox would report the kind and the exact byte
     * size of a file the model is not allowed to open. Such an entry is reported as
     * {@code type: 'link'} with no size — its name is sandbox-internal and hiding it would only
     * make the refusal that {@code read} then returns inexplicable.
     */
    private static String fsList(Path resolved, String pathStr, DispatchContext ctx) {
        if( !Files.exists(resolved) )
            return error("directory not found: ${pathStr}")
        if( !Files.isDirectory(resolved) )
            return error("path is not a directory: ${pathStr}")
        final entries = new ArrayList<Map>()
        final Stream<Path> stream = Files.list(resolved)
        try {
            for( final Iterator<Path> it=stream.sorted().iterator(); it.hasNext(); ) {
                final Path p = it.next()
                final entry = new LinkedHashMap<String,Object>()
                entry.put('name', p.getFileName().toString())
                if( !SandboxGuard.isAllowed(p, ctx.workDir, ctx.readablePaths, false) ) {
                    entry.put('type', 'link')
                    entries.add(entry)
                    continue
                }
                final boolean dir = Files.isDirectory(p)
                entry.put('type', dir ? 'dir' : 'file')
                // the size lets the model decide whether reading the file is worthwhile
                if( !dir )
                    entry.put('size', safeSize(p))
                entries.add(entry)
            }
        }
        finally {
            stream.close()
        }
        return JsonOutput.toJson([path: resolved.toAbsolutePath().toString(), entries: entries])
    }

    /**
     * Content search, pure JVM. Deliberately NOT a shell-out to {@code grep}: the driver may be
     * macOS, where BSD grep takes different flags from the GNU one an LLM will have been trained
     * to write, and a flag mismatch surfaces as an unexplainable empty result.
     *
     * <p>Every visited file is re-checked against the sandbox even though the search ROOT was
     * already checked: {@code Files.walk} does not follow directory symlinks, but it does report
     * a symlink to a regular file, and reading through it would leak a file outside the sandbox.
     * Such an entry is skipped silently — it is not the model's error to correct.
     */
    private static String fsGrep(Path root, String pathStr, Map args, DispatchContext ctx) {
        final String patternStr = args.pattern?.toString()
        if( !patternStr )
            return error('`grep` tool: missing required argument: pattern')
        Pattern regex
        try {
            regex = Pattern.compile(patternStr, asBoolean(args.case_insensitive) ? Pattern.CASE_INSENSITIVE : 0)
        }
        catch( PatternSyntaxException e ) {
            return error("`grep` tool: invalid regular expression `${patternStr}` - ${e.description}")
        }
        if( !Files.exists(root) )
            return error("path not found: ${pathStr}")
        PathMatcher include
        try {
            include = globMatcher(args.include?.toString())
        }
        catch( Exception e ) {
            return error("`grep` tool: invalid include glob `${args.include}` - ${e.message}")
        }
        final int limit = boundedArg(args.max_results, DEFAULT_MAX_RESULTS, MAX_MAX_RESULTS)
        final int depth = boundedArg(args.max_depth, DEFAULT_MAX_DEPTH, MAX_MAX_DEPTH)

        final matches = new ArrayList<Map>()
        int filesScanned = 0
        int visited = 0
        String truncatedReason = null
        final Stream<Path> stream = Files.walk(root, depth)
        try {
            for( final Iterator<Path> it=stream.iterator(); it.hasNext() && truncatedReason == null; ) {
                final Path p = it.next()
                // bound the WALK, not just the answer: a pattern matching nothing would otherwise
                // stat + realpath + sniff + read every entry of an arbitrarily large tree on the
                // agent dispatch thread. Reported distinctly so an incomplete search is not read
                // by the model as an exhaustive one that found nothing.
                if( ++visited > MAX_VISITED_ENTRIES ) {
                    truncatedReason = TRUNCATED_SEARCH
                    break
                }
                if( !Files.isRegularFile(p) )
                    continue
                if( !SandboxGuard.isAllowed(p, ctx.workDir, ctx.readablePaths, false) )
                    continue
                if( include != null && !matchesGlob(include, root, p) )
                    continue
                // a match inside a multi-gigabyte data file is noise, and reading it into the
                // driver heap to find out would be worse than useless
                if( safeSize(p) > MAX_GREP_FILE_BYTES || ToolOutputReader.looksBinary(p) )
                    continue
                if( filesScanned >= MAX_GREP_FILES ) {
                    truncatedReason = TRUNCATED_SEARCH
                    break
                }
                filesScanned++
                // decode with replacement rather than readAllLines, which throws on a malformed
                // byte sequence and would abort the whole search over one bad file
                final String text = new String(Files.readAllBytes(p), StandardCharsets.UTF_8)
                final String[] lines = text.split('\n', -1)
                for( int i=0; i<lines.length; i++ ) {
                    if( !regex.matcher(lines[i]).find() )
                        continue
                    if( matches.size() >= limit ) {
                        truncatedReason = TRUNCATED_RESULTS
                        break
                    }
                    matches.add([
                        file: p.toAbsolutePath().toString(),
                        line: i + 1,
                        text: truncate(lines[i], MAX_LINE_CHARS) ] as Map)
                }
            }
        }
        finally {
            stream.close()
        }
        final result = new LinkedHashMap<String,Object>()
        result.put('matches', matches)
        result.put('count', matches.size())
        result.put('files_scanned', filesScanned)
        result.put('truncated', truncatedReason != null)
        if( truncatedReason != null )
            result.put('truncated_reason', truncatedReason)
        result.put('limit', limit)
        result.put('max_depth', depth)
        return JsonOutput.toJson(result)
    }

    /**
     * Name search, pure JVM {@link PathMatcher} — same reasoning as {@link #fsGrep}: no shelling
     * out to {@code find}, whose predicate syntax differs between BSD and GNU.
     */
    private static String fsFind(Path root, String pathStr, Map args, DispatchContext ctx) {
        final String patternStr = args.pattern?.toString()
        if( !patternStr )
            return error('`find` tool: missing required argument: pattern')
        PathMatcher matcher
        try {
            matcher = globMatcher(patternStr)
        }
        catch( Exception e ) {
            return error("`find` tool: invalid glob pattern `${patternStr}` - ${e.message}")
        }
        if( !Files.exists(root) )
            return error("path not found: ${pathStr}")
        if( !Files.isDirectory(root) )
            return error("path is not a directory: ${pathStr}")
        final String type = args.type?.toString() ?: 'any'
        if( !(type in ['file','dir','any']) )
            return error("`find` tool: unknown type '${type}'; supported: file, dir, any")
        final int limit = boundedArg(args.max_results, DEFAULT_MAX_RESULTS, MAX_MAX_RESULTS)
        final int depth = boundedArg(args.max_depth, DEFAULT_MAX_DEPTH, MAX_MAX_DEPTH)

        final paths = new ArrayList<String>()
        int visited = 0
        String truncatedReason = null
        final Stream<Path> stream = Files.walk(root, depth)
        try {
            for( final Iterator<Path> it=stream.iterator(); it.hasNext(); ) {
                final Path p = it.next()
                // bound the WALK as well as the answer -- see the same guard in fsGrep
                if( ++visited > MAX_VISITED_ENTRIES ) {
                    truncatedReason = TRUNCATED_SEARCH
                    break
                }
                // the root itself is the thing being searched, never a result
                if( p == root )
                    continue
                final boolean dir = Files.isDirectory(p)
                if( type=='file' && dir )
                    continue
                if( type=='dir' && !dir )
                    continue
                // as in grep: a symlink to a target outside the sandbox must not be reported
                if( !SandboxGuard.isAllowed(p, ctx.workDir, ctx.readablePaths, false) )
                    continue
                if( !matchesGlob(matcher, root, p) )
                    continue
                if( paths.size() >= limit ) {
                    truncatedReason = TRUNCATED_RESULTS
                    break
                }
                paths.add(p.toAbsolutePath().toString())
            }
        }
        finally {
            stream.close()
        }
        Collections.sort(paths)
        final result = new LinkedHashMap<String,Object>()
        result.put('paths', paths)
        result.put('count', paths.size())
        result.put('truncated', truncatedReason != null)
        if( truncatedReason != null )
            result.put('truncated_reason', truncatedReason)
        result.put('limit', limit)
        result.put('max_depth', depth)
        return JsonOutput.toJson(result)
    }

    /** Compile a glob into a matcher; {@code null} for a null/blank pattern (i.e. no filtering). */
    private static PathMatcher globMatcher(String glob) {
        if( !glob )
            return null
        return FileSystems.getDefault().getPathMatcher("glob:${glob}")
    }

    /**
     * A path-shaped glob (one containing a {@code /}) is matched against the path RELATIVE to the
     * search root, so a pattern is portable between search roots; a bare glob is matched against
     * the file name alone, which is what {@code *.json} is expected to mean.
     */
    private static boolean matchesGlob(PathMatcher matcher, Path root, Path candidate) {
        final Path name = candidate.getFileName()
        if( name != null && matcher.matches(name) )
            return true
        try {
            return matcher.matches(root.relativize(candidate))
        }
        catch( IllegalArgumentException e ) {
            return false
        }
    }

    /** An LLM-supplied numeric bound, defaulted when absent/unparseable and clamped to {@code [1,max]}. */
    private static int boundedArg(Object value, int defaultValue, int max) {
        int result = defaultValue
        if( value instanceof Number )
            result = ((Number) value).intValue()
        else if( value != null ) {
            try { result = Integer.parseInt(value.toString().trim()) }
            catch( NumberFormatException e ) { result = defaultValue }
        }
        return Math.max(1, Math.min(result, max))
    }

    /** JSON booleans arrive as Boolean, but a model that stringifies them should still be understood. */
    private static boolean asBoolean(Object value) {
        if( value instanceof Boolean )
            return ((Boolean) value).booleanValue()
        return value != null && Boolean.parseBoolean(value.toString().trim())
    }

    private static String truncate(String line, int max) {
        return line.length() <= max ? line : line.substring(0, max) + '...'
    }

    /** File size, or {@code -1} when it cannot be read (a broken symlink, a vanished file). */
    private static long safeSize(Path path) {
        try { return Files.size(path) }
        catch( IOException e ) { return -1 }
    }
}

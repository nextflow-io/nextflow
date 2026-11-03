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

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.npr.api.schema.v1.ModuleMetadata
import io.seqera.npr.client.RegistryClient
import nextflow.Global
import nextflow.Session
import nextflow.config.RegistryConfig
import nextflow.module.ModuleInfo
import nextflow.module.ModuleReference
import nextflow.module.ModuleSpec
import nextflow.module.ModuleSpecFactory
import nextflow.module.RegistryClientFactory
import nextflow.script.BaseScript
import nextflow.script.ProcessDef
import nextflow.script.ScriptMeta

/**
 * Turns the brokered half of an agent's resolved {@code tools} selection into a
 * {@link ModuleToolBridge}: for every {@code nf:module_run} process in scope it finds the sibling
 * {@code meta.yml} spec and, when the module came from the registry, its public
 * {@link ModuleMetadata} — the three facts a tool is wired from, collected into one
 * {@link WiredModuleTool} each.
 *
 * <p>Owning this keeps the lookup chain (module dir -> registry install marker -> registry client
 * -> metadata) out of the agent definition, which only needs the bridge that comes out of it.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ModuleToolResolver {

    /** The agent the tools are wired for; named by the debug/warn messages below. */
    private final String agentName

    /** The script the agent is defined in, i.e. the scope the processes are looked up in. */
    private final BaseScript owner

    ModuleToolResolver(String agentName, BaseScript owner) {
        this.agentName = agentName
        this.owner = owner
    }

    /**
     * Wire the <b>brokered</b> half of the resolved selection — the {@code nf:module_run}
     * processes the driver executes as real Nextflow tasks — into a {@link ModuleToolBridge},
     * and, for an in-JVM runner only, the {@code fs:} leaves it must serve itself.
     *
     * <p>Which {@code fs:} leaves the bridge serves is decided by WHO EXECUTES the agent, not by
     * whether any were declared:
     * <ul>
     *   <li><b>containerized runner</b> — none. The {@code fs:} tools are the runner's own
     *       builtins, rooted at the session cwd inside the container, enabled by the names on
     *       {@code AgentRunnerRequest.nativeToolNames}. Handing them to the bridge would relocate
     *       a container-side tool into the driver JVM, where the sandbox is rooted at a
     *       driver-side work-dir path that on a cloud executor is not even a local file.</li>
     *   <li><b>in-JVM runner</b> — exactly the selected leaves, never the whole family: the
     *       resolver already expanded {@code fs:*} where that is what was written, so collapsing
     *       a partial selection here would hand a read-only agent {@code write} and {@code edit}.
     *       They are served through {@link ModuleToolBridge#call} but stay OUT of
     *       {@code toolSpecs} (§5); the runner advertises them from the same native names.</li>
     * </ul>
     *
     * @param selection the resolved {@code tools} selection, already partitioned into its brokered
     *      and runner-native halves; {@code null} or empty when the agent declared no tools
     * @param containerized whether the selected runner executes the agent in its own container
     * @return a {@link ModuleToolBridge}, or {@code null} when nothing was declared
     */
    ModuleToolBridge createToolBridge(ToolRefResolver.Selection selection, boolean containerized) {
        if( selection == null || selection.isEmpty() )
            return null
        final meta = ScriptMeta.get(owner)
        final List<WiredModuleTool> wired = wireModuleRunTools(meta, selection.brokeredNames)
        final List<String> inJvmFsTools = containerized
                ? Collections.<String>emptyList()
                : selection.runnerNative
                        .findAll { it.ref.startsWith(ToolRefResolver.FS_FAMILY + ':') }
                        .collect { it.name }
        return new ModuleToolBridge(wired, inJvmFsTools)
    }

    /**
     * Wire each selected in-scope process as its own tool with an enforced per-module schema.
     * The descriptor is sourced from the public registry
     * {@link io.seqera.npr.api.schema.v1.ModuleMetadata} when the module is a registry install
     * (richer), else from the sibling {@code meta.yml} {@link nextflow.module.ModuleSpec}. One
     * registry client is built for the whole selection, and only when at least one module needs it,
     * so a purely-local selection never constructs one.
     *
     * @param procNames the process names selected through {@code nf:module_run}, already
     *                  de-duplicated and verified to exist by the resolver
     * @return one {@link WiredModuleTool} per resolvable process, in {@code procNames} order
     */
    private List<WiredModuleTool> wireModuleRunTools(ScriptMeta meta, Collection<String> procNames) {
        if( meta == null || !procNames )
            return new ArrayList<WiredModuleTool>()
        // resolve the module refs FIRST: that is what decides whether a registry client is needed
        // at all, so a purely-local selection never builds one, and a mixed selection builds one
        final Map<String,ProcessDef> procs = resolvableProcesses(meta, procNames)
        final Map<String,ModuleReference> refs = procs.collectEntries { name, proc ->
            [(name): recoverModuleRef(resolveIncludedModuleDir(proc))]
        } as Map<String,ModuleReference>
        final RegistryClient client = refs.values().any { it != null }
                ? newRegistryClient(Global.session as Session)
                : null
        return procs.collect { name, proc -> wireOne(name, proc, refs.get(name), client) }
    }

    /** The selected names that name a process in scope, in selection order. */
    private static Map<String,ProcessDef> resolvableProcesses(ScriptMeta meta, Collection<String> procNames) {
        final Map<String,ProcessDef> result = new LinkedHashMap<String,ProcessDef>()
        for( final name : procNames ) {
            final proc = meta.getProcess(name)
            if( proc != null )
                result.put(name, proc)
        }
        return result
    }

    /**
     * One process, wired as one tool: its sibling {@code meta.yml} spec plus — when the module came
     * from the registry and the metadata could be fetched — the richer registry descriptor source.
     */
    private WiredModuleTool wireOne(String procName, ProcessDef proc, ModuleReference moduleRef, RegistryClient client) {
        // sibling meta.yml for spec-driven marshalling; null for locally-defined processes
        final spec = loadSiblingSpec(proc)
        final metadata = moduleRef != null ? registryMetadata(client, moduleRef) : null
        // the nf-core meta.id convention is a property of the SOURCE, so it only applies when the
        // registry metadata is what the descriptor is built from
        final nfCore = metadata != null && moduleRef.scope == 'nf-core'
        return new WiredModuleTool(procName, proc, spec, metadata, nfCore)
    }

    /**
     * The registry metadata for an installed module, or {@code null} when it cannot be had — the
     * descriptor then falls back to the sibling {@code meta.yml} spec, which is always present for
     * a registry install.
     */
    private ModuleMetadata registryMetadata(RegistryClient client, ModuleReference moduleRef) {
        try {
            return fetchModuleMetadata(client, moduleRef, null)
        }
        catch( Exception e ) {
            log.debug("Agent `${agentName}` nf:module_run: could not fetch registry metadata for `${moduleRef.fullName}` (${e.message}) — falling back to meta.yml spec")
            return null
        }
    }

    /** Build a {@link io.seqera.npr.client.RegistryClient} from the session's {@code registry} config scope. */
    private static RegistryClient newRegistryClient(Session session) {
        final registryConfig = new RegistryConfig((session?.config?.registry as Map) ?: Collections.emptyMap())
        return RegistryClientFactory.forConfig(registryConfig)
    }

    /**
     * Look for a sibling {@code meta.yml} / {@code meta.yaml} in the module dir of the script
     * that defines the process and, when present, load it into a {@link nextflow.module.ModuleSpec}
     * to drive spec-driven tool schema and tuple/path/map marshalling (Phase 3.2). Returns
     * {@code null} when no sibling spec is found.
     */
    private static ModuleSpec loadSiblingSpec(ProcessDef proc) {
        final dir = resolveIncludedModuleDir(proc)
        if( dir == null )
            return null
        for( final candidate : ['meta.yml', 'meta.yaml'] ) {
            final specPath = dir.resolve(candidate)
            if( specPath.toFile().exists() )
                return ModuleSpecFactory.fromYaml(specPath)
        }
        return null
    }

    /**
     * The module dir of an included process (brought in via {@code include { X } from '...'}),
     * derived from the {@link ScriptMeta} of the script that defines it. Returns {@code null}
     * when the owner, its ScriptMeta, or its module dir cannot be resolved.
     */
    private static Path resolveIncludedModuleDir(ProcessDef proc) {
        final BaseScript owner = proc.getOwner()
        if( owner == null )
            return null
        final ScriptMeta scriptMeta = ScriptMeta.get(owner)
        if( scriptMeta == null )
            return null
        return scriptMeta.getModuleDir()
    }

    /**
     * Recover a {@link nextflow.module.ModuleReference} from an included module's install dir
     * when it is a registry install. A registry install is identified by the presence of the
     * {@link nextflow.module.ModuleInfo#MODULE_INFO_FILE} marker ({@code .module-info}) inside the
     * module dir, AND the directory layout {@code <base>/modules/<scope>/<name>} (depth ≥ 2 relative
     * to the markers parent's parent so both {@code scope} and {@code name} components exist).
     *
     * <p>Returns {@code null} for any case where recovery is infeasible or unsafe:
     * <ul>
     *   <li>the {@code moduleDir} is null;</li>
     *   <li>the {@code .module-info} marker is not present;</li>
     *   <li>the parent chain does not have a grandparent named {@code modules};</li>
     *   <li>any exception during path inspection.</li>
     * </ul>
     * Local-file {@code include} statements (no marker file) fall through to {@code null}.
     */
    private static ModuleReference recoverModuleRef(Path moduleDir) {
        if( moduleDir == null )
            return null
        try {
            // must have the .module-info marker to be a registry install
            if( !Files.exists(moduleDir.resolve(ModuleInfo.MODULE_INFO_FILE)) )
                return null
            // layout: <anything>/modules/<scope>/<name>
            // dir.fileName = <name>, dir.parent.fileName = <scope>, dir.parent.parent.fileName = modules
            final nameComp = moduleDir.fileName
            final scopeDir = moduleDir.parent
            if( nameComp == null || scopeDir == null )
                return null
            final scopeComp = scopeDir.fileName
            final modulesDir = scopeDir.parent
            if( scopeComp == null || modulesDir == null )
                return null
            if( modulesDir.fileName?.toString() != 'modules' )
                return null
            return new ModuleReference(scopeComp.toString(), nameComp.toString())
        }
        catch( Exception e ) {
            log.debug("recoverModuleRef: could not recover module reference from `${moduleDir}`: ${e.message}")
            return null
        }
    }

    /**
     * Fetch the public {@link io.seqera.npr.api.schema.v1.ModuleMetadata} for a resolved registry
     * module via the SAME {@link io.seqera.npr.client.RegistryClient} already built to resolve the
     * module ({@code GET /api/v1/modules/{name}} is anonymous/public). Any failure (network hiccup,
     * private module without metadata, missing release) is logged and yields {@code null} so the
     * caller degrades gracefully to the sibling {@code meta.yml} {@link nextflow.module.ModuleSpec}.
     */
    private ModuleMetadata fetchModuleMetadata(RegistryClient client, ModuleReference moduleRef, String version) {
        try {
            if( version )
                return client.getModuleRelease(moduleRef.fullName, version)?.metadata
            return client.getModule(moduleRef.fullName)?.latest?.metadata
        }
        catch( Exception e ) {
            log.warn("Agent `${agentName}` tool `${moduleRef.fullName}`: could not fetch registry metadata (${e.message}) - falling back to the local module spec for the tool descriptor")
            return null
        }
    }

}

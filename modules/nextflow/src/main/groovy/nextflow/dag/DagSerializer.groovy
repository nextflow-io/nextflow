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

package nextflow.dag

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

/**
 * Serialize a {@link DAG} into a plain {@code Map} that mirrors the DAG model
 * one-to-one — every vertex and every edge, with no filtering, reduction or
 * reinterpretation.
 *
 * Unlike the file-oriented renderers in this package ({@link DotRenderer},
 * {@link GexfRenderer}, {@link MermaidRenderer}, ...) which emit a layout-ready
 * document, this serializer produces a faithful machine-readable dump intended
 * for transport to external consumers (e.g. Seqera Platform). Any interpretation
 * — collapsing operator/channel vertices to a process-only graph, building a
 * subworkflow hierarchy, choosing a layout — is left to the consumer, so it can
 * evolve without a Nextflow release.
 *
 * Output shape:
 * <pre>
 * {
 *   vertices: [ { id, type, label, scope } ... ],
 *   edges:    [ { source, target, label } ... ]
 * }
 * </pre>
 *
 * Notes:
 * <ul>
 *   <li>{@code id} is the payload-scoped per-vertex id ({@link DAG.Vertex#getId()})
 *       used to join edges to vertices within this payload. It comes from a static
 *       JVM-global counter, so it is <em>not</em> deterministic across runs and will
 *       not match between a run and its {@code -resume} — treat it as identity only
 *       within one serialized payload, never as cross-run identity. It is preferred
 *       over the positional {@code name} ({@code v0}, {@code v1}, ...) because
 *       {@link DAG#normalize()} inserts boundary vertices and shifts positions.</li>
 *   <li>for {@code PROCESS} vertices {@code label} carries the fully-qualified
 *       process name with the scope prefix included (e.g. {@code RNASEQ:ALIGN:STAR}),
 *       while {@code scope} carries the enclosing (sub)workflow name — so the FQ
 *       name remains reconstructable without guessing.</li>
 *   <li>every vertex type ({@code PROCESS}, {@code OPERATOR}, {@code ORIGIN},
 *       {@code NODE}) and every edge is emitted as-is; operator/channel vertices
 *       between processes are preserved, never dropped.</li>
 * </ul>
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
@Slf4j
@CompileStatic
class DagSerializer {

    /**
     * Convert the given DAG into a faithful {@code Map} representation.
     *
     * The DAG is {@link DAG#normalize() normalized} first so that dangling edges
     * get their boundary ({@code ORIGIN}/{@code NODE}) vertices and channel-derived
     * edge labels are resolved — matching the behaviour of the file renderers.
     *
     * <p><b>Side effect:</b> this mutates {@code dag} in place via
     * {@link DAG#normalize()}. {@code normalize()} is idempotent, so calling this
     * alongside other normalizers (e.g. {@code GraphObserver} at completion) is safe.
     *
     * @param dag The DAG to serialize, or {@code null}
     * @return A map with {@code vertices} and {@code edges} lists, or {@code null}
     *         when the given DAG is {@code null}
     */
    static Map<String,Object> toMap(DAG dag) {
        if( dag == null )
            return null

        // resolve dangling edge endpoints and channel-derived edge labels (mutates dag)
        dag.normalize()

        final vertices = new ArrayList<Map<String,Object>>(dag.vertices.size())
        for( DAG.Vertex v : dag.vertices ) {
            final entry = new LinkedHashMap<String,Object>(4)
            entry.id = String.valueOf(v.id)
            entry.type = v.type?.toString()
            entry.label = v.label
            // createVertex defaults `workflow` to "" — coerce the empty-string root
            // scope back to null so consumers see a clean "no enclosing subworkflow"
            entry.scope = v.workflow ?: null
            vertices.add(entry)
        }

        final edges = new ArrayList<Map<String,Object>>(dag.edges.size())
        for( DAG.Edge e : dag.edges ) {
            final entry = new LinkedHashMap<String,Object>(3)
            entry.source = e.from != null ? String.valueOf(e.from.id) : null
            entry.target = e.to != null ? String.valueOf(e.to.id) : null
            entry.label = e.label
            edges.add(entry)
        }

        final result = new LinkedHashMap<String,Object>(2)
        result.vertices = vertices
        result.edges = edges
        return result
    }

}

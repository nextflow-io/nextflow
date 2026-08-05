/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.lineage.test

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import nextflow.lineage.LinNormalizer
import nextflow.lineage.LinStore

/**
 * Spock integration for lineage validation using snapshot testing.
 * 
 * This class allows pipeline tests to compare workflow runs against saved
 * baseline snapshots. On first run (or when UPDATE_SNAPSHOTS=true), the
 * normalized lineage is saved. On subsequent runs, it's compared against
 * the saved snapshot.
 *
 * Usage in Spock tests:
 * <pre>
 * class MyPipelineTest extends Specification {
 *     @Shared LineageSnapshotter snapshotter
 *     
 *     def setupSpec() {
 *         snapshotter = new LineageSnapshotter(store)
 *             .withSnapshotDir('src/test/resources/snapshots')
 *     }
 *     
 *     def "pipeline produces expected outputs"() {
 *         when:
 *         def lid = runPipeline('main.nf', params)
 *         
 *         then:
 *         snapshotter.assertMatchesSnapshot(lid, 'baseline-v1')
 *     }
 * }
 * </pre>
 *
 * Set UPDATE_SNAPSHOTS=true environment variable to update snapshots.
 *
 * @author Edmund Miller <edmund.a.miller@gmail.com>
 */
@CompileStatic
class LineageSnapshotter {

    private final LinStore store
    private final LinNormalizer normalizer
    private Path snapshotDir
    private boolean updateSnapshots

    /**
     * Create a new LineageSnapshotter
     *
     * @param store The lineage store to read workflow data from
     */
    LineageSnapshotter(LinStore store) {
        this.store = store
        this.normalizer = new LinNormalizer()
        this.updateSnapshots = System.getenv('UPDATE_SNAPSHOTS')?.toLowerCase() in ['true', '1', 'yes']
    }

    /**
     * Set the directory where snapshots are stored
     */
    LineageSnapshotter withSnapshotDir(String dir) {
        this.snapshotDir = Paths.get(dir)
        return this
    }

    /**
     * Set the directory where snapshots are stored
     */
    LineageSnapshotter withSnapshotDir(Path dir) {
        this.snapshotDir = dir
        return this
    }

    /**
     * Set the base path for relativizing file paths in outputs
     */
    LineageSnapshotter withOutputBase(String base) {
        this.normalizer.withOutputBase(base)
        return this
    }

    /**
     * Set the base path for relativizing file paths in outputs
     */
    LineageSnapshotter withOutputBase(Path base) {
        this.normalizer.withOutputBase(base)
        return this
    }

    /**
     * Add additional fields to ignore during comparison
     */
    LineageSnapshotter withIgnoreFields(Collection<String> fields) {
        this.normalizer.withIgnoreFields(fields)
        return this
    }

    /**
     * Force update mode (save snapshots instead of comparing)
     */
    LineageSnapshotter withUpdateSnapshots(boolean update) {
        this.updateSnapshots = update
        return this
    }

    /**
     * Assert that a workflow run matches a saved snapshot.
     * 
     * On first run (or when UPDATE_SNAPSHOTS=true), saves the normalized
     * lineage as the snapshot. On subsequent runs, compares against the
     * saved snapshot and throws AssertionError if they differ.
     *
     * @param lid The LID of the workflow run to validate
     * @param snapshotId Identifier for the snapshot (used as filename)
     * @throws AssertionError if the workflow run doesn't match the snapshot
     */
    void assertMatchesSnapshot(String lid, String snapshotId) {
        if (!snapshotDir) {
            throw new IllegalStateException("Snapshot directory not configured. Call withSnapshotDir() first.")
        }

        // Normalize the workflow tree
        def normalized = normalizer.normalizeWorkflowTree(store, lid)
        def actualJson = JsonOutput.prettyPrint(JsonOutput.toJson(normalized))

        // Get snapshot file path
        def snapshotFile = snapshotDir.resolve("${snapshotId}.json")

        if (updateSnapshots || !Files.exists(snapshotFile)) {
            // Save snapshot
            Files.createDirectories(snapshotFile.parent)
            snapshotFile.text = actualJson
            
            if (updateSnapshots) {
                println "Updated snapshot: ${snapshotFile}"
            } else {
                println "Created snapshot: ${snapshotFile}"
            }
            return
        }

        // Compare against saved snapshot
        def expectedJson = snapshotFile.text
        def expected = new JsonSlurper().parseText(expectedJson) as Map
        def actual = new JsonSlurper().parseText(actualJson) as Map

        def differences = LinNormalizer.compare(expected, actual)
        
        if (!differences.isEmpty()) {
            // Save actual output for debugging
            def actualFile = snapshotDir.resolve("${snapshotId}.actual.json")
            actualFile.text = actualJson
            
            def diffReport = buildDiffReport(differences)
            throw new AssertionError(
                "Workflow run does not match snapshot '${snapshotId}'.\n" +
                "Differences:\n${diffReport}\n" +
                "Expected: ${snapshotFile}\n" +
                "Actual:   ${actualFile}\n" +
                "Set UPDATE_SNAPSHOTS=true to update the snapshot."
            )
        }
    }

    /**
     * Assert that two workflow runs are semantically equivalent.
     * 
     * This compares two runs directly without using saved snapshots.
     *
     * @param lid1 LID of the first workflow run
     * @param lid2 LID of the second workflow run
     * @throws AssertionError if the workflow runs differ
     */
    void assertEquivalent(String lid1, String lid2) {
        def tree1 = normalizer.normalizeWorkflowTree(store, lid1)
        def tree2 = normalizer.normalizeWorkflowTree(store, lid2)

        def differences = LinNormalizer.compare(tree1, tree2)
        
        if (!differences.isEmpty()) {
            def diffReport = buildDiffReport(differences)
            throw new AssertionError(
                "Workflow runs are not semantically equivalent.\n" +
                "LID 1: ${lid1}\n" +
                "LID 2: ${lid2}\n" +
                "Differences:\n${diffReport}"
            )
        }
    }

    /**
     * Build a human-readable diff report from differences map
     */
    private String buildDiffReport(Map<String, Object> differences) {
        def sb = new StringBuilder()
        differences.each { String path, Object diff ->
            if (diff instanceof Map) {
                def diffMap = diff as Map
                sb.append("  ${path}:\n")
                sb.append("    expected: ${diffMap['expected']}\n")
                sb.append("    actual:   ${diffMap['actual']}\n")
                if (diffMap['reason']) {
                    sb.append("    reason:   ${diffMap['reason']}\n")
                }
            } else {
                sb.append("  ${path}: ${diff}\n")
            }
        }
        return sb.toString()
    }

    /**
     * Get the normalized representation of a workflow run.
     * Useful for debugging or manual inspection.
     *
     * @param lid The LID of the workflow run
     * @return Normalized Map representation
     */
    Map<String, Object> getNormalized(String lid) {
        return normalizer.normalizeWorkflowTree(store, lid)
    }

    /**
     * Get the normalized JSON of a workflow run.
     * Useful for debugging or manual inspection.
     *
     * @param lid The LID of the workflow run
     * @return Pretty-printed JSON string
     */
    String getNormalizedJson(String lid) {
        def normalized = normalizer.normalizeWorkflowTree(store, lid)
        return JsonOutput.prettyPrint(JsonOutput.toJson(normalized))
    }
}

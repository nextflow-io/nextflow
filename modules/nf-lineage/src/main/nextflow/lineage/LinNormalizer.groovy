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

package nextflow.lineage

import java.nio.file.Path
import java.nio.file.Paths

import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import nextflow.lineage.model.v1beta1.FileOutput
import nextflow.lineage.model.v1beta1.TaskRun
import nextflow.lineage.model.v1beta1.WorkflowOutput
import nextflow.lineage.model.v1beta1.WorkflowRun
import nextflow.lineage.serde.LinEncoder
import nextflow.lineage.serde.LinSerializable

/**
 * Normalizes lineage data for comparison by stripping ephemeral fields
 * (timestamps, session IDs, absolute paths) that vary between runs but
 * don't affect semantic equivalence.
 *
 * @author Edmund Miller <edmund.a.miller@gmail.com>
 */
@CompileStatic
class LinNormalizer {

    /**
     * Fields to remove during normalization as they vary between runs
     */
    static final Set<String> EPHEMERAL_FIELDS = [
        'sessionId',
        'createdAt', 
        'modifiedAt',
        'name',         // workflow run name is auto-generated
        'source',       // LID reference - varies per run
        'workflowRun',  // LID reference - varies per run
        'taskRun',      // LID reference - varies per run
        'path'          // absolute path - varies between runs, checksum is what matters
    ] as Set<String>

    /**
     * Base path for relativizing absolute paths in FileOutput
     */
    private Path outputBase

    /**
     * Additional fields to ignore during comparison
     */
    private Set<String> additionalIgnoreFields = [] as Set<String>

    LinNormalizer() {
    }

    /**
     * Set the base path for relativizing file paths
     */
    LinNormalizer withOutputBase(Path base) {
        this.outputBase = base
        return this
    }

    /**
     * Set the base path for relativizing file paths
     */
    LinNormalizer withOutputBase(String base) {
        this.outputBase = base ? Paths.get(base) : null
        return this
    }

    /**
     * Add additional fields to ignore during comparison
     */
    LinNormalizer withIgnoreFields(Collection<String> fields) {
        this.additionalIgnoreFields.addAll(fields)
        return this
    }

    /**
     * Normalize a lineage record by removing ephemeral fields
     * and relativizing paths
     *
     * @param record The lineage record to normalize
     * @return A normalized Map representation
     */
    Map<String, Object> normalize(LinSerializable record) {
        // Encode to JSON and parse back as a mutable Map
        def encoder = new LinEncoder()
        def json = encoder.encode(record)
        def slurper = new JsonSlurper()
        def map = slurper.parseText(json) as Map<String, Object>
        
        return normalizeMap(map, record)
    }

    /**
     * Normalize a Map representation of a lineage record
     */
    protected Map<String, Object> normalizeMap(Map<String, Object> map, LinSerializable original) {
        def result = [:] as Map<String, Object>
        def fieldsToIgnore = EPHEMERAL_FIELDS + additionalIgnoreFields
        
        // Get the spec section which contains the actual data
        def spec = map['spec'] as Map<String, Object>
        if (!spec) {
            return map
        }

        result['version'] = map['version']
        result['kind'] = map['kind']
        
        def normalizedSpec = [:] as Map<String, Object>
        spec.each { String key, Object value ->
            if (key in fieldsToIgnore) {
                return // skip ephemeral fields
            }
            
            if (key == 'path' && original instanceof FileOutput) {
                normalizedSpec[key] = relativizePath(value as String)
            } else if (value instanceof Map) {
                normalizedSpec[key] = normalizeNestedMap(value as Map<String, Object>)
            } else if (value instanceof List) {
                normalizedSpec[key] = normalizeList(value as List)
            } else {
                normalizedSpec[key] = value
            }
        }
        
        result['spec'] = normalizedSpec
        return result
    }

    /**
     * Normalize nested maps recursively
     */
    protected Map<String, Object> normalizeNestedMap(Map<String, Object> map) {
        def result = [:] as Map<String, Object>
        def fieldsToIgnore = EPHEMERAL_FIELDS + additionalIgnoreFields
        
        map.each { String key, Object value ->
            if (key in fieldsToIgnore) {
                return
            }
            if (value instanceof Map) {
                result[key] = normalizeNestedMap(value as Map<String, Object>)
            } else if (value instanceof List) {
                result[key] = normalizeList(value as List)
            } else {
                result[key] = value
            }
        }
        return result
    }

    /**
     * Normalize lists recursively
     */
    protected List normalizeList(List list) {
        return list.collect { item ->
            if (item instanceof Map) {
                return normalizeNestedMap(item as Map<String, Object>)
            } else if (item instanceof List) {
                return normalizeList(item as List)
            } else {
                return item
            }
        }
    }

    /**
     * Relativize an absolute path to the output base
     */
    String relativizePath(String absolutePath) {
        if (!absolutePath || !outputBase) {
            return absolutePath
        }
        
        try {
            def path = Paths.get(absolutePath)
            if (path.isAbsolute() && path.startsWith(outputBase)) {
                return outputBase.relativize(path).toString()
            }
        } catch (Exception e) {
            // If path parsing fails, return as-is
        }
        return absolutePath
    }

    /**
     * Collect and normalize an entire workflow tree including all outputs
     *
     * @param store The lineage store
     * @param workflowLid The LID of the workflow run
     * @return Normalized representation of the entire workflow tree
     */
    Map<String, Object> normalizeWorkflowTree(LinStore store, String workflowLid) {
        def result = [:] as Map<String, Object>
        
        // Strip lid:// prefix if present
        def key = workflowLid.startsWith('lid://') ? workflowLid.substring(6) : workflowLid
        
        // Load and normalize the workflow run
        def workflowRun = store.load(key)
        if (!workflowRun) {
            throw new IllegalArgumentException("Workflow run not found: ${workflowLid}")
        }
        
        result['workflowRun'] = normalize(workflowRun)
        
        // Collect all outputs (sub-keys)
        def outputs = [:] as Map<String, Object>
        store.getSubKeys(key).each { String subKey ->
            def record = store.load(subKey)
            if (record) {
                // Use relative key for output identification
                def relativeKey = subKey.startsWith(key) ? subKey.substring(key.length()) : subKey
                if (relativeKey.startsWith('/')) {
                    relativeKey = relativeKey.substring(1)
                }
                outputs[relativeKey] = normalize(record)
            }
        }
        
        if (outputs) {
            result['outputs'] = outputs
        }
        
        return result
    }

    /**
     * Compare two normalized workflow trees and return differences
     *
     * @param tree1 First normalized tree
     * @param tree2 Second normalized tree
     * @return Map of differences, empty if trees are equivalent
     */
    static Map<String, Object> compare(Map<String, Object> tree1, Map<String, Object> tree2) {
        def differences = [:] as Map<String, Object>
        
        compareRecursive('', tree1, tree2, differences)
        
        return differences
    }

    /**
     * Recursively compare two maps and collect differences
     */
    protected static void compareRecursive(String path, Object obj1, Object obj2, Map<String, Object> differences) {
        if (obj1 == obj2) {
            return
        }
        
        if (obj1 == null || obj2 == null) {
            differences[path ?: 'root'] = [expected: obj1, actual: obj2]
            return
        }
        
        if (obj1.class != obj2.class) {
            differences[path ?: 'root'] = [expected: obj1, actual: obj2]
            return
        }
        
        if (obj1 instanceof Map) {
            def map1 = obj1 as Map<String, Object>
            def map2 = obj2 as Map<String, Object>
            def allKeys = (map1.keySet() + map2.keySet()) as Set<String>
            
            allKeys.each { key ->
                String newPath = path ? "${path}.${key}".toString() : key
                compareRecursive(newPath, map1[key], map2[key], differences)
            }
        } else if (obj1 instanceof List) {
            def list1 = obj1 as List
            def list2 = obj2 as List
            
            if (list1.size() != list2.size()) {
                differences[path] = [expected: list1, actual: list2, reason: 'size mismatch']
                return
            }
            
            list1.eachWithIndex { item, int idx ->
                String idxPath = "${path}[${idx}]".toString()
                compareRecursive(idxPath, item, list2[idx], differences)
            }
        } else if (obj1 != obj2) {
            differences[path] = [expected: obj1, actual: obj2]
        }
    }
}

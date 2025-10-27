/*
 * Copyright 2013-2024, Seqera Labs
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
 *
 */

package io.seqera.tower.plugin

import groovy.transform.CompileStatic
import nextflow.config.spec.ConfigOption
import nextflow.script.dsl.Description

/**
 * Model Seqera Platform dataset upload configuration
 *
 * @author Edmund Miller <edmund.a.miller@gmail.com>
 */
@CompileStatic
class DatasetConfig {

    @ConfigOption
    @Description("""
        Enable automatic upload of workflow outputs to Seqera Platform datasets (default: `false`).
    """)
    final boolean enabled

    @ConfigOption
    @Description("""
        Dataset creation mode: `auto` to automatically create datasets, `existing` to only use existing datasets (default: `auto`).
    """)
    final String createMode

    @ConfigOption
    @Description("""
        Name pattern for auto-created datasets. Supports variables: `workflow.runName`, `workflow.sessionId` (default: `\${workflow.runName}-outputs`).
    """)
    final String namePattern

    @ConfigOption
    @Description("""
        Per-output dataset configuration. Each output can specify `datasetId` and `enabled` properties.
    """)
    final Map<String, Map> perOutput

    DatasetConfig() {
        this(Collections.emptyMap())
    }

    DatasetConfig(Map opts) {
        this.enabled = opts.enabled != null ? opts.enabled as boolean : false
        this.createMode = opts.createMode as String ?: 'auto'
        this.namePattern = opts.namePattern as String ?: '${workflow.runName}-outputs'
        this.perOutput = opts.perOutput as Map ?: Collections.emptyMap()
    }

    /**
     * Get configuration for a specific output
     *
     * @param outputName The name of the workflow output
     * @return Configuration map for the output, or empty map if not configured
     */
    Map getOutputConfig(String outputName) {
        return perOutput?.get(outputName) as Map ?: Collections.emptyMap()
    }

    /**
     * Check if dataset upload is enabled for a specific output
     *
     * @param outputName The name of the workflow output
     * @return true if enabled, false otherwise
     */
    boolean isEnabledForOutput(String outputName) {
        if (!enabled)
            return false

        final outputConfig = getOutputConfig(outputName)
        if (outputConfig.containsKey('enabled'))
            return outputConfig.enabled as boolean

        return true
    }

    /**
     * Get the dataset ID for a specific output, if configured
     *
     * @param outputName The name of the workflow output
     * @return The dataset ID, or null if not configured
     */
    String getDatasetId(String outputName) {
        final outputConfig = getOutputConfig(outputName)
        return outputConfig.datasetId as String
    }

    /**
     * Check if the configuration allows auto-creating datasets
     *
     * @return true if auto-create is enabled
     */
    boolean isAutoCreateEnabled() {
        return createMode == 'auto'
    }

}

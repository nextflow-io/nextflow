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

package nextflow.lineage.cli

import static nextflow.lineage.fs.LinPath.*

import java.nio.charset.StandardCharsets
import java.nio.file.Path

import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.cli.CmdLineage
import nextflow.config.ConfigMap
import nextflow.exception.AbortOperationException
import nextflow.lineage.LinHistoryRecord
import nextflow.lineage.LinNormalizer
import nextflow.lineage.LinPropertyValidator
import nextflow.lineage.LinStore
import nextflow.lineage.LinStoreFactory
import nextflow.lineage.LinUtils
import nextflow.lineage.fs.LinPathFactory
import nextflow.lineage.serde.LinEncoder
import nextflow.ui.TableBuilder
import org.eclipse.jgit.diff.DiffAlgorithm
import org.eclipse.jgit.diff.DiffFormatter
import org.eclipse.jgit.diff.RawText
import org.eclipse.jgit.diff.RawTextComparator
/**
 * Implements lineage command line operations
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class LinCommandImpl implements CmdLineage.LinCommand {

    private static final Path DEFAULT_HTML_FILE = Path.of("lineage.html")

    private static final String ERR_NOT_LOADED = 'Error lineage store not loaded - Check Nextflow configuration'

    @Override
    void list(ConfigMap config) {
        final session = new Session(config)
        final store = LinStoreFactory.getOrCreate(session)
        if (store) {
            printHistory(store)
        } else {
            println ERR_NOT_LOADED
        }
    }

    private void printHistory(LinStore store) {
        final records = store.historyLog?.records
        if( !records ) {
            println("No workflow runs found in lineage history log")
            return
        }
        def table = new TableBuilder(cellSeparator: '\t')
            .head('TIMESTAMP')
            .head('RUN NAME')
            .head('SESSION ID')
            .head('LINEAGE ID')
        for (LinHistoryRecord record : records) {
            table.append(record.toList())
        }
        println table.toString()
    }

    @Override
    void view(ConfigMap config, List<String> args) {
        if( !isLidUri(args[0]) )
            throw new Exception("Identifier is not a lineage URL")
        final store = LinStoreFactory.getOrCreate(new Session(config))
        if ( !store ) {
            println ERR_NOT_LOADED
            return
        }
        try {
            def entry = LinUtils.getMetadataObject(store, new URI(args[0]))
            if( entry == null ) {
                println "No entry found for ${args[0]}"
                return
            }
            println LinUtils.encodeSearchOutputs(entry, true)
        } catch (Throwable e) {
            println "Error loading ${args[0]} - ${e.message}"
        }
    }

    @Override
    void render(ConfigMap config, List<String> args) {
        final store = LinStoreFactory.getOrCreate(new Session(config))
        if( !store ) {
            println ERR_NOT_LOADED
            return
        }
        try {
            final renderFile = args.size() > 1 ? Path.of(args[1]) : DEFAULT_HTML_FILE
            new LinDagRenderer(store).render(args[0], renderFile)
            println("Rendered lineage graph for ${args[0]} to $renderFile")
        } catch (Throwable e) {
            println("ERROR: rendering lineage graph - ${e.message}")
        }
    }

    @Override
    void diff(ConfigMap config, List<String> args) {
        if (!isLidUri(args[0]) || !isLidUri(args[1]))
            throw new Exception("Identifier is not a lineage URL")

        final store = LinStoreFactory.getOrCreate(new Session(config))
        if (!store) {
            println ERR_NOT_LOADED
            return
        }
        try {
            final key1 = args[0].substring(LID_PROT.size())
            final entry1 = store.load(key1)
            if (!entry1) {
                println "No entry found for ${args[0]}."
                return
            }
            final key2 = args[1].substring(LID_PROT.size())
            final entry2 = store.load(key2)
            if (!entry2) {
                println "No entry found for ${args[1]}."
                return
            }
            final encoder = new LinEncoder().withPrettyPrint(true)
            generateDiff(encoder.encode(entry1), key1, encoder.encode(entry2), key2)
        } catch (Throwable e) {
            println "Error generating diff between ${args[0]}: $e.message"
        }
    }

    private static void generateDiff(String entry1, String key1, String entry2, String key2) {
        // Convert strings to JGit RawText format
        final text1 = new RawText(entry1.getBytes(StandardCharsets.UTF_8))
        final text2 = new RawText(entry2.getBytes(StandardCharsets.UTF_8))

        // Set up the diff algorithm (Git-style diff)
        final diffAlgorithm = DiffAlgorithm.getAlgorithm(DiffAlgorithm.SupportedAlgorithm.MYERS)
        final diffComparator = RawTextComparator.DEFAULT

        // Compute the differences
        final editList = diffAlgorithm.diff(diffComparator, text1, text2)

        final output = new StringBuilder()
        // Add header
        output.append("diff --git ${key1} ${key2}\n")
        output.append("--- ${key1}\n")
        output.append("+++ ${key2}\n")

        // Use DiffFormatter to display results in Git-style format
        final outputStream = new ByteArrayOutputStream()
        final diffFormatter = new DiffFormatter(outputStream)
        diffFormatter.setOldPrefix(key1)
        diffFormatter.setNewPrefix(key2)
        diffFormatter.format(editList, text1, text2)
        output.append(outputStream.toString(StandardCharsets.UTF_8))

        println output.toString()
    }

    @Override
    void find(ConfigMap config, List<String> args) {
        final store = LinStoreFactory.getOrCreate(new Session(config))
        if (!store) {
            println ERR_NOT_LOADED
            return
        }
        try {
            final params = parseFindArgs(args)
            new LinPropertyValidator().validateQueryParams(params.keySet())
            try (def stream = store.search(params) ) {
                stream.forEach { println asUriString(it) }
            }
        } catch (Throwable e){
            println "Error searching for ${args[0]}. ${e.message}"
        }
    }

    @Override
    void check(ConfigMap config, List<String> args) {
        final store = LinStoreFactory.getOrCreate(new Session(config))
        if (!store) {
            println ERR_NOT_LOADED
            return
        }
        final valid = LinPathFactory.create(args[0]).validate()
        if( !valid )
            throw new AbortOperationException(valid.error)
        else
            println("Checksum validation succeed")
    }

    private Map<String, List<String>> parseFindArgs(List<String> args){
        Map<String, List<String>> params = [:].withDefault { [] }

        args.collectEntries { pair ->
            final idx = pair.indexOf('=')
            if( idx < 0 )
                throw new IllegalArgumentException("Parameter $pair doesn't contain '=' separator")
            def key = URLDecoder.decode(pair[0..<idx], 'UTF-8')
            final value = URLDecoder.decode(pair[(idx + 1)..<pair.length()], 'UTF-8')

            // Convert 'label' key in the CLI to 'labels' property in the model
            if( key == 'label' )
                key = 'labels'

            params[key] << value
        }
        return params
    }

    @Override
    void validate(ConfigMap config, List<String> args) {
        final store = LinStoreFactory.getOrCreate(new Session(config))
        if (!store) {
            println ERR_NOT_LOADED
            return
        }

        try {
            // Parse arguments
            final parsedArgs = parseValidateArgs(args)
            final String newLid = parsedArgs.newLid as String
            final String baselineLid = parsedArgs.baselineLid as String
            final List<String> ignoreFields = parsedArgs.ignoreFields as List<String>
            final String outputBase = parsedArgs.outputBase as String

            if (!isLidUri(newLid) || !isLidUri(baselineLid)) {
                throw new AbortOperationException("Both arguments must be lineage URLs (lid://...)")
            }

            // Create normalizer with configuration
            final normalizer = new LinNormalizer()
            if (outputBase) {
                normalizer.withOutputBase(outputBase)
            }
            if (ignoreFields) {
                normalizer.withIgnoreFields(ignoreFields)
            }

            // Normalize both workflow trees
            final newTree = normalizer.normalizeWorkflowTree(store, newLid)
            final baselineTree = normalizer.normalizeWorkflowTree(store, baselineLid)

            // Compare
            final differences = LinNormalizer.compare(newTree, baselineTree)

            if (differences.isEmpty()) {
                println "Workflow runs are semantically equivalent"
            } else {
                println "Workflow runs differ:"
                println ""
                
                // Generate a readable diff using JGit
                final String newJson = JsonOutput.prettyPrint(JsonOutput.toJson(newTree))
                final String baselineJson = JsonOutput.prettyPrint(JsonOutput.toJson(baselineTree))
                
                final String newKey = newLid.substring(LID_PROT.size())
                final String baselineKey = baselineLid.substring(LID_PROT.size())
                generateDiff(baselineJson, baselineKey, newJson, newKey)
                
                throw new AbortOperationException("Validation failed: workflow runs are not equivalent")
            }
        } catch (AbortOperationException e) {
            throw e
        } catch (IllegalArgumentException e) {
            throw new AbortOperationException(e.message)
        } catch (Throwable e) {
            throw new AbortOperationException("Error validating workflow runs: ${e.message}")
        }
    }

    /**
     * Parse validate command arguments
     */
    private Map parseValidateArgs(List<String> args) {
        String newLid = null
        String baselineLid = null
        List<String> ignoreFields = []
        String outputBase = null

        def iter = args.iterator()
        while (iter.hasNext()) {
            def arg = iter.next()
            
            if (arg == '--against') {
                if (!iter.hasNext()) {
                    throw new IllegalArgumentException("--against requires a value")
                }
                baselineLid = iter.next()
            } else if (arg.startsWith('--against=')) {
                baselineLid = arg.substring('--against='.length())
            } else if (arg == '--ignore-fields') {
                if (!iter.hasNext()) {
                    throw new IllegalArgumentException("--ignore-fields requires a value")
                }
                ignoreFields = iter.next().split(',').toList()
            } else if (arg.startsWith('--ignore-fields=')) {
                ignoreFields = arg.substring('--ignore-fields='.length()).split(',').toList()
            } else if (arg == '--output-base') {
                if (!iter.hasNext()) {
                    throw new IllegalArgumentException("--output-base requires a value")
                }
                outputBase = iter.next()
            } else if (arg.startsWith('--output-base=')) {
                outputBase = arg.substring('--output-base='.length())
            } else if (!arg.startsWith('-') && !newLid) {
                newLid = arg
            } else if (!arg.startsWith('-')) {
                throw new IllegalArgumentException("Unexpected argument: ${arg}")
            }
        }

        if (!newLid) {
            throw new IllegalArgumentException("Missing workflow run LID")
        }
        if (!baselineLid) {
            throw new IllegalArgumentException("Missing --against baseline LID")
        }

        return [
            newLid: newLid,
            baselineLid: baselineLid,
            ignoreFields: ignoreFields,
            outputBase: outputBase
        ]
    }
}

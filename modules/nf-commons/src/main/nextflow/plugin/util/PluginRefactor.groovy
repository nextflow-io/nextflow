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

package nextflow.plugin.util

import java.util.regex.Matcher
import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import nextflow.extension.FilesEx

@Slf4j
@CompileStatic
class PluginRefactor {

    private String pluginName

    private String providerName

    private String packageName

    private File pluginDir

    private String pluginClassPrefix

    private File gradleSettingsFile

    private File gradleBuildFile

    private Map<String,String> tokenMapping = new HashMap<>()

    String getPluginName() {
        return pluginName
    }

    String getPackageName() {
        return packageName
    }

    String getProviderName() {
        return providerName
    }

    File getPluginDir() {
        return pluginDir
    }

    PluginRefactor withPluginDir(File directory) {
        this.pluginDir = directory.absoluteFile.canonicalFile
        return this
    }

    PluginRefactor withPluginName(String name) {
        this.pluginName = normalizeToKebabCase(name)
        this.pluginClassPrefix = normalizeToClassName(name)
        if( pluginName.toLowerCase()=='plugin' )
            throw new IllegalStateException("Invalid plugin name: '$name'")
        if( !pluginClassPrefix )
            throw new IllegalStateException("Invalid plugin name: '$name'")
        return this
    }

    PluginRefactor withProviderName(String name) {
        this.providerName = name?.trim()
        if( !providerName )
            throw new AbortOperationException("Provider name cannot be empty or blank")
        this.packageName = normalizeToPackageNameSegment(name)
        if( !packageName )
            throw new AbortOperationException("Invalid provider name: $name")
        return this
    }

    protected void init() {
        if( !pluginName )
            throw new IllegalStateException("Missing plugin name")
        if( !providerName )
            throw new IllegalStateException("Missing provider name")
        // initial
        this.gradleBuildFile = new File(pluginDir, 'build.gradle')
        this.gradleSettingsFile = new File(pluginDir, 'settings.gradle')
        if( !gradleBuildFile.exists() )
            throw new AbortOperationException("Plugin file does not exist: $gradleBuildFile")
        if( !gradleSettingsFile.exists() )
            throw new AbortOperationException("Plugin file does not exist: $gradleSettingsFile")
        if( !providerName )
            throw new AbortOperationException("Plugin provider name is missing")
        if( !packageName )
            throw new AbortOperationException("Plugin package name is missing")
        // packages to be updates
        tokenMapping.put('acme', packageName)
        tokenMapping.put('provider-name', providerName)
        tokenMapping.put('nf-plugin-template', pluginName)
    }

    void apply() {
        init()
        replacePrefixInFiles(pluginDir, pluginClassPrefix)
        renameDirectory(new File(pluginDir, "src/main/groovy/acme"), new File(pluginDir, "src/main/groovy/${packageName}"))
        renameDirectory(new File(pluginDir, "src/test/groovy/acme"), new File(pluginDir, "src/test/groovy/${packageName}"))
        updateClassNamesAndSymbols(pluginDir)
    }

    protected void replacePrefixInFiles(File rootDir, String newPrefix) {
        if (!rootDir.exists() || !rootDir.isDirectory()) {
            throw new IllegalStateException("Invalid directory: $rootDir")
        }

        rootDir.eachFileRecurse { file ->
            if (file.isFile() && file.name.startsWith('My') && FilesEx.getExtension(file) in ['groovy']) {
                final newName = file.name.replaceFirst(/^My/, newPrefix)
                final renamedFile = new File(file.parentFile, newName)
                if (file.renameTo(renamedFile)) {
                    log.debug "Renamed: ${file.name} -> ${renamedFile.name}"
                    final source = FilesEx.getBaseName(file)
                    final target = FilesEx.getBaseName(renamedFile)
                    tokenMapping.put(source, target)
                }
                else {
                    throw new IllegalStateException("Failed to rename: ${file.name}")
                }
            }
        }
    }

    protected void updateClassNamesAndSymbols(File rootDir) {
        rootDir.eachFileRecurse { file ->
            if (file.isFile() && FilesEx.getExtension(file) in ['groovy','gradle','md']) {
                replaceTokensInFile(file, tokenMapping)
            }
        }
    }

    protected void replaceTokensInFile(File inputFile, Map<String, String> replacements, File outputFile = inputFile) {
        def content = inputFile.text

        // Replace each key with its corresponding value
        for( Map.Entry<String,String> entry : replacements ) {
            content = content.replaceAll(Pattern.quote(entry.key), Matcher.quoteReplacement(entry.value))
        }

        outputFile.text = content
        log.debug "Replacements done in: ${outputFile.path}"
    }

    protected void renameDirectory(File oldDir, File newDir) {
        if (!oldDir.exists() || !oldDir.isDirectory()) {
            throw new AbortOperationException("Plugin template directory to rename does not exist: $oldDir")
        }

        if( oldDir==newDir ) {
            log.debug "Unneeded path rename: $oldDir -> $newDir"
        }

        if (newDir.exists()) {
            throw new AbortOperationException("Plugin target directory already exists: $newDir")
        }

        if (oldDir.renameTo(newDir)) {
            log.debug "Successfully renamed: $oldDir -> $newDir"
        }
        else {
            throw new AbortOperationException("Unable to replace plugin template path: $oldDir -> $newDir")
        }
    }

    static String normalizeToClassName(String input) {
        // Replace non-alphanumeric characters with spaces (except underscores)
        final cleaned = input.replaceAll(/[^a-zA-Z0-9_]/, ' ')
            .replaceAll(/_/, ' ')
            .trim()
        // Split by whitespace, capitalize each word, join them
        final parts = cleaned.split(/\s+/).collect { it.capitalize() }
        final result = parts.join('').replace('Plugin','')
        // Remove "Nf" prefix only
        return result.startsWith('Nf') ? result.substring(2) : result
    }

    static String normalizeToKebabCase(String input) {
        // Insert spaces before capital letters (handles CamelCase)
        def spaced = input.replaceAll(/([a-z])([A-Z])/, '$1 $2')
            .replaceAll(/([A-Z]+)([A-Z][a-z])/, '$1 $2')
        // Replace non-alphanumeric characters and underscores with spaces
        def cleaned = spaced.replaceAll(/[^a-zA-Z0-9]/, ' ')
            .trim()
        // Split, lowercase, and join with hyphens
        def parts = cleaned.split(/\s+/).collect { it.toLowerCase() }
        return parts.join('-')
    }

    static String normalizeToPackageNameSegment(String input) {
        // Replace non-alphanumeric characters with spaces
        def cleaned = input.replaceAll(/[^a-zA-Z0-9]/, ' ')
            .trim()
        // Split into lowercase words and join
        def parts = cleaned.split(/\s+/).collect { it.toLowerCase() }
        def name = parts.join('')

        // Strip leading digits
        name = name.replaceFirst(/^\d+/, '')
        return name ?: null
    }

}

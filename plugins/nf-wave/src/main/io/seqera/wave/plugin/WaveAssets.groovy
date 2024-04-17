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

package io.seqera.wave.plugin


import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import io.seqera.wave.api.PackagesSpec
import nextflow.script.bundle.ResourcesBundle
import nextflow.util.CacheHelper
import nextflow.util.StringUtils
/**
 * Hold assets required to fulfill wave container image build
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Canonical
@CompileStatic
class WaveAssets {
    final String containerImage
    final String containerPlatform
    final ResourcesBundle moduleResources
    final ContainerConfig containerConfig
    final String containerFile
    final PackagesSpec packagesSpec
    final ResourcesBundle projectResources
    final boolean singularity

    static fromImage(String containerImage,String containerPlatform=null) {
        new WaveAssets(containerImage, containerPlatform)
    }

    static fromDockerfile(String dockerfile, String containerPlatform=null) {
        new WaveAssets(null, containerPlatform, null, null, dockerfile)
    }

    String dockerFileEncoded() {
        return containerFile
                ? containerFile.bytes.encodeBase64()
                : null
    }

    @Memoized
    String fingerprint() {
        final allMeta = new ArrayList(10)
        allMeta.add( this.containerImage )
        allMeta.add( this.moduleResources?.fingerprint() )
        allMeta.add( this.containerConfig?.fingerprint() )
        allMeta.add( this.containerFile )
        allMeta.add( this.packagesSpec ? fingerprint(this.packagesSpec) : null )
        allMeta.add( this.projectResources?.fingerprint() )
        allMeta.add( this.containerPlatform )
        return CacheHelper.hasher(allMeta).hash().toString()
    }

    protected String fingerprint(PackagesSpec spec) {
        final allMeta = new ArrayList(10)
        allMeta.add( spec.type.toString() )
        allMeta.add( spec.environment )
        allMeta.add( spec.entries )
        allMeta.add( spec.condaOpts?.mambaImage )
        allMeta.add( spec.condaOpts?.commands )
        allMeta.add( spec.condaOpts?.basePackages )
        allMeta.add( spec.spackOpts?.commands )
        allMeta.add( spec.spackOpts?.basePackages )
        allMeta.add( spec.channels )
        return CacheHelper.hasher(allMeta).hash().toString()
    }

    static void validateContainerName(String name) {
        if( !name )
            return
        final scheme = StringUtils.getUrlProtocol(name)
        if( scheme )
            throw new IllegalArgumentException("Wave container request image cannot start with URL like prefix - offending value: $name")
    }
}

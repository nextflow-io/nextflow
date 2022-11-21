/*
 * Copyright 2020-2022, Seqera Labs
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

import java.nio.file.Path

import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import nextflow.script.bundle.ResourcesBundle
import nextflow.util.CacheHelper
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
    final String dockerFileContent
    final Path condaFile
    final ResourcesBundle projectResources

    static fromImage(String containerImage,String containerPlatform=null) {
        new WaveAssets(containerImage, containerPlatform)
    }

    static fromDockerfile(String dockerfile, String containerPlatform=null) {
        new WaveAssets(null, containerPlatform, null, null, dockerfile)
    }

    String dockerFileEncoded() {
        return dockerFileContent
                ? dockerFileContent.bytes.encodeBase64()
                : null
    }

    String condaFileEncoded() {
        return condaFile
                ? condaFile.text.bytes.encodeBase64()
                : null
    }

    @Memoized
    String fingerprint() {
        final allMeta = new ArrayList(10)
        allMeta.add( this.containerImage )
        allMeta.add( this.moduleResources?.fingerprint() )
        allMeta.add( this.containerConfig?.fingerprint() )
        allMeta.add( this.dockerFileContent )
        allMeta.add( this.condaFile )
        allMeta.add( this.projectResources?.fingerprint() )
        allMeta.add( this.containerPlatform )
        return CacheHelper.hasher(allMeta).hash().toString()
    }
}

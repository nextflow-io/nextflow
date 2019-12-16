/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.container
import java.nio.file.Path
import java.nio.file.Paths

import groovy.transform.PackageScope
import nextflow.util.Escape

/**
 * Helper class to normalise a container image name depending
 * the the current select container engine
 *
 * @author Emilio Palumbo <emilio.palumbo@crg.eu>
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ContainerHandler {

    final private static Path CWD = Paths.get('.').toAbsolutePath()

    @PackageScope ContainerConfig config

    @PackageScope Path baseDir

    ContainerHandler(Map containerConfig) {
        this(containerConfig, CWD)
    }

    ContainerHandler(Map containerConfig, Path dir) {
        this.config = containerConfig as ContainerConfig
        this.baseDir = dir
    }

    String normalizeImageName(String imageName) {
        final engine = config.getEngine()
        if( engine == 'shifter' ) {
            normalizeShifterImageName(imageName)
        }
        else if( engine == 'udocker' ) {
            normalizeUdockerImageName(imageName)
        }
        else if( engine == 'singularity' ) {
            final normalizedImageName = normalizeSingularityImageName(imageName)
            if( !config.isEnabled() || !normalizedImageName )
                return normalizedImageName
            final formats = ['docker', 'docker-daemon', 'shub', 'library']
            final requiresCaching =  formats.any { normalizedImageName.startsWith(it) }
            final result = requiresCaching ? createCache(this.config, normalizedImageName) : normalizedImageName
            Escape.path(result)
        }
        else {
            normalizeDockerImageName(imageName)
        }
    }

    @PackageScope
    String createCache(Map config, String imageName) {
        new SingularityCache(new ContainerConfig(config)) .getCachePathFor(imageName) .toString()
    }

    /**
     * Normalize Shifter image name adding `docker:` prefix or `:latest`
     * when required
     *
     * @param imageName The container image name
     * @return Image name in Shifter canonical format
     */
     @PackageScope
     String normalizeShifterImageName( String imageName ) {

        if( !imageName )
            return null

        def items = imageName.tokenize(':')
        if( items.size()==3 ) {
            // it is in the canonical form i.e. `type:image:tag`
            return imageName
        }

        if( items.size()==1 ) {
            return "docker:$imageName:latest"
        }

        return !imageName.startsWith("docker:") ? "docker:$imageName" : "$imageName:latest"
    }

    /**
     * Normalize Docker image name adding the docker registry
     * when required
     *
     * @param imageName The container image name
     * @return Image name in Docker canonical format
     */
     @PackageScope
     String normalizeDockerImageName( String imageName) {

        if( !imageName )
            return null

        String reg = this.config?.registry
        if( !reg )
            return imageName

        if( isAbsoluteDockerName(imageName) )
            return imageName

        if( !reg.endsWith('/') )
            reg += '/'

        return reg + imageName
    }

     static boolean isAbsoluteDockerName(String image) {
        def p = image.indexOf('/')
        if( p==-1 )
            return false

        image = image.substring(0,p)
        image.contains('.') || image.contains(':')
    }

    /**
     * Normalize Udocker image name adding `:latest`
     * when required
     *
     * @param imageName The container image name
     * @return Image name in Udocker canonical format
     */
     @PackageScope
     String normalizeUdockerImageName( String imageName ) {

        if( !imageName )
            return null

        if( !imageName.contains(':') )
            imageName += ':latest'

        return imageName
    }

    /**
     * Normalize Singularity image name resolving the absolute path or
     * adding `docker://` prefix when required
     *
     * @param imageName The container image name
     * @return Image name in Singularity canonical format
     */
     @PackageScope
     String normalizeSingularityImageName(String img) {
        if( !img )
            return null

        // when starts with `/` it's an absolute image file path, just return it
        if( img.startsWith("/") )
            return img

         // when starts with `file://` it's an image file path, resolve it against the current path
        if (img.startsWith("file://")) {
             return baseDir.resolve(img.substring(7)).toString()
        }

        // check if matches a protocol scheme such as `docker://xxx`
        if( img =~ '^[^/:\\. ]+://(.*)' ) {
            return img
        }

        // if it's the path of an existing image file return it
        def imagePath = baseDir.resolve(img)
        if( imagePath.exists() ) {
            return imagePath.toString()
        }

        // in all other case it's supposed to be the name of an image in the docker hub
        // prefix it with the `docker://` pseudo protocol used by singularity to download it
        return "docker://${img}"
    }
}

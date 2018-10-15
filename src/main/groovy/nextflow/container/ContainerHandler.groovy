/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.container
import java.nio.file.Path
import java.nio.file.Paths

import groovy.transform.PackageScope
/**
 * Helper class to normalise a container image name depending
 * the the current select container engine
 *
 * @author Emilio Palumbo <emiliopalumbo@gmail.com>
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
            final formats = ['docker', 'shub', 'library']
            final requiresCaching =  formats.any { normalizedImageName.startsWith(it) }
            return requiresCaching ? createCache(this.config, normalizedImageName) : normalizedImageName
        }
        else {
            normalizeDockerImageName(imageName)
        }
    }

    @PackageScope
    String createCache(Map config, String imageName) {
        new SingularityCache(config) .getCachePathFor(imageName) .toString()
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

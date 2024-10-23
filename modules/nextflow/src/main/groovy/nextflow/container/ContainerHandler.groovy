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
 */

package nextflow.container

import java.nio.file.Path
import java.nio.file.Paths
import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.container.inspect.ContainerInspectMode
import nextflow.util.Escape
/**
 * Helper class to normalise a container image name depending
 * the the current select container engine
 *
 * @author Emilio Palumbo <emilio.palumbo@crg.eu>
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ContainerHandler {

    final private static Path CWD = Paths.get('.').toAbsolutePath()

    private ContainerConfig config

    private Path baseDir

    ContainerHandler(Map containerConfig) {
        this(containerConfig, CWD)
    }

    ContainerHandler(Map containerConfig, Path dir) {
        this.config = containerConfig as ContainerConfig
        this.baseDir = dir
    }

    ContainerConfig getConfig() { config }
    
    Path getBaseDir() { baseDir }

    String normalizeImageName(String imageName) {
        final engine = config.getEngine()
        if( engine == 'shifter' ) {
            return normalizeShifterImageName(imageName)
        }
        if( engine == 'udocker' ) {
            return normalizeUdockerImageName(imageName)
        }
        if( engine == 'singularity' ) {
            final normalizedImageName = normalizeSingularityImageName(imageName)
            if( !config.isEnabled() || !normalizedImageName )
                return normalizedImageName
            if( normalizedImageName.startsWith('docker://') && config.canRunOciImage() )
                return normalizedImageName
            final requiresCaching = normalizedImageName =~ IMAGE_URL_PREFIX
            if( ContainerInspectMode.active() && requiresCaching )
                return imageName
            final result = requiresCaching ? createSingularityCache(this.config, normalizedImageName) : normalizedImageName
            return Escape.path(result)
        }
        if( engine == 'apptainer' ) {
            final normalizedImageName = normalizeApptainerImageName(imageName)
            if( !config.isEnabled() || !normalizedImageName )
                return normalizedImageName
            if( normalizedImageName.startsWith('docker://') && config.canRunOciImage() )
                return normalizedImageName
            final requiresCaching = normalizedImageName =~ IMAGE_URL_PREFIX
            if( ContainerInspectMode.active() && requiresCaching )
                return imageName
            final result = requiresCaching ? createApptainerCache(this.config, normalizedImageName) : normalizedImageName
            return Escape.path(result)
        }
        if( engine == 'charliecloud' ) {
            final normalizedImageName = normalizeCharliecloudImageName(imageName)
            if( !config.isEnabled() || !normalizedImageName )
                return normalizedImageName
            // if the imagename starts with '/' it's an absolute path
            // otherwise we assume it's in a remote registry and pull it from there
            final requiresCaching = !imageName.startsWith('/')
            if( ContainerInspectMode.active() && requiresCaching )
                return imageName
            final result = requiresCaching ? createCharliecloudCache(this.config, normalizedImageName) : normalizedImageName
            return Escape.path(result)
        }
        // fallback to docker
        return normalizeDockerImageName(imageName)
    }

    @PackageScope
    String createSingularityCache(Map config, String imageName) {
        new SingularityCache(new ContainerConfig(config)) .getCachePathFor(imageName) .toString()
    }

    @PackageScope
    String createApptainerCache(Map config, String imageName) {
        new ApptainerCache(new ContainerConfig(config)) .getCachePathFor(imageName) .toString()
    }

    @PackageScope
    String createCharliecloudCache(Map config, String imageName) {
        new CharliecloudCache(new ContainerConfig(config)) .getCachePathFor(imageName) .toString()
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
     String normalizeDockerImageName(String imageName) {

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


    public static final Pattern IMAGE_URL_PREFIX = ~/^[^\/:. ]+:\/\/(.*)/

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
        if( img =~ IMAGE_URL_PREFIX ) {
            return img
        }

        // if it's the path of an existing image file return it
        def imagePath = baseDir.resolve(img)
        if( imagePath.exists() ) {
            return imagePath.toString()
        }

        // in all other case it's supposed to be the name of an image in the docker hub
        // prefix it with the `docker://` pseudo protocol used by singularity to download it
        return "docker://${normalizeDockerImageName(img)}"
    }

    /**
     * Normalize Apptainer image name resolving the absolute path or
     * adding `docker://` prefix when required
     *
     * @param imageName The container image name
     * @return Image name in Apptainer canonical format
     */
     @PackageScope
     String normalizeApptainerImageName(String img) {
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
        if( img =~ IMAGE_URL_PREFIX ) {
            return img
        }

        // if it's the path of an existing image file return it
        def imagePath = baseDir.resolve(img)
        if( imagePath.exists() ) {
            return imagePath.toString()
        }

        // in all other case it's supposed to be the name of an image in the docker hub
        // prefix it with the `docker://` pseudo protocol used by apptainer to download it
        return "docker://${normalizeDockerImageName(img)}"
    }

    /**
     * Normalize charliecloud image name resolving the absolute path
     *
     * @param imageName The container image name
     * @return Image name in canonical format
     */
     @PackageScope
     String normalizeCharliecloudImageName(String img) {
        if( !img )
            return null

        // when starts with `/` it's an absolute image file path, just return it
        if( img.startsWith("/") ) {
            return img
        }
        // remove docker:// if present
        if( img.startsWith("docker://") ) {
            img = img.minus("docker://")
        }
        // if no tag, add :latest
        if( !img.contains(':') ) {
            img += ':latest'
        }

        // if it's the path of an existing image file return it
        def imagePath = baseDir.resolve(img)
        if( imagePath.exists() ) {
            return imagePath.toString()
        }

        // in all other case it's supposed to be the name of an image
        return "${normalizeDockerImageName(img)}"
    }
}
package nextflow.container
import java.nio.file.Paths

import groovy.transform.PackageScope

/**
 * Helper class to normalise a container image name depending
 * the the current select container engine
 *
 * @author Emilio Palumbo <emiliopalumbo@gmail.com>
 */
class ContainerHandler {

    private ContainerConfig cfg
    private String baseDir

    ContainerHandler(Map containerConfig) {
        this(containerConfig, null)
    }

    ContainerHandler(Map containerConfig, String dir) {
        this.cfg = containerConfig as ContainerConfig
        this.baseDir = dir
    }

    String normalizeImageName(String imageName) {
        final engine = cfg.getEngine()
        if( engine == 'shifter' ) {
            normalizeShifterImageName(imageName)
        }
        else if( engine == 'udocker' ) {
            normalizeUdockerImageName(imageName)
        }
        else if( engine == 'singularity' ) {
            def normalizedImageName = normalizeSingularityImageName(imageName)
            if (normalizedImageName && (normalizedImageName.startsWith("docker://") || normalizedImageName.startsWith("shub://"))) {
                return createCache(this.cfg, normalizedImageName)
            }
            return normalizedImageName
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

        String reg = this.cfg?.registry
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

        if (img.startsWith("/") || img.startsWith("docker://") || img.startsWith("shub://")) {
            return img
        }

        if (img.startsWith("file://")) {
            return Paths.get(new URI(img)).toAbsolutePath().toString()
        }

        def imagePath = Paths.get(baseDir, img)

        if( imagePath.exists() ) {
            return imagePath.toAbsolutePath().toString()
        }

        img = "docker://${img}"

        return img
    }
}

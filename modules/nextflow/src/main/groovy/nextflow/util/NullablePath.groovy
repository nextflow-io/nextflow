package nextflow.util

import groovy.transform.EqualsAndHashCode
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import java.nio.file.Path

/**
 * @author : jorge <jorge.aguilera@seqera.io>
 *
 */
@EqualsAndHashCode
@Slf4j
class NullablePath implements Path {

    @PackageScope
    @Delegate
    Path delegate

    NullablePath(String path) {
        delegate = of(path)
    }
}

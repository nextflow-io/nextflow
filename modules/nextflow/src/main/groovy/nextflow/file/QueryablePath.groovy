package nextflow.file

import java.nio.file.Path

/**
 * Interface to indicate a Path could contain a query that is resolved to several real paths.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
interface QueryablePath {
    boolean hasQuery();
    List<Path> resolveQuery();
}
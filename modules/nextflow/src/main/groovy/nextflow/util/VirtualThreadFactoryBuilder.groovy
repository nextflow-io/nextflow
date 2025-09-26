package nextflow.util

import java.util.concurrent.ThreadFactory

/**
 * Builder class to create a named virtual thread factory.
 * This class is required to avoid failures when using with Java version < 21
 *
 * @author jorge.ejarque@seqera.io
 */
class VirtualThreadFactoryBuilder {

        static ThreadFactory create(String name){
            return Thread.ofVirtual().name(name).factory()
        }
}

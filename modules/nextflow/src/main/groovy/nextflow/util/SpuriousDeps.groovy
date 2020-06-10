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

package nextflow.util

import java.lang.reflect.Method

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cli.CmdBase
import static nextflow.Const.S3_UPLOADER_CLASS
/**
 * This class is used to resolve at runtime some spurious dependencies
 * with optional modules
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class SpuriousDeps {

    static CmdBase cmdCloud() {
        try {
            final clazz = Class.forName('nextflow.cli.CmdCloud')
            return (CmdBase)clazz.newInstance()
        }
        catch (ClassNotFoundException e) {
            return null
        }
    }

    static String getS3UploaderScript() {
        try {
            final clazz = Class.forName('nextflow.cloud.aws.batch.S3Helper')
            final m = clazz.getMethod('getUploaderScript')
            return m.invoke(null)
        }
        catch (ClassNotFoundException e) {
            return null
        }
    }

    static void shutdownS3Uploader() {
        if( classWasLoaded(S3_UPLOADER_CLASS) ) {
            log.debug "AWS S3 uploader shutdown"
            final s3 = Class.forName(S3_UPLOADER_CLASS)
            s3.getMethod('shutdownExecutor').invoke(null)
        }
    }

    static private boolean classWasLoaded(String className) {
        Method find = ClassLoader.class.getDeclaredMethod("findLoadedClass", [String.class] as Class[] );
        find.setAccessible(true)
        return find.invoke(ClassLoader.getSystemClassLoader(), className)
    }

}

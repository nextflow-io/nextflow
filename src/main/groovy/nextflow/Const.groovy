/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow
import java.nio.file.Path
import java.nio.file.Paths
import java.text.SimpleDateFormat

import static nextflow.extension.Bolts.DATETIME_FORMAT

/**
 * Application main constants
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class Const {

    static final public transient BOOL_YES = ['true','yes','on']

    static final public transient BOOL_NO = ['false','no','off']

    /**
     * The application main package name
     */
    static public final String MAIN_PACKAGE = Const.name.split('\\.')[0]

    /**
     * The application main name
     */
    static public final String APP_NAME = MAIN_PACKAGE

    /**
     * The application home folder
     */
    static public final Path APP_HOME_DIR = getHomeDir(APP_NAME)

    /**
     * The application version
     */
    static public final String APP_VER = "18.11.0-edge"

    /**
     * The app build time as linux/unix timestamp
     */
    static public final long APP_TIMESTAMP = 1542042294500

    /**
     * The app build number
     */
    static public final int APP_BUILDNUM = 5016


    /**
     * The app build time string relative to UTC timezone
     */
    static public final String APP_TIMESTAMP_UTC = {

        def tz = TimeZone.getTimeZone('UTC')
        def fmt = new SimpleDateFormat(DATETIME_FORMAT)
        fmt.setTimeZone(tz)
        fmt.format(new Date(APP_TIMESTAMP)) + ' ' + tz.getDisplayName( true, TimeZone.SHORT )

    } ()


    /**
     * The app build time string relative to local timezone
     */
    static public final String APP_TIMESTAMP_LOCAL = {

        def tz = TimeZone.getDefault()
        def fmt = new SimpleDateFormat(DATETIME_FORMAT)
        fmt.setTimeZone(tz)
        fmt.format(new Date(APP_TIMESTAMP)) + ' ' + tz.getDisplayName( true, TimeZone.SHORT )

    } ()

    static String deltaLocal() {
        def utc = APP_TIMESTAMP_UTC.split(' ')
        def loc = APP_TIMESTAMP_LOCAL.split(' ')

        if( APP_TIMESTAMP_UTC == APP_TIMESTAMP_LOCAL ) {
            return ''
        }

        def result = utc[0] == loc[0] ? loc[1,-1].join(' ') : loc.join(' ')
        return "($result)"
    }


    private static Path getHomeDir(String appname) {
        def home = System.getenv('NXF_HOME')
        def result = home ? Paths.get(home) : Paths.get(System.getProperty("user.home")).resolve(".$appname")

        if( !result.exists() && !result.mkdir() ) {
            throw new IllegalStateException("Cannot create path '${result}' -- check file system access permission")
        }

        return result
    }

    /*
     * The application 'logo'
     */
    static public final String SPLASH =

"""
      N E X T F L O W
      version ${APP_VER} build ${APP_BUILDNUM}
      last modified ${APP_TIMESTAMP_UTC} ${deltaLocal()}
      cite doi:10.1038/nbt.3820
      http://nextflow.io
"""

    static public final String S3_UPLOADER_CLASS = 'com.upplication.s3fs.S3OutputStream'

    static public final String ROLE_WORKER = 'worker'

    static public final String ROLE_MASTER = 'master'

    static public final String MANIFEST_FILE_NAME = 'nextflow.config'

    static public final String DEFAULT_MAIN_FILE_NAME = 'main.nf'

    static public final String DEFAULT_ORGANIZATION = System.getenv('NXF_ORG') ?: 'nextflow-io'

    static public final String DEFAULT_HUB = System.getenv('NXF_HUB') ?: 'github'

    static public final File DEFAULT_ROOT = System.getenv('NXF_ASSETS') ? new File(System.getenv('NXF_ASSETS')) : Const.APP_HOME_DIR.resolve('assets').toFile()

    static public final String DEFAULT_BRANCH = 'master'

}

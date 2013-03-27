/*
 * Copyright (c) 2012, the authors.
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

package nextflow

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class Const {

    static final File APP_HOME_DIR

    static final File APP_TMP_DIR

    static {
        APP_HOME_DIR = new File( System.getProperty("user.home"), ".${Const.APP_NAME}" )
        if( !APP_HOME_DIR.exists() && !APP_HOME_DIR.mkdir() ) {
            throw new IllegalStateException("Cannot create path '${APP_HOME_DIR}' -- check file system access permission")
        }

        APP_TMP_DIR = new File(APP_HOME_DIR, 'tmp')
        if( !APP_TMP_DIR.exists() && !APP_TMP_DIR.mkdirs()) {
            throw new IllegalStateException("Cannot create path '${APP_TMP_DIR}' -- check file system access permission")
        }
    }


    static final String MAIN_PACKAGE = Const.class.name.split('\\.')[0]

    static final String APP_NAME = MAIN_PACKAGE

    static final String APP_VER = "0.1"

    static final long APP_TIMESTAMP = 1364382426901

    static final int APP_BUILDNUM = 100

    static final String LOGO =

        """
         _  _         _    __ _
        | \\| |_____ _| |_ / _| |_____ __ __
        | .` / -_) \\ /  _|  _| / _ \\ V  V /
        |_|\\_\\___/_\\_\\\\__|_| |_\\___/\\_/\\_/   ver. $APP_VER
        """
        .stripIndent()




    static DATETIME_FORMAT = "dd/MMM/yyyy HH:mm"

    static SHORT_DATETIME_FORMAT = "HH:mm dd/MMM"

    static TIME_FORMAT = "HH:mm:ss"

}

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

package nextflow.util

import groovy.util.logging.Slf4j
import nextflow.Const

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class FileHelper {

    static private Random rndGen = new Random()

    static final char[] ALPHA = ('a'..'z') as char[]

    static final char[] NUMERIC = ('0'..'9') as char[]

    static final char[] ALPHANUM = (('a'..'z')+('0'..'9')) as char[]


    static String normalizePath( String path ) {

        if( path == '~' ) {
            path = System.properties['user.home']
        }
        else if ( path.startsWith('~/') ) {
            path = System.properties['user.home'] + path[1..-1]
        }


        return path
    }


    /**
     * Creates a random string with the number of character specified
     * e.g. {@code s8hm2nxt3}
     *
     * @param len The len of the final random string
     * @param alphabet The set of characters allowed in the random string
     * @return The generated random string
     */
    static String randomString( int len, char[] alphabet ) {
        assert alphabet.size() > 0

        StringBuilder result = new StringBuilder()
        final max = alphabet.size()

        len.times {
            def index = rndGen.nextInt(max)
            result.append( alphabet[index] )
        }

        return result.toString()
    }

    static String randomString( int len ) {
        randomString(len, ALPHA)
    }


    /**
     * The process scratch folder
     * @param seed
     * @return
     */
    static File createScratchDir( final File baseDir = Const.APP_TMP_DIR ) {

        long timestamp = System.currentTimeMillis()
        while( true ) {

            String rnd1 = randomString(2, NUMERIC)
            String rnd2 = randomString(4, ALPHANUM)
            String rnd3 = randomString(4, ALPHANUM)

            File tempDir = new File(baseDir, "$rnd1/$rnd2-$rnd3".toString());

            if (tempDir.mkdirs()) {
                return tempDir;
            }
            else if ( !tempDir.exists() ) {
                // when 'mkdirs' failed because it was unable to create the folder
                // (since it does not exist) throw an exception
                throw new IllegalStateException("Cannot create scratch folder: '${tempDir}' -- verify access permissions" )
            }

            if( System.currentTimeMillis() - timestamp > 1_000 ) {
                throw new IllegalStateException("Unable to create scratch folder: '${tempDir}' after multiple attempts -- verify access permissions" )
            }

            Thread.sleep(50)
        }
    }

    /**
     * Try to create the create specified folder. If the folder already exists 'increment' the folder name
     * by adding +1 to the name itself. For example
     * <p>
     *     {@code /some/path/name1}  become {@code /some/path/name2}
     *
     *
     * @param path
     * @return
     */
    static File tryCreateDir( File path ) {
        assert path

        path = path.canonicalFile
        assert path != new File('/')


        int count=0
        while( true ) {

            if( path.exists() ) {
                def (name, index) = nameParts( path.name )
                path = new File( path.parentFile, "${name}${index+1}")
                count++

                if( count>100 ) {
                    throw new RuntimeException("Unable to create a work directory using path: '$path' -- delete unused directories '${path.parentFile}/${name}*'")
                }
                else if( count == 20 ) {
                    log.warn "Too many temporary work directory -- delete unused directories '${path.parentFile}/${name}*' "
                }

                continue
            }

            if( !path.mkdirs() ) {
                throw new IOException("Unable to create path: ${path} -- Verify the access permission")
            }

            return path

        }

    }


    def static List nameParts( String name ) {
        assert name

        if( name.isLong() )  {
            return ['', name.toLong()]
        }

        def matcher = name =~ /^(\S*)?(\D)(\d+)$/
        if( matcher.matches() ) {
            return [ matcher[0][1] + matcher[0][2], matcher[0][3].toString().toInteger() ]
        }
        else {
            return [name,0]
        }

    }

    /**
     * Check whenever a file or a directory is empty
     *
     * @param file The file path to
     */
    def static boolean isEmpty( File file ) {
        assert file
        if ( file.isDirectory() ) {
            file.list()?.size()==0 }
        else {
            return file.size()==0
        }
    }

    def static isNotEmpty( File path ) {
        !isEmpty(path)
    }

    def delete( File file )  {

        if ( file.isFile() ) {
            file.deleteDir()
        }
    }

}

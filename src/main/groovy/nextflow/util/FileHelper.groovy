/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.spi.FileSystemProvider

import embed.com.google.common.hash.HashCode
import groovy.util.logging.Slf4j
import nextflow.extension.FilesExtensions
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
            path = System.properties['user.home'].toString() + path[1..-1]
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
    @Deprecated
    static File createTempDir( final File baseDir = null ) {

        long timestamp = System.currentTimeMillis()
        while( true ) {

            String rnd1 = randomString(2, NUMERIC)
            String rnd2 = randomString(4, ALPHANUM)
            String rnd3 = randomString(4, ALPHANUM)
            String path = "$rnd1/$rnd2-$rnd3"

            File tempDir = baseDir ? new File(baseDir,path) : new File(path)

            if (tempDir.mkdirs()) {
                return tempDir.absoluteFile;
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



    def static List nameParts( String name ) {
        assert name

        if( name.isLong() )  {
            return ['', name.toLong()]
        }

        def matcher = name =~ /^(\S*)?(\D)(\d+)$/
        if( matcher.matches() ) {
            def entries = (String[])matcher[0]
            return [ entries[1] + entries[2], entries[3].toString().toInteger() ]
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
    def static boolean empty( File file ) {
        assert file

        if( !file.exists() ) {
            return true
        }

        if ( file.isDirectory() ) {
            file.list()?.size()==0
        }
        else {
            file.size()==0
        }
    }

    def static boolean empty( Path path ) {
        if( !Files.exists(path) ) {
            return true
        }

        if ( Files.isDirectory(path) ) {
            def stream = Files.newDirectoryStream(path)
            try {
                Iterator<Path> itr = stream.iterator()
                return !itr.hasNext()
            }
            finally {
                stream.close()
            }
        }
        else {
            Files.size(path)==0
        }
    }

    /**
     * Defines a cacheable path for the specified {@code HashCode} instance.
     *
     * @param hash
     * @return
     */
    final static Path getWorkFolder(Path bashPath, HashCode hash) {
        assert bashPath
        assert hash

        def str = hash.toString()
        def bucket = str.substring(0,2)
        def folder = bashPath.resolve("${bucket}/${str.substring(2)}")

        return folder.toAbsolutePath()
    }

    final static Path createTempFolder(Path basePath) {
        assert basePath

        basePath.delete()
        int count = 0
        while( true ) {
            def hash = CacheHelper.hasher(rndGen.nextLong()).hash().toString()
            def bucket = hash.substring(0,2)
            def result = basePath.resolve( "tmp/$bucket/${hash.substring(2)}" )

            if( FilesExtensions.exists(result) ) {
                if( ++count > 100 ) { throw new IOException("Unable to create a unique temporary path: $result") }
                continue
            }
            if( !FilesExtensions.mkdirs(result) ) {
                throw new IOException("Unable to create temporary parth: $result -- Verify file system access permission")
            }

            return result.toAbsolutePath()
        }

    }

    static Path asPath( String str ) {
        assert str

        int p = str.indexOf('://')
        if( p == -1  ) {
            return Paths.get(str)
        }

        String scheme = str.substring(0, p).trim()
        def provider = getProviderByScheme(scheme)
        if( !provider ) {
            throw new IllegalArgumentException("Unknown file scheme: $scheme");
        }

        return provider.getPath( URI.create(str) )

    }


    private static FileSystemProvider getProviderByScheme( String scheme ) {

        for (FileSystemProvider provider : FileSystemProvider.installedProviders()) {
            if ( scheme == provider.getScheme() ) {
                return provider;
            }
        }
        return null;
    }
}

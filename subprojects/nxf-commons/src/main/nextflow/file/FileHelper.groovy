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

package nextflow.file
import java.lang.reflect.Field
import java.nio.file.FileSystem
import java.nio.file.FileSystemNotFoundException
import java.nio.file.FileSystems
import java.nio.file.FileVisitOption
import java.nio.file.FileVisitResult
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.PathMatcher
import java.nio.file.Paths
import java.nio.file.SimpleFileVisitor
import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.spi.FileSystemProvider
import java.util.concurrent.locks.Lock
import java.util.concurrent.locks.ReentrantLock
import java.util.regex.Pattern

import com.google.common.hash.HashCode
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.extension.Bolts
import nextflow.extension.FilesEx
import nextflow.util.CacheHelper
/**
 * Provides some helper method handling files
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class FileHelper {

    static final Path localTempBasePath

    static final Pattern GLOB_CURLY_BRACKETS = Pattern.compile(/(.*)(\{.+,.+\})(.*)/)

    static final Pattern GLOB_SQUARE_BRACKETS = Pattern.compile(/(.*)(\[.+\])(.*)/)

    static private Random rndGen = new Random()

    static final char[] ALPHA = ('a'..'z') as char[]

    static {
        def temp = System.getenv('NXF_TEMP')
        if( temp ) {
            localTempBasePath = Paths.get(temp)
            log.debug "Using NXF_TEMP=${localTempBasePath}"
            // make sure the path exists
            if( !Files.exists(localTempBasePath) )
                Files.createDirectories(localTempBasePath)
        }
    }

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
        empty(file.toPath())
    }

    def static boolean empty( Path path ) {

        def attrs
        try {
            attrs = Files.readAttributes(path, BasicFileAttributes.class)
        }
        catch (IOException e) {
            return true;
        }

        if ( attrs.isDirectory() ) {
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
            return attrs.size() == 0
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

        int count = 0
        while( true ) {
            def hash = CacheHelper.hasher(rndGen.nextLong()).hash().toString()
            def bucket = hash.substring(0,2)
            def result = basePath.resolve( "tmp/$bucket/${hash.substring(2)}" )

            if( FilesEx.exists(result) ) {
                if( ++count > 100 ) { throw new IOException("Unable to create a unique temporary path: $result") }
                continue
            }
            if( !FilesEx.mkdirs(result) ) {
                throw new IOException("Unable to create temporary parth: $result -- Verify file system access permission")
            }

            return result.toAbsolutePath()
        }

    }


    final static Path createLocalDir(String prefix = 'nxf') {
        if( localTempBasePath )
            Files.createTempDirectory(localTempBasePath, prefix)
        else
            Files.createTempDirectory(prefix)
    }

    static Path asPath( String str ) {
        assert str

        int p = str.indexOf('://')
        if( p == -1 ) {
            return FileSystems.getDefault().getPath(str)
        }

        final uri = URI.create(str)
        if( uri.scheme == 'file' )
            return FileSystems.getDefault().getPath(uri.path)

        getOrCreateFileSystemFor(uri).provider().getPath(uri)
    }

    /**
     * Lazy create and memoize this object since it will never change
     * @return A map holding the current session
     */
    @Memoized
    static protected Map getEnvMap(String scheme) {
        getEnvMap0(scheme, System.getenv())
    }

    @PackageScope
    static Map getEnvMap0(String scheme, Map env) {
        def result = [:]
        if( scheme?.toLowerCase() == 's3' ) {

            List credentials = getAwsCredentials( env, (Map)Global.config )
            if( credentials ) {
                // S3FS expect the access - secret keys pair in lower notation
                result.access_key = credentials[0]
                result.secret_key = credentials[1]
            }
        }
        else {
            assert Global.session, "Session is not available -- make sure to call this after Session object has been created"
            result.session = Global.session
        }
        return result
    }


    /**
     *  Caches File system providers
     */
    @PackageScope
    static final Map<String, FileSystemProvider> providersMap = [:]

    /**
     * Returns a {@link FileSystemProvider} for a file scheme
     * @param scheme A file system scheme e.g. {@code file}, {@code s3}, {@code dxfs}, etc.
     * @return A {@link FileSystemProvider} instance for the specified scheme
     */
    static FileSystemProvider getProviderFor( String scheme ) {

        if( providersMap.containsKey(scheme) )
            return providersMap[scheme]

        synchronized (providersMap) {
            if( providersMap.containsKey(scheme) )
                return providersMap[scheme]

            for (FileSystemProvider provider : FileSystemProvider.installedProviders()) {
                if ( scheme == provider.getScheme() ) {
                    return providersMap[scheme] = provider;
                }
            }

            return null;
        }

    }

    /**
     * Get the instance of the specified {@code FileSystemProvider} class. If the provider is not
     * in the list of installed provided, it creates a new instance and add it to the list
     * <p>
     * This method has been deprecated use {@link #getOrCreateFileSystemFor(java.lang.String)} instead
     *
     * @see {@code FileSystemProvider#installedProviders}
     *
     * @param clazz A class extending {@code FileSystemProvider}
     * @return An instance of the specified class
     */
    @Deprecated
    synchronized static <T extends FileSystemProvider> T getOrInstallProvider( Class<T> clazz ) {

        FileSystemProvider result = FileSystemProvider.installedProviders().find { it.class == clazz }
        if( result )
            return (T)result

        // try to load DnaNexus file system provider dynamically
        result = (T)clazz.newInstance()

        // add it manually
        Field field = FileSystemProvider.class.getDeclaredField('installedProviders')
        field.setAccessible(true)
        List installedProviders = new ArrayList((List)field.get(null))
        installedProviders.add( result )
        field.set(this, Collections.unmodifiableList(installedProviders))
        log.debug "> Added '${clazz.simpleName}' to list of installed providers [${result.scheme}]"
        return (T)result
    }

    private static Lock _fs_lock = new ReentrantLock()

    /**
     * Acquire or create the file system for the given {@link URI}
     *
     * @param uri A {@link URI} locating a file into a file system
     * @param env An option environment specification that may be used to instantiate the underlying file system.
     *          As defined by {@link FileSystemProvider#newFileSystem(java.net.URI, java.util.Map)}
     * @return The corresponding {@link FileSystem} object
     * @throws IllegalArgumentException if does not exist a valid provider for the given URI scheme
     */
    static FileSystem getOrCreateFileSystemFor( URI uri, Map env = null ) {
        assert uri

        /*
         * get the provider for the specified URI
         */
        def provider = getProviderFor(uri.scheme)
        if( !provider )
            throw new IllegalArgumentException("Cannot a find a file system provider for scheme: ${uri.scheme}")

        /*
         * check if already exists a file system for it
         */
        FileSystem fs
        try { fs = provider.getFileSystem(uri) }
        catch( FileSystemNotFoundException e ) { fs=null }
        if( fs )
            return fs

        /*
         * since the file system does not exist, create it a protected block
         */
        Bolts.withLock(_fs_lock) {

            try { fs = provider.getFileSystem(uri) }
            catch( FileSystemNotFoundException e ) { fs=null }
            if( !fs ) {
                log.debug "Creating a file system instance for provider: ${provider.class.simpleName}"
                fs = provider.newFileSystem(uri, env ?: getEnvMap(uri.scheme))
            }
            fs
        }

        return fs
    }

    static FileSystem getOrCreateFileSystemFor( String scheme, Map env = null ) {
        getOrCreateFileSystemFor(URI.create("$scheme:///"), env)
    }

    /**
     * Given a path look for an existing file with that name in the {@code cacheDir} folder.
     * If exists will return the path to it, otherwise make a copy of it in the cache folder
     * and returns to the path to it.
     *
     * @param sourcePath A {@link Path} object to an existing file or directory
     * @param cacheDir A {@link Path} to the folder containing the cached objects
     * @param sessionId An option {@link UUID} used to create the object invalidation key
     * @return A {@link Path} to the cached version of the specified object
     */
    static Path getLocalCachePath( Path sourcePath, Path cacheDir, UUID sessionId = null ) {
        assert sourcePath
        assert cacheDir

        final key = sessionId ? [sessionId, sourcePath] : sourcePath
        final hash = CacheHelper.hasher(key).hash()

        final cached = getWorkFolder(cacheDir, hash).resolve(sourcePath.getFileName().toString())
        if( Files.exists(cached) )
            return cached

        synchronized (cacheDir)  {
            if( Files.exists(cached) )
                return cached

            return FilesEx.copyTo(sourcePath, cached)
        }
    }

    /**
     * Whenever the specified string is a glob file pattern
     *
     * @link  http://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
     *
     * @param filePattern
     * @return
     */
    static isGlobPattern( String filePattern ) {
        assert filePattern
        boolean glob  = false
        glob |= filePattern.contains('*')
        glob |= filePattern.contains('?')
        glob |= GLOB_SQUARE_BRACKETS.matcher(filePattern).matches()
        return glob || GLOB_CURLY_BRACKETS.matcher(filePattern).matches()
    }


    /**
     * Returns a {@code PathMatcher} that performs match operations on the
     * {@code String} representation of {@link Path} objects by interpreting a
     * given pattern.
     *
     * @see FileSystem#getPathMatcher(java.lang.String)
     *
     * @param syntaxAndInput
     * @return
     */
    static PathMatcher getDefaultPathMatcher(String syntaxAndInput) {

        int pos = syntaxAndInput.indexOf(':');
        if (pos <= 0 || pos == syntaxAndInput.length())
            throw new IllegalArgumentException();

        String syntax = syntaxAndInput.substring(0, pos);
        String input = syntaxAndInput.substring(pos+1);

        String expr;
        if (syntax == 'glob') {
            expr = Globs.toUnixRegexPattern(input);
        }
        else if (syntax == 'regex' ) {
            expr = input;
        }
        else {
            throw new UnsupportedOperationException("Syntax '$syntax' not recognized");
        }

        // return matcher
        final Pattern pattern = Pattern.compile(expr);
        return new PathMatcher() {
            @Override
            public boolean matches(Path path) {
                return pattern.matcher(path.toString()).matches();
            }
        };
    }

    /**
     * Get a path matcher for the specified file system and file pattern.
     *
     * It tries to get the matcher returned by {@link FileSystem#getPathMatcher} and
     * falling back to {@link #getDefaultPathMatcher(java.lang.String)} when it return a null or an {@link UnsupportedOperationException}
     *
     * @param fileSystem
     * @param syntaxAndPattern
     * @return A {@link PathMatcher} instance for the specified file pattern
     */
    static PathMatcher getPathMatcherFor(String syntaxAndPattern, FileSystem fileSystem ) {
        assert fileSystem
        assert syntaxAndPattern

        PathMatcher matcher
        try {
            matcher = fileSystem.getPathMatcher(syntaxAndPattern)
        }
        catch( UnsupportedOperationException e ) {
            matcher = null
        }

        if( !matcher ) {
            log.debug "Path matcher not defined by '${fileSystem.class.simpleName}' file system -- using default default strategy"
            matcher = getDefaultPathMatcher(syntaxAndPattern)
        }

        return matcher
    }


    /**
     * Applies the specified action on one or more files and directories matching the specified glob pattern
     *
     * @param folder
     * @param filePattern
     * @param action
     *
     * @throws java.nio.file.NoSuchFileException if the specified {@code folder} does not exist
     */
    static void visitFiles( Map options = null, Path folder, String filePattern, Closure action ) {
        assert folder
        assert filePattern
        assert action

        final type = options?.type ?: 'any'
        final walkOptions = options?.followLinks == false ? EnumSet.noneOf(FileVisitOption.class) : EnumSet.of(FileVisitOption.FOLLOW_LINKS)
        final int maxDepth = getMaxDepth(options?.maxDepth, filePattern)
        final includeHidden = options?.hidden as Boolean ?: filePattern.startsWith('.')
        final includeDir = type in ['dir','any']
        final includeFile = type in ['file','any']
        final syntax = options?.syntax ?: 'glob'
        final relative = options?.relative == true

        final matcher = getPathMatcherFor("$syntax:${folder.resolve(filePattern)}", folder.fileSystem)
        final singleParam = action.getMaximumNumberOfParameters() == 1

        Files.walkFileTree(folder, walkOptions, Integer.MAX_VALUE, new SimpleFileVisitor<Path>() {

            @Override
            public FileVisitResult preVisitDirectory(Path path, BasicFileAttributes attrs) throws IOException {
                int depth = path.nameCount - folder.nameCount
                log.trace "visit dir ($depth) > $path; includeDir: $includeDir; matches: ${matcher.matches(path)}; isDir: ${Files.isDirectory(path)}"

                if (includeDir && matcher.matches(path) && Files.isDirectory(path) && (includeHidden || !Files.isHidden(path))) {
                    def result = relative ? folder.relativize(path) : path
                    singleParam ? action.call(result) : action.call(result,attrs)
                }

                return depth > maxDepth ? FileVisitResult.SKIP_SUBTREE : FileVisitResult.CONTINUE
            }

            @Override
            public FileVisitResult visitFile(Path path, BasicFileAttributes attrs) throws IOException {
                log.trace "visit dir > $path; includeFile: $includeFile; matches: ${matcher.matches(path)}; isDir: ${Files.isDirectory(path)}"

                if (includeFile && matcher.matches(path) && !Files.isDirectory(path) && (includeHidden || !Files.isHidden(path))) {
                    def result = relative ? folder.relativize(path) : path
                    singleParam ? action.call(result) : action.call(result,attrs)
                }

                return FileVisitResult.CONTINUE;
            }
      })

    }

    @PackageScope
    static int getMaxDepth( value, String filePattern ) {

        if( value != null )
            return value as int

        if( filePattern?.contains('**') )
            return Integer.MAX_VALUE

        if( filePattern?.contains('/') )
            return filePattern.split('/').findAll { it }.size()-1

        return 0
    }

    /**
     * Retrieve the AWS credentials from the given context. It look for AWS credential in the following order
     * 1) Nextflow config {@code aws.accessKey} and {@code aws.secretKey} pair
     * 2) System env {@code AWS_ACCESS_KEY} and {@code AWS_SECRET_KEY} pair
     * 3) System env {@code AWS_ACCESS_KEY_ID} and {@code AWS_SECRET_ACCESS_KEY} pair
     *
     *
     * @param env The system environment map
     * @param config The nextflow config object map
     * @return A pair where the first element is the access key and the second the secret key or
     *      {@code null} if the credentials are missing
     */
    static List<String> getAwsCredentials( Map env, Map config ) {

        String a
        String b

        if( config && config.aws instanceof Map ) {
            a = ((Map)config.aws).accessKey
            b = ((Map)config.aws).secretKey

            if( a && b )
                return [a, b]
        }

        if( env && (a=env.AWS_ACCESS_KEY) && (b=env.AWS_SECRET_KEY) ) {
            return [a, b]
        }

        // as define by amazon doc
        // http://docs.aws.amazon.com/cli/latest/userguide/cli-chap-getting-started.html
        if( env && (a=env.AWS_ACCESS_KEY_ID) && (b=env.AWS_SECRET_ACCESS_KEY) )  {
            return [a, b]
        }

        return null
    }

    static List<String> getFolderAndPattern( String filePattern ) {

        def scheme = null;
        int i = filePattern.indexOf('://')
        if( i != -1 ) {
            scheme = filePattern.substring(0, i)
            filePattern = filePattern.substring(i+3)
        }

        def folder
        def pattern
        def matcher
        int p = filePattern.indexOf('*')
        if( p != -1 ) {
            i = filePattern.substring(0,p).lastIndexOf('/')
        }
        else if( (matcher=FileHelper.GLOB_CURLY_BRACKETS.matcher(filePattern)).matches() ) {
            def prefix = matcher.group(1)
            if( prefix ) {
                i = prefix.contains('/') ? prefix.lastIndexOf('/') : -1
            }
            else {
                i = matcher.start(2) -1
            }
        }
        else {
            i = filePattern.lastIndexOf('/')
        }

        if( i != -1 ) {
            folder = filePattern.substring(0,i+1)
            pattern = filePattern.substring(i+1)
        }
        else {
            folder = './'
            pattern = filePattern
        }

        return [ folder, pattern, scheme ]
    }


    static FileSystem fileSystemForScheme(String scheme) {
        ( !scheme
                ? FileSystems.getDefault()
                : getOrCreateFileSystemFor(scheme) )
    }


}

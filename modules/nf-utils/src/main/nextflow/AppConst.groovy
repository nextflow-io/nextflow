package nextflow

import nextflow.util.BasicFileHelper

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths

class AppConst {
    /**
     * The application main package name
     */
    static public final String MAIN_PACKAGE = System.getenv('NXF_MAIN_PACKAGE') ?: 'nextflow'

    /**
     * The application main name
     */
    static public final String APP_NAME = MAIN_PACKAGE

    /**
     * The application home folder
     */
    static public final Path APP_HOME_DIR = getHomeDir(APP_NAME)

    static public final File DEFAULT_ROOT = System.getenv('NXF_ASSETS') ? new File(System.getenv('NXF_ASSETS')) : APP_HOME_DIR.resolve('assets').toFile()

    static public final String MANIFEST_FILE_NAME = 'nextflow.config'

    static public final String DEFAULT_MAIN_FILE_NAME = 'main.nf'

    static Path sysHome() {
        def home = System.getProperty("user.home")
        if( !home || home=='?' )
            home = System.getenv('HOME')
        if( !home )
            throw new IllegalStateException("Unable to detect system home path - Make sure the variable HOME or NXF_HOME is defined in your environment")
        return Path.of(home)
    }

    private static Path getHomeDir(String appname) {
        final home = System.getenv('NXF_HOME')
        final result = home ? Paths.get(home) : sysHome().resolve(".$appname")

        if( !Files.exists(result) && !BasicFileHelper.mkdir(result) ) {
            throw new IllegalStateException("Cannot create path '${result}' -- check file system access permission")
        }

        return result
    }

}

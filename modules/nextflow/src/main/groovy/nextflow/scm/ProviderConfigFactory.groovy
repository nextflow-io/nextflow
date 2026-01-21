package nextflow.scm

import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Const
import nextflow.config.ConfigParserFactory
import nextflow.exception.AbortOperationException
import nextflow.exception.ConfigParseException
import nextflow.file.FileHelper
import nextflow.util.StringUtils

import java.nio.file.NoSuchFileException
import java.nio.file.Path

@Slf4j
class ProviderConfigFactory {

    @PackageScope
    static Path DEFAULT_SCM_FILE = Const.APP_HOME_DIR.resolve('scm')

    @PackageScope
    static Map<String,String> env = new HashMap<>(System.getenv())

    static Path getScmConfigPath() {
        def cfg = env.get('NXF_SCM_FILE')
        if( !cfg ) {
            log.debug "Using SCM config path: ${DEFAULT_SCM_FILE}"
            return DEFAULT_SCM_FILE
        }
        log.debug "Detected SCM custom path: $cfg"
        // check and return it if valid
        return FileHelper.asPath(cfg)
    }

    @PackageScope
    static Map parse(String text) {
        def slurper = ConfigParserFactory.create()
        slurper.setBinding(env)
        return slurper.parse(text)
    }

    @PackageScope
    static List<ProviderConfig> createFromText(String text) {
        def config = parse(text)
        def result = NextflowRepositoryFactoryLoader.instance.createFromMap(config)
        return result
    }

    @PackageScope
    static Map getFromFile(Path file) {
        try {
            // note: since this can be a remote file via e.g. read via HTTP
            // it should not read more than once
            final content = file.text
            dumpScmContent(file, content)
            final result = parse(content)
            dumpConfig(result)
            return result
        }
        catch ( NoSuchFileException | FileNotFoundException e) {
            if( file == DEFAULT_SCM_FILE ) {
                return new LinkedHashMap()
            }
            else
                throw new AbortOperationException("Missing SCM config file: ${file.toUriString()} - Check the env variable NXF_SCM_FILE")
        }
        catch ( UnknownHostException e ) {
            final message = "Unable to access config file '${file?.toUriString()}' -- Unknown host: ${e}"
            throw new ConfigParseException(message,e)
        }
        catch( IOException e ) {
            final message = "Unable to access config file '${file?.toUriString()}' -- Cause: ${e.message?:e.toString()}"
            throw new ConfigParseException(message,e)
        }
        catch( Exception e ) {
            final message = "Failed to parse config file '${file?.toUriString()}' -- Cause: ${e.message?:e.toString()}"
            throw new ConfigParseException(message,e)
        }
    }

    static private void dumpScmContent(Path file, String content) {
        try {
            log.trace "Parsing SCM config path: ${file.toUriString()}\n${StringUtils.stripSecrets(content)}\n"
        }catch(Exception e){
            log.debug "Error dumping configuration ${file.toUriString()}", e
        }
    }

    static private void dumpConfig(Map config) {
        try {
            log.debug "Detected SCM config: ${StringUtils.stripSecrets(config).toMapString()}"
        }
        catch (Throwable e) {
            log.debug "Failed to dump SCM config: ${e.message ?: e}", e
        }
    }

    static Map getDefault() {
        final file = getScmConfigPath()
        return getFromFile(file)
    }
}

package nextflow.container

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ShifterBuilder implements ContainerBuilder {

    static private String SHIFTER_HELPERS = '''
        function shifter_img() {
          local cmd=$1
          local image=$2
          shifterimg -v $cmd $image |  awk -F: '$0~/status/{gsub("[\\", ]","",$2);print $2}\'
        }

        function shifter_pull() {
          local image=$1
          local STATUS=$(shifter_img lookup $image)
          if [[ $STATUS != READY && $STATUS != '' ]]; then
           STATUS=$(shifter_img pull $image)
           while [[ $STATUS != READY && $STATUS != FAILURE && $STATUS != '' ]]; do
             sleep 5
             STATUS=$(shifter_img pull $image)
           done
          fi

          [[ $STATUS == FAILURE || $STATUS == '' ]] && echo "Shifter failed to pull image \\`$image\\`" >&2  && exit 1
        }
        '''.stripIndent()


    private String entryPoint

    private String image

    private boolean verbose

    private String runCommand

    ShifterBuilder( String image ) {
        assert image
        this.image = image
    }

    String getRunCommand() { runCommand }

    @Override
    String build(StringBuilder result) {
        assert image

        result << 'shifter '

        if( verbose )
            result << '--verbose '

        result << '--image ' << image

        if( entryPoint )
            result << ' ' << entryPoint

        runCommand = result.toString()
    }

    ShifterBuilder params( Map params ) {

        if( params.containsKey('verbose') )
            this.verbose = params.verbose.toString() == 'true'

        if( params.containsKey('entry') )
            this.entryPoint = params.entry

        return this
    }

    StringBuilder appendHelpers( StringBuilder wrapper ) {
        wrapper << SHIFTER_HELPERS << '\n'
    }

    @Override
    StringBuilder appendRunCommand( StringBuilder wrapper ) {
        wrapper << 'shifter_pull ' << image << '\n'
        wrapper << runCommand
    }

    /**
     * Normalize Shifter image name adding `docker:` prefix or `:latest`
     * when required
     *
     * @param imageName The container image name
     * @param shifterConfig Shifter configuration map
     * @return Image name in Shifter canonical format
     */
    static String normalizeImageName( String imageName, Map shifterConfig ) {

        if( !imageName )
            return null

        def items = imageName.tokenize(':')
        if( items.size()==3 ) {
            // it is in the canonical form i.e. `type:image:tag`
            return imageName
        }

        if( items.size()==1 ) {
            return "docker:$imageName:latest"
        }

        return !imageName.startsWith("docker:") ? "docker:$imageName" : "$imageName:latest"
    }
}

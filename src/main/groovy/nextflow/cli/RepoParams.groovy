package nextflow.cli

import com.beust.jcommander.Parameter

/**
 * Created by mchatzou on 8/19/14.
 */
trait RepoParams {

    @Parameter(names=['-rep', '-repository'], description = 'Repository code provider to clone - It can be either github or bitbucket. If it is a private repository please consider using options -repID and -repPWD')
    String repository="github"

    @Parameter(names='-user', description = 'Private repository username')
    String repID

    @Parameter(names='-pwd', description = 'Private repository password')
    String repPWD


    String getPassword() {

        if( !repID ) {
            return null
        }

        if( repPWD ) {
            return repPWD
        }

        print "Enter your $repository password: "
        char[] pwd = System.console().readPassword()
        new String(pwd)
    }

}
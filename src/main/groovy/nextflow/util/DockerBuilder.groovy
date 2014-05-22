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
import java.nio.file.Path

import groovy.transform.PackageScope
import nextflow.processor.FileHolder

/**
 * Helper methods to handle Docker containers
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class DockerBuilder {

    final String containerName

    final List env = []

    List<Path> mounts = []

    DockerBuilder( String name ) {
        this.containerName = name
    }

    DockerBuilder addEnv( entry ) {
        env.add(entry)
        return this
    }

    DockerBuilder addMount( Path path ) {
        if( path )
            mounts.add(path)
        return this
    }

    DockerBuilder addMountForInputs( Map<?,List<FileHolder>> inputFiles ) {

        mounts.addAll( inputFilesToPaths(inputFiles) )
        return this
    }


    CharSequence build(StringBuilder result = new StringBuilder()) {
        assert containerName

        result << 'docker run --rm -u $(id -u) '

        // add the environment
        for( def entry : env ) {
            result << makeEnv(entry) << ' '
        }

        // mount the input folders
        result << makeVolumes(mounts)
        result << ' -w $PWD '

        // finally the container name
        result << (containerName)

        return result

    }


    /**
     * Get the volumes command line options for the given list of input files
     *
     * @param mountPaths
     * @param binDir
     * @param result
     * @return
     */
    @PackageScope
    static CharSequence makeVolumes(List<Path> mountPaths, StringBuilder result = new StringBuilder() ) {

        // always provide external /tmp
        result << '-v ${NXF_SCRATCH:-$(mktemp -d)}:/tmp '

        // find the longest commons paths and mount only them
        def trie = new PathTrie()
        mountPaths.each { trie.add(it) }

        def paths = trie.longest()
        paths.each{ if(it) result << "-v $it:$it " }

        // -- append by default the current path
        result << '-v $PWD:$PWD'

        return result
    }

    /**
     * Get the env command line option for the give environment definition
     *
     * @param env
     * @param result
     * @return
     */
    @PackageScope
    static CharSequence makeEnv( env, StringBuilder result = new StringBuilder() ) {
        // append the environment configuration
        if( env instanceof File ) {
            env = env.toPath()
        }
        if( env instanceof Path ) {
            result << '-e "BASH_ENV=' << env.getName() << '"'
        }
        else if( env instanceof Map ) {
            env.eachWithIndex { key, value, index ->
                if( index ) result << ' '
                result << ("-e \"$key=$value\"")
            }
        }
        else if( env instanceof String && env.contains('=') ) {
            result << '-e "' << env << '"'
        }
        else if( env ) {
            throw new IllegalArgumentException("Not a valid Docker environment value: $env [${env.class.name}]")
        }

        return result
    }


    @PackageScope
    static List<Path> inputFilesToPaths( Map<?,List<FileHolder>> inputFiles ) {

        def List<Path> files = []
        inputFiles.each { param, holders ->

            holders.each { FileHolder item ->
                def path = item.storePath?.getParent()
                if( path ) files << path
            }

        }
        return files

    }

}

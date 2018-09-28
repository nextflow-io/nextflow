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

package nextflow.scm

import java.nio.file.Files
import java.nio.file.Path

import org.eclipse.jgit.api.Git
import spock.lang.Shared
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class UpdateModuleTest extends Specification {

    @Shared
    Path baseFolder

    def setupSpec() {
        baseFolder = Files.createTempDirectory('test').toAbsolutePath()

        // create module a
        baseFolder.resolve('prj_aaa').mkdir()
        def dir = baseFolder.resolve('prj_aaa').toFile()
        def init = Git.init()
        def repo = init.setDirectory( dir ).call()
        new File(dir, 'file1.txt').text = 'Hello'
        new File(dir, 'file2.log').text = 'World'
        repo.add().addFilepattern('.').call()
        repo.commit().setAll(true).setMessage('First commit').call()
        repo.close()

        // create module b
        baseFolder.resolve('prj_bbb').mkdir()
        dir = baseFolder.resolve('prj_bbb').toFile()
        init = Git.init()
        repo = init.setDirectory( dir ).call()
        new File(dir, 'file1.txt').text = 'Ciao'
        new File(dir, 'file2.log').text = 'Mondo'
        repo.add().addFilepattern('.').call()
        repo.commit().setAll(true).setMessage('First commit').call()
        repo.close()

        // create module c
        baseFolder.resolve('prj_ccc').mkdir()
        dir = baseFolder.resolve('prj_ccc').toFile()
        init = Git.init()
        repo = init.setDirectory( dir ).call()
        new File(dir, 'file-x.txt').text = 'x'
        new File(dir, 'file-y.txt').text = 'y'
        repo.add().addFilepattern('.').call()
        repo.commit().setAll(true).setMessage('First commit').call()
        repo.close()

    }

    def cleanupSpec() {
        baseFolder?.deleteDir()
    }

    def 'should clone and update submodules' () {

        setup:
        // create the main project
        baseFolder.resolve('pipe_x').mkdirs()
        def dir = baseFolder.resolve('pipe_x').toFile()
        def init = Git.init()
        def repo = init.setDirectory( dir ).call()
        // create two git submodules in this repository
        repo.submoduleAdd().setPath('prj_aaa').setURI( baseFolder.resolve('prj_aaa/.git').toString() ).call()
        repo.submoduleAdd().setPath('prj_bbb').setURI( baseFolder.resolve('prj_bbb/.git').toString() ).call()
        // add the main file and commit all
        new File(dir, 'main.nf').text = 'main script'
        repo.add().addFilepattern('.').call()
        repo.commit().setAll(true).setMessage('First commit').call()

        // set the root folder for the asset manager
        def target = baseFolder.resolve('target')
        AssetManager.root = target.toFile()

        when:
        def manager = new AssetManager("file:${baseFolder}/pipe_x")
        manager.download()
        manager.updateModules()

        then:
        target.resolve('local/pipe_x').exists()
        target.resolve('local/pipe_x/.git').exists()
        target.resolve('local/pipe_x/main.nf').exists()

        target.resolve('local/pipe_x/prj_aaa').exists()
        target.resolve('local/pipe_x/prj_aaa/file1.txt').text == 'Hello'
        target.resolve('local/pipe_x/prj_aaa/file2.log').text == 'World'

        target.resolve('local/pipe_x/prj_bbb').exists()
        target.resolve('local/pipe_x/prj_bbb/file1.txt').text == 'Ciao'
        target.resolve('local/pipe_x/prj_bbb/file2.log').text == 'Mondo'
    }


    def 'should not clone submodules' () {

        setup:
        // create the main project
        baseFolder.resolve('pipe_2').mkdirs()
        def dir = baseFolder.resolve('pipe_2').toFile()
        def init = Git.init()
        def repo = init.setDirectory( dir ).call()
        repo.submoduleAdd().setPath('prj_aaa').setURI( baseFolder.resolve('prj_aaa/.git').toString() ).call()
        repo.submoduleAdd().setPath('prj_bbb').setURI( baseFolder.resolve('prj_bbb/.git').toString() ).call()
        new File(dir, 'main.nf').text = 'main script'
        new File(dir, 'nextflow.config').text = 'manifest.gitmodules = false'  // note: <-- this setting switch-off the submodule update
        repo.add().addFilepattern('.').call()
        repo.commit().setAll(true).setMessage('First commit').call()

        def target = baseFolder.resolve('target2')
        AssetManager.root = target.toFile()

        when:
        def manager = new AssetManager( "file:${baseFolder}/pipe_2" )
        manager.download()
        manager.updateModules()

        then:
        target.resolve('local/pipe_2').exists()
        target.resolve('local/pipe_2/.git').exists()
        target.resolve('local/pipe_2/main.nf').exists()

        !target.resolve('local/pipe_2/prj_aaa').exists()
        !target.resolve('local/pipe_2/prj_bbb').exists()
    }

    def 'should selected submodules' () {

        setup:
        // create the main project
        baseFolder.resolve('pipe_3').mkdirs()
        def dir = baseFolder.resolve('pipe_3').toFile()
        def init = Git.init()
        def repo = init.setDirectory( dir ).call()
        repo.submoduleAdd().setPath('prj_aaa').setURI( baseFolder.resolve('prj_aaa/.git').toString() ).call()
        repo.submoduleAdd().setPath('prj_bbb').setURI( baseFolder.resolve('prj_bbb/.git').toString() ).call()
        repo.submoduleAdd().setPath('prj_ccc').setURI( baseFolder.resolve('prj_ccc/.git').toString() ).call()
        new File(dir, 'main.nf').text = 'main script'
        new File(dir, 'nextflow.config').text = "manifest.gitmodules = 'prj_bbb,prj_ccc'"  // note: <-- this setting switch-off the submodule update
        repo.add().addFilepattern('.').call()
        repo.commit().setAll(true).setMessage('First commit').call()

        def target = baseFolder.resolve('target3')
        AssetManager.root = target.toFile()

        when:
        def manager = new AssetManager( "file:${baseFolder}/pipe_3" )
        manager.download()
        manager.updateModules()

        then:
        target.resolve('local/pipe_3').exists()
        target.resolve('local/pipe_3/.git').exists()
        target.resolve('local/pipe_3/main.nf').exists()

        !target.resolve('local/pipe_3/prj_aaa').exists()

        target.resolve('local/pipe_3/prj_bbb').exists()
        target.resolve('local/pipe_3/prj_bbb/file1.txt').text == 'Ciao'

        target.resolve('local/pipe_3/prj_ccc').exists()
        target.resolve('local/pipe_3/prj_ccc/file-x.txt').text == 'x'

    }


}

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

package nextflow.scm
import java.nio.file.Files

import nextflow.exception.AbortOperationException
import spock.lang.Requires
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AssetManagerTest extends Specification {

    def testList() {

        given:
        def folder = Files.createTempDirectory('test')
        folder.resolve('cbcrg/pipe1').mkdirs()
        folder.resolve('cbcrg/pipe2').mkdirs()
        folder.resolve('ncbi/blast').mkdirs()

        def manager = new AssetManager().setRoot(folder.toFile())

        when:
        def list = manager.list()
        then:
        list.sort() == ['cbcrg/pipe1','cbcrg/pipe2','ncbi/blast']

        expect:
        manager.find('blast') == 'ncbi/blast'
        manager.find('pipe1') == 'cbcrg/pipe1'
        manager.find('pipe') as Set == ['cbcrg/pipe1', 'cbcrg/pipe2'] as Set

        cleanup:
        folder?.deleteDir()

    }


    def testResolveName() {

        given:
        def folder = Files.createTempDirectory('test')
        folder.resolve('cbcrg/pipe1').mkdirs()
        folder.resolve('cbcrg/pipe2').mkdirs()
        folder.resolve('ncbi/blast').mkdirs()

        def manager = new AssetManager().setRoot(folder.toFile())

        when:
        def result = manager.resolveName('x/y')
        then:
        result == 'x/y'

        when:
        result = manager.resolveName('blast')
        then:
        result == 'ncbi/blast'

        when:
        result = manager.resolveName('ncbi/blast/script.nf')
        then:
        result == 'ncbi/blast'

        when:
        result = manager.resolveName('blast/script.nf')
        then:
        result == 'ncbi/blast'

        when:
        manager.resolveName('pipe')
        then:
        thrown(AbortOperationException)

        when:
        manager.resolveName('pipe/alpha/beta')
        then:
        thrown(AbortOperationException)

        cleanup:
        folder?.deleteDir()

    }


    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def testPull() {

        given:
        def (user,pwd) = System.getenv('NXF_GITHUB_ACCESS_TOKEN')?.tokenize(':')
        def folder = Files.createTempDirectory('test')
        def manager = new AssetManager(user:user, pwd: pwd).setRoot(folder.toFile()).setPipeline('nextflow-io/hello')

        when:
        manager.download()
        then:
        folder.resolve('nextflow-io/hello/.git').isDirectory()

        when:
        manager.download()
        then:
        noExceptionThrown()

        cleanup:
        folder?.deleteDir()

    }


    @Requires({System.getenv('NXF_GITHUB_ACCESS_TOKEN')})
    def testClone() {

        given:
        def (user,pwd) = System.getenv('NXF_GITHUB_ACCESS_TOKEN')?.tokenize(':')
        def dir = Files.createTempDirectory('test')
        def manager = new AssetManager(user:user, pwd: pwd, pipeline: 'pditommaso/awesome-pipeline')

        when:
        manager.clone(dir.toFile())

        then:
        dir.resolve('README.md').exists()
        dir.resolve('.git').isDirectory()

        cleanup:
        dir?.deleteDir()

    }

    def testGetScriptName() {

        given:
        def dir = Files.createTempDirectory('test')
        dir.resolve('sub1').mkdir()
        dir.resolve('sub1/nextflow.config').text = "manifest.mainScript = 'pippo.nf'"
        dir.resolve('sub2').mkdir()


        when:
        def holder = new AssetManager()
        holder.localPath = dir.resolve('sub1').toFile()
        then:
        holder.getMainScriptName() == 'pippo.nf'

        when:
        holder = new AssetManager()
        holder.localPath = dir.resolve('sub2').toFile()
        then:
        holder.getMainScriptName() == 'main.nf'

        when:
        holder = new AssetManager()
        holder.localPath = dir.resolve('sub1').toFile()
        holder.resolveName('nextflow/hello')
        then:
        holder.getMainScriptName() == 'pippo.nf'

        when:
        holder = new AssetManager()
        holder.localPath = dir.resolve('sub1').toFile()
        holder.resolveName('nextflow/hello/my-script.nf')
        then:
        holder.getMainScriptName() == 'my-script.nf'

        when:
        holder = new AssetManager()
        holder.localPath = dir.resolve('sub1').toFile()
        then:
        holder.resolveName('nextflow-io/hello/x/y/z/my-script.nf') == 'nextflow-io/hello'
        holder.getMainScriptName() == 'x/y/z/my-script.nf'

        when:
        holder = new AssetManager()
        holder.localPath = dir.resolve('sub1').toFile()
        then:
        holder.resolveName('nextflow-io/hello/my-script.nf') == 'nextflow-io/hello'
        holder.getMainScriptName() == 'my-script.nf'

        when:
        holder = new AssetManager()
        holder.localPath = dir.resolve('sub1').toFile()
        then:
        holder.resolveName('hello/my-script.nf') == 'nextflow-io/hello'
        holder.getMainScriptName() == 'my-script.nf'

        cleanup:
        dir?.deleteDir()

    }

    def testCreateProviderFor(){

        when:
        def manager= [ pipeline:'x/y', user:'maria', pwd: 'whatever' ] as AssetManager
        def repo=manager.createHubProviderFor('github')
        then:
        repo instanceof GithubRepositoryProvider
        repo.user=="maria"
        repo.pwd=="whatever"
        repo.pipeline=="x/y"

        when:
        manager= [ pipeline:'x/y', user:'maria', pwd: 'whatever' ] as AssetManager
        repo=manager.createHubProviderFor('bitbucket')
        then:
        repo instanceof BitbucketRepositoryProvider
        repo.user=="maria"
        repo.pwd=="whatever"
        repo.pipeline=="x/y"

        when:
        manager = [ pipeline: 'project/foo', user: 'paolo', pwd: '12345'] as AssetManager
        repo = manager.createHubProviderFor('gitlab')
        then:
        repo instanceof GitlabRepositoryProvider
        repo.user == 'paolo'
        repo.pwd == '12345'
        repo.pipeline == 'project/foo'

        when:
        manager = [ ] as AssetManager
        manager.createHubProviderFor('xxx')
        then:
        thrown(AbortOperationException)


    }


    def testReadManifest () {

        given:
        def config =
                '''
                manifest {
                    homePage = 'http://foo.com'
                    mainScript = 'hello.nf'
                    defaultBranch = 'super-stuff'
                    description = 'This pipeline do this and that'
                }
                '''
        def dir = Files.createTempDirectory('test')
        dir.resolve('foo/bar').mkdirs()
        dir.resolve('foo/bar/nextflow.config').text = config

        when:
        def holder = new AssetManager()
        holder.root = dir.toFile()
        holder.setPipeline('foo/bar')
        then:
        holder.getMainScriptName() == 'hello.nf'
        holder.getDefaultBranch() == 'super-stuff'
        holder.getHomePage() == 'http://foo.com'
        holder.getDescription() == 'This pipeline do this and that'


        cleanup:
        dir.deleteDir()

    }

    def testReadManifest2 () {

        given:
        def dir = Files.createTempDirectory('test')
        dir.resolve('foo/bar').mkdirs()
        dir.resolve('foo/bar/nextflow.config').text = 'empty: 1'

        when:
        def holder = new AssetManager()
        holder.root = dir.toFile()
        holder.setPipeline('foo/bar')

        then:
        holder.getMainScriptName() == 'main.nf'
        holder.getDefaultBranch() == 'master'
        holder.getHomePage() == 'https://github.com/foo/bar'
        holder.getDescription() == null

        cleanup:
        dir.deleteDir()

    }
}

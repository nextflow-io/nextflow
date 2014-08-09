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

package nextflow.script
import java.nio.file.Files

import nextflow.exception.AbortOperationException
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
        manager.find('pipe') == ['cbcrg/pipe1', 'cbcrg/pipe2']

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


    def testPull() {

        given:
        def folder = Files.createTempDirectory('test')
        def manager = new AssetManager().setRoot(folder.toFile()).setPipeline('nextflow-io/hello')

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

    def testClone() {

        given:
        def dir = Files.createTempDirectory('test')
        def manager = new AssetManager('pditommaso/awesome-pipeline')

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
        dir.resolve('sub1/.PIPELINE-INF').text = 'main-script = pippo.nf'
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

    def testReadManifestFrom() {

        given:
        def dir = Files.createTempDirectory('test').toFile()
        new File(dir,'response').text = '''{
              "name": ".PIPELINE-INF",
              "path": ".PIPELINE-INF",
              "sha": "318cb862eda4a1782cd1e1333d913a8e3c61f0a1",
              "size": 21,
              "url": "https://api.github.com/repos/cbcrg/piper-nf/contents/.PIPELINE-INF?ref=master",
              "html_url": "https://github.com/cbcrg/piper-nf/blob/master/.PIPELINE-INF",
              "git_url": "https://api.github.com/repos/cbcrg/piper-nf/git/blobs/318cb862eda4a1782cd1e1333d913a8e3c61f0a1",
              "type": "file",
              "content": "bWFpbi1zY3JpcHQ6IHBpcGVyLm5m\\n",
              "encoding": "base64",
              "_links": {
                "self": "https://api.github.com/repos/cbcrg/piper-nf/contents/.PIPELINE-INF?ref=master",
                "git": "https://api.github.com/repos/cbcrg/piper-nf/git/blobs/318cb862eda4a1782cd1e1333d913a8e3c61f0a1",
                "html": "https://github.com/cbcrg/piper-nf/blob/master/.PIPELINE-INF"
              }
            }'''

        when:
        def manifest = AssetManager.readManifestFrom("file://$dir/response")

        then:
        manifest.size()==1
        manifest.get('main-script') == 'piper.nf'

        cleanup:
        dir?.deleteDir()

    }

}

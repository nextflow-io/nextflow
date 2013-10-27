package nextflow

import java.nio.file.Files

import spock.lang.Specification


/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ChannelTest extends Specification {

    def testGetPathAndPattern () {

        expect:
        Channel.getFolderAndPattern( '/some/file/name.txt' ) == ['/some/file/', 'name.txt']
        Channel.getFolderAndPattern( '/some/file/*.txt' ) == ['/some/file/', '*.txt']
        Channel.getFolderAndPattern( '/some/file/*' ) == ['/some/file/', '*']
        Channel.getFolderAndPattern( '/some/file/' ) == ['/some/file/', '']

        Channel.getFolderAndPattern( '/some/file/**/*.txt' ) == ['/some/file/', '**/*.txt']

        Channel.getFolderAndPattern( 'dxfs:///some/file/**/*.txt' ) == ['dxfs:///some/file/', '**/*.txt']
        Channel.getFolderAndPattern( 'dxfs://some/file/**/*.txt' ) == ['dxfs://some/file/', '**/*.txt']
        Channel.getFolderAndPattern( 'dxfs://*.txt' ) == ['dxfs://.', '*.txt']
        Channel.getFolderAndPattern( 'dxfs:///*.txt' ) == ['dxfs:///', '*.txt']
        Channel.getFolderAndPattern( 'dxfs:///**/*.txt' ) == ['dxfs:///', '**/*.txt']
    }


    def testFiles() {

        setup:
        def folder = Files.createTempDirectory('testFiles')
        def file1 = Files.createFile(folder.resolve('file1.txt'))
        def file2 = Files.createFile(folder.resolve('file2.txt'))
        def file3 = Files.createFile(folder.resolve('file3.txt'))
        def file4 = Files.createFile(folder.resolve('file4.log'))
        def sub1 = Files.createDirectories(folder.resolve('sub1'))
        def file5 = Files.createFile(sub1.resolve('file5.log'))
        def file6 = Files.createFile(sub1.resolve('file6.txt'))

        when:
        def channel = Channel.files("$folder/*.txt")
        then:
        channel.val == file1
        channel.val == file2
        channel.val == file3
        channel.val == Channel.STOP


        when:
        def channel2 = Channel.files("$folder/**.txt")
        then:
        channel2.val.toString() == file1.toString()
        channel2.val.toString() == file2.toString()
        channel2.val.toString() == file3.toString()
        channel2.val.toString() == file6.toString()
        channel2.val == Channel.STOP

        when:
        def channel3 = Channel.files("$folder/sub1/**.log")
        then:
        channel3.val.toString() == file5.toString()
        channel3.val == Channel.STOP

        cleanup:
        folder.deleteDir()



    }


}
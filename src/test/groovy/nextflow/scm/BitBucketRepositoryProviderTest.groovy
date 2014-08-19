package nextflow.scm

import spock.lang.Specification

/**
 * Created by mchatzou on 8/14/14.
 */
class BitBucketRepositoryProviderTest extends Specification {

    def testBitbucketCloneURL() {
        when:
        def url = new BitbucketRepositoryProvider(pipeline: 'pditommaso/tutorial').getCloneURL()
        then:
        url == 'https://bitbucket.org/pditommaso/tutorial.git'
    }


    def testGetHomePage() {
        expect:
        new BitbucketRepositoryProvider(pipeline: 'mariach/tutorial').getHomePage() == "https://bitbucket.org/mariach/tutorial"
    }

    def testGitCloneUrlPrivate() {

        when:
        def url = new BitbucketRepositoryProvider(pipeline: 'mariach/tutorial', user: "mariach",password: "t_coffee1").getCloneURL()
        then:
        url == 'https://mariach@bitbucket.org/mariach/tutorial.git'

    }


    def testReadContent() {

        when:
        def repo = new BitbucketRepositoryProvider(pipeline: 'pditommaso/tutorial')
        def result = repo.readContent('main.nf')
        then:
        result.trim().startsWith('#!/usr/bin/env nextflow')

    }
}

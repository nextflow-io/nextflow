package nextflow.scm

import spock.lang.Specification

/**
 * Created by mchatzou on 8/14/14.
 */
class GithubRepositoryProviderTest extends Specification {

    def testGitCloneUrl() {

        when:
        def url = new GithubRepositoryProvider(pipeline: 'nextflow-io/hello').getCloneURL()
        then:
        url == 'https://github.com/nextflow-io/hello.git'

    }

    def testGetHomePage() {
        expect:
        new GithubRepositoryProvider(pipeline: 'nextflow-io/hello').getHomePage() == "https://github.com/nextflow-io/hello"
    }


    def testReadContent() {

        when:
        def repo = new GithubRepositoryProvider(pipeline: 'nextflow-io/hello')
        def result = repo.readContent('main.nf')
        then:
        result.trim().startsWith('#!/usr/bin/env nextflow')

    }
}


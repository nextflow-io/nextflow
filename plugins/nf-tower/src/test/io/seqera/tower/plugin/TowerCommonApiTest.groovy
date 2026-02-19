package io.seqera.tower.plugin

import spock.lang.Specification

class TowerCommonApiTest extends Specification{


    def 'should build URL without query params'() {
        given:
        def api = new TowerCommonApi()

        when:
        def url = api.buildUrl('https://api.cloud.seqera.io', '/workflow/launch', [:])

        then:
        url == 'https://api.cloud.seqera.io/workflow/launch'
    }

    def 'should build URL with query params'() {
        given:
        def api = new TowerCommonApi()

        when:
        def url = api.buildUrl('https://api.cloud.seqera.io', '/workflow/launch', [workspaceId: '12345'])

        then:
        url.contains('https://api.cloud.seqera.io/workflow/launch?')
        url.contains('workspaceId=12345')
    }

    def 'should URL encode query params'() {
        given:
        def api = new TowerCommonApi()

        when:
        def url = api.buildUrl('https://api.cloud.seqera.io', '/workflow', [name: 'test workflow'])

        then:
        url.contains('name=test+workflow')
    }
}

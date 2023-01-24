/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
 */

package nextflow.executor

import java.nio.file.Path

import nextflow.Global
import nextflow.SysEnv
import nextflow.cloud.aws.util.S3PathFactory
import nextflow.fusion.FusionScriptLauncher
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FusionScriptLauncherS3Test extends Specification {

    def 'should get container mount' () {
        given:
        def fusion = new FusionScriptLauncher(scheme: 's3')

        when:
        def result = fusion.toContainerMount(S3PathFactory.parse('s3://foo/a/b/c.txt'))
        then:
        result == Path.of('/fusion/s3/foo/a/b/c.txt')

        when:
        result = fusion.toContainerMount(S3PathFactory.parse('s3://foo/a/x/y.txt'))
        then:
        result == Path.of('/fusion/s3/foo/a/x/y.txt')

        when:
        result = fusion.toContainerMount(S3PathFactory.parse('s3://bar/z.txt'))
        then:
        result == Path.of('/fusion/s3/bar/z.txt')


        expect:
        fusion.fusionBuckets() == [ 'foo', 'bar' ] as Set

    }


    def 'should get fusion env with s3 endpoint' () {
        given:
        Global.config = [:]
        and:
        SysEnv.push([AWS_S3_ENDPOINT: 'http://foo.com'])
        and:
        def fusion = new FusionScriptLauncher(
                scheme: 's3',
                buckets: ['foo'] as Set,
                remoteWorkDir: S3PathFactory.parse('s3://foo/work'))

        expect:
        fusion.fusionEnv() == [AWS_S3_ENDPOINT: 'http://foo.com',
                               NXF_FUSION_WORK: '/fusion/s3/foo/work']

        cleanup:
        SysEnv.pop()
    }

    def 'should get fusion env with aws credentials' () {
        given:
        SysEnv.push([AWS_ACCESS_KEY_ID: 'xxx', AWS_SECRET_ACCESS_KEY: 'zzz'])
        and:
        Global.config = [fusion: [exportAwsAccessKeys: true]]
        and:
        def fusion = new FusionScriptLauncher(
                scheme: 's3',
                buckets: ['foo'] as Set,
                remoteWorkDir: S3PathFactory.parse('s3://foo/work'))

        expect:
        fusion.fusionEnv() == [AWS_ACCESS_KEY_ID: 'xxx',
                               AWS_SECRET_ACCESS_KEY: 'zzz',
                               NXF_FUSION_WORK: '/fusion/s3/foo/work']

        cleanup:
        Global.config = null
        SysEnv.pop()
    }

    def 'should get fusion env with aws credentials in nextflow config' () {
        given:
        SysEnv.push([:])
        and:
        def CONFIG = [fusion: [exportAwsAccessKeys: true], aws: [accessKey: 'k1', secretKey: 's1', client: [endpoint: 'http://minio.com']]]
        Global.config = CONFIG
        and:
        def fusion = new FusionScriptLauncher(
                scheme: 's3',
                buckets: ['foo'] as Set,
                remoteWorkDir: S3PathFactory.parse('s3://foo/work'))

        expect:
        fusion.fusionEnv() == [AWS_ACCESS_KEY_ID: 'k1',
                               AWS_SECRET_ACCESS_KEY: 's1',
                               AWS_S3_ENDPOINT: 'http://minio.com',
                               NXF_FUSION_WORK: '/fusion/s3/foo/work']

        cleanup:
        Global.config = null
        SysEnv.pop()
    }

}

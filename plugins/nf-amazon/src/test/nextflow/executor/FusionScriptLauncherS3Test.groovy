/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
 */

package nextflow.executor

import java.nio.file.Path

import nextflow.cloud.aws.util.S3PathFactory
import nextflow.executor.fusion.FusionScriptLauncher
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

}

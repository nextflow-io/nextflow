package nextflow

import static nextflow.extension.Bolts.*

import java.text.SimpleDateFormat

import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class NextflowMetaTest extends Specification {

    private String dateToString(Date date) {
        def tz = TimeZone.getTimeZone('UTC')
        def fmt = new SimpleDateFormat(DATETIME_FORMAT)
        fmt.setTimeZone(tz)
        fmt.format(date) + ' ' + tz.getDisplayName(true, TimeZone.SHORT)
    }

    def 'should convert to map'() {
        given:
        def meta = new NextflowMeta('10.12.0', 123, BuildInfo.timestampUTC)
        meta.enableDsl2()

        when:
        def map = meta.toJsonMap()
        then:
        map.version == '10.12.0'
        map.build == 123
        map.enable.dsl == 2
        dateToString((Date) map.timestamp) == BuildInfo.timestampUTC

    }

    @Unroll
    def 'should find dsl2 declaration: #SCRIPT'() {
        when:
        def result = NextflowMeta.checkDslMode(SCRIPT)
        then:
        result == RESULT

        where:
        RESULT | SCRIPT
        null   | 'hello'
        '1'    | 'nextflow.enable.dsl=1'
        and:
        '2'    | 'nextflow.enable.dsl=2'
        '2'    | 'nextflow.enable.dsl = 2'
        '2'    | 'nextflow.enable.dsl =  2;'
        '2'    | '#!/bin/env nextflow\nnextflow.enable.dsl=2\nprintln foo'
        '2'    | 'nextflow.enable.dsl =  2 // a comment'
        '2'    | 'nextflow.enable.dsl = 2;// another comment'
        '2'    | 'nextflow.enable.dsl=2 // another comment'
        and:
        '1'    | 'nextflow.enable.dsl = 1'
    }

    def 'should probe dsl1' () {

        expect:
        NextflowMeta.probeDsl1('''
            process foo {
                input:
                file transcriptome from transcriptome_file
            }
            ''')
        and:
        NextflowMeta.probeDsl1('''
            process foo {
                input:
                file transcriptome 
                    from transcriptome_file
            }
            ''')

        NextflowMeta.probeDsl1('''
            process foo {
                input: file transcriptome from transcriptome_file
            }
            ''')


        and:
        NextflowMeta.probeDsl1('''
            process foo {
                output:
                file transcriptome into transcriptome_file
            }
            ''')
        and:
        NextflowMeta.probeDsl1('''
            process foo {
                output:
                file transcriptome 
                    into transcriptome_file
            }
            ''')

        NextflowMeta.probeDsl1('''
            process foo {
                output: file transcriptome into transcriptome_file
            }
            ''')

        and:
        NextflowMeta.probeDsl1('''
            process foo {
                input: 
                file transcriptome from transcriptome_file
                output:
                file transcriptome into transcriptome_file
            }
            ''')

        and:
        !NextflowMeta.probeDsl1('''
            process foo {
                input: 
                file transcriptome 
                output:
                file transcriptome 
            }
            ''')

        and:
        !NextflowMeta.probeDsl1('''
            process foo {
                input: 
                file transcriptome from transcriptome_file
                output:
                file transcriptome into transcriptome_file
            }
            
            workflow {
                foo()
            }
            ''')

        and:
        !NextflowMeta.probeDsl1('''
            process foo {
                input: 
                file transcriptome from transcriptome_file
                output:
                file transcriptome into transcriptome_file
            }
            
            workflow{ foo() }
            ''')

        and:
        !NextflowMeta.probeDsl1('''
            process foo {
                input: 
                file transcriptome from transcriptome_file
                output:
                file transcriptome into transcriptome_file
            }
            
            workflow bar { 
                foo() 
            }
            ''')
    }

    def 'should check has workfloW'() {
        expect:
        NextflowMeta.hasWorkflowDef('''
            workflow bar {
                foo()
            }
            ''')

        and:
        NextflowMeta.hasWorkflowDef('''
            workflow bar { foo }
            ''')

        and:
        NextflowMeta.hasWorkflowDef('''
            workflow { foo }
            ''')

        and:
        NextflowMeta.hasWorkflowDef('''
            workflow{ foo }
            ''')

        and:
        NextflowMeta.hasWorkflowDef('''
            workflow{
                foo
            }
            ''')

        and:
        !NextflowMeta.hasWorkflowDef('''
            if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
                custom_runName = workflow.runName
            }

            ''')

        and:
        !NextflowMeta.hasWorkflowDef('''
            workflow.onComplete {

                // Set up the e-mail variables
                def subject = "[nf-core/chipseq] Successful: $workflow.runName"
                if (!workflow.success) {
                    subject = "[nf-core/chipseq] FAILED: $workflow.runName"
                }
                def email_fields = [:]
                email_fields['version'] = workflow.manifest.version

            ''')

        and:
        !NextflowMeta.hasWorkflowDef('''
            workflow . onComplete {
                def subject = "[nf-core/chipseq] Successful: $workflow.runName"
            ''')
    }
}

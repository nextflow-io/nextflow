/*
 * Fusion snapshots e2e: a counter task interrupts its own spot instance
 * halfway through via AWS FIS. If Fusion restores the snapshot, the counter
 * continues on a new instance instead of starting again from 1.
 */

process COUNT {
    publishDir params.outdir, mode: 'copy'
    conda 'conda-forge::awscli conda-forge::curl'
    cpus 1
    memory '2 GB'

    output:
    path 'count.txt'

    script:
    """
    # 'local x=\$(cmd)' hides a failing cmd from 'set -e', so declare and assign separately
    imds() {
        local t
        t=\$(curl -sf -X PUT http://169.254.169.254/latest/api/token -H 'X-aws-ec2-metadata-token-ttl-seconds: 60')
        curl -sf -H "X-aws-ec2-metadata-token: \$t" http://169.254.169.254/latest/meta-data/\$1
    }

    interrupt_self() {
        local iid eid state
        iid=\$(imds instance-id)
        [ "\$(imds instance-life-cycle)" = spot ] || { echo "instance \$iid is not spot" >&2; exit 1; }
        echo "tagging \$iid and starting FIS template ${params.fis_template}"
        aws ec2 create-tags --resources \$iid --tags Key=fis-target,Value=nf-snapshot-e2e
        # a tag left behind makes the next experiment interrupt this instance too
        trap "aws ec2 delete-tags --resources \$iid --tags Key=fis-target || true" EXIT
        eid=\$(aws fis start-experiment --experiment-template-id ${params.fis_template} --query experiment.id --output text)
        echo "FIS experiment \$eid started"
        # untag once FIS has resolved its targets, so a reused instance is never hit by a later run
        # fis_wait <tries> <states>: returns 0 once the experiment reaches one of <states>, exits on failure
        fis_wait() {
            for _ in \$(seq \$1); do
                state=\$(aws fis get-experiment --id \$eid --query experiment.state.status --output text)
                case \$state in
                    failed|stopped)
                        aws fis get-experiment --id \$eid --query experiment.state.reason --output text >&2
                        exit 1 ;;
                esac
                [[ \$state =~ ^(\$2)\$ ]] && return 0
                sleep 2
            done
            return 1
        }
        fis_wait 60 'running|completed' || { echo "FIS experiment \$eid still '\$state' after 2 min" >&2; exit 1; }
        aws ec2 delete-tags --resources \$iid --tags Key=fis-target
        trap - EXIT
        echo "FIS experiment \$eid \$state"
        # catch a failing action while the instance is still alive, well inside the 2 min warning;
        # still 'running' after 30 s is fine, the action may last until the interruption
        fis_wait 15 completed || true
        echo "FIS experiment \$eid \$state, waiting for the interruption"
    }

    export AWS_DEFAULT_REGION=\$(imds placement/region)
    for i in \$(seq 1 ${params.steps}); do
        echo "\$i \$(imds instance-id)" >> count.txt
        if [ \$((i % 30)) -eq 0 ]; then echo "step \$i on \$(imds instance-id)"; fi
        # a restore resumes after this point; a restart from scratch is a new Batch attempt and must not interrupt again
        if [ \$i -eq ${params.steps.intdiv(2)} ] && [ "\${AWS_BATCH_JOB_ATTEMPT:-1}" -eq 1 ]; then
            interrupt_self
        fi
        sleep 1
    done
    """
}

workflow {
    COUNT()
}

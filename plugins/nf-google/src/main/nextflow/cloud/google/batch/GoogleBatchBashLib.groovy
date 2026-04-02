/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package nextflow.cloud.google.batch

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import nextflow.cloud.google.batch.client.BatchConfig
import nextflow.executor.BashFunLib
import nextflow.util.Escape

/**
 * Bash helpers for Google Cloud Storage transfers in Google Batch tasks.
 * Order within each transport is controlled by {@code NXF_STAGE_IN_COPY_TRANSPORT} / {@code NXF_STAGE_OUT_COPY_TRANSPORT}.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class GoogleBatchBashLib extends BashFunLib<GoogleBatchBashLib> {

    private String gcloudCli = 'gcloud'

    private String gsutilCli = 'gsutil'

    private BatchConfig batchConfig

    GoogleBatchBashLib withGcloudCli(String value) {
        if( value )
            this.gcloudCli = value
        return this
    }

    GoogleBatchBashLib withGsutilCli(String value) {
        if( value )
            this.gsutilCli = value
        return this
    }

    GoogleBatchBashLib withBatchConfig(BatchConfig config) {
        this.batchConfig = config
        return this
    }

    protected String gsLib() {
        final gcloud = Escape.path(gcloudCli ?: 'gcloud')
        final gsutil = Escape.path(gsutilCli ?: 'gsutil')
        final mountRoot = GoogleBatchScriptLauncher.MOUNT_ROOT
        """
        export NXF_GCLOUD=$gcloud
        export NXF_GSUTIL=$gsutil
        export NXF_GS_MOUNT_ROOT="$mountRoot"

        nxf_gs_uri_to_mount() {
            case "\$1" in
            gs://*)
                local p=\${1#gs://}
                local bucket=\${p%%/*}
                if [[ "\$bucket" == "\$p" ]]; then
                    echo "\$NXF_GS_MOUNT_ROOT/\$bucket"
                else
                    local rest=\${p#*/}
                    echo "\$NXF_GS_MOUNT_ROOT/\$bucket/\$rest"
                fi
                ;;
            *)
                echo ""
                ;;
            esac
        }

        nxf_gs_download_gcloud() {
            command -v "\$NXF_GCLOUD" >/dev/null 2>&1 || return 1
            local source=\$1
            local target=\$2
            local basedir=\$(dirname "\$2")
            mkdir -p "\$basedir"
            if \$NXF_GCLOUD storage cp "\$source" "\$target" 2>/dev/null; then
                return 0
            fi
            mkdir -p "\$target"
            if \$NXF_GCLOUD storage cp --recursive "\$source" "\$target" >/dev/null 2>&1; then
                return 0
            fi
            rm -rf "\$target"
            return 1
        }

        nxf_gs_download_gsutil() {
            command -v "\$NXF_GSUTIL" >/dev/null 2>&1 || return 1
            local source=\$1
            local target=\$2
            local basedir=\$(dirname "\$2")
            mkdir -p "\$basedir"
            local ret
            ret=\$(\$NXF_GSUTIL cp "\$source" "\$target" 2>&1) || {
                mkdir -p "\$target"
                \$NXF_GSUTIL -m cp -r "\$source" "\$target" >/dev/null || {
                    rm -rf "\$target"
                    return 1
                }
            }
            return 0
        }

        nxf_gs_download_mount() {
            local source=\$1
            local target=\$2
            local ms
            ms=\$(nxf_gs_uri_to_mount "\$source")
            [[ -n "\$ms" ]] || return 1
            [[ -e "\$ms" || -d "\$ms" ]] || return 1
            mkdir -p "\$(dirname "\$target")"
            cp -fRL "\$ms" "\$target"
        }

        nxf_gs_download() {
            local source=\$1
            local target=\$2
            case "\${NXF_STAGE_IN_COPY_TRANSPORT:-posix}" in
            gsutil)
                if nxf_gs_download_gsutil "\$source" "\$target"; then return 0; fi
                if nxf_gs_download_gcloud "\$source" "\$target"; then return 0; fi
                if nxf_gs_download_mount "\$source" "\$target"; then return 0; fi
                ;;
            gcloud)
                if nxf_gs_download_gcloud "\$source" "\$target"; then return 0; fi
                if nxf_gs_download_gsutil "\$source" "\$target"; then return 0; fi
                if nxf_gs_download_mount "\$source" "\$target"; then return 0; fi
                ;;
            posix|*)
                if nxf_gs_download_mount "\$source" "\$target"; then return 0; fi
                ;;
            esac
            >&2 echo "Unable to download path: \$source"
            exit 1
        }

        nxf_gs_upload_try_gcloud() {
            command -v "\$NXF_GCLOUD" >/dev/null 2>&1 || return 1
            local name=\$1
            local gspath=\$2
            if [[ "\$name" == '-' ]]; then
                return 1
            fi
            if [[ -d "\$name" ]]; then
                \$NXF_GCLOUD storage cp --recursive "\$name" "\$gspath/\$name"
            else
                \$NXF_GCLOUD storage cp "\$name" "\$gspath/\$name"
            fi
        }

        nxf_gs_upload_try_gsutil() {
            command -v "\$NXF_GSUTIL" >/dev/null 2>&1 || return 1
            local name=\$1
            local gspath=\$2
            if [[ "\$name" == '-' ]]; then
                return 1
            fi
            if [[ -d "\$name" ]]; then
                \$NXF_GSUTIL -m cp -r "\$name" "\$gspath/\$name"
            else
                \$NXF_GSUTIL cp "\$name" "\$gspath/\$name"
            fi
        }

        nxf_gs_upload_try_mount() {
            local name=\$1
            local gspath=\$2
            if [[ "\$name" == '-' ]]; then
                return 1
            fi
            local dest
            dest=\$(nxf_gs_uri_to_mount "\$gspath/\$name")
            [[ -n "\$dest" ]] || return 1
            mkdir -p "\$(dirname "\$dest")"
            cp -fRL "\$name" "\$dest"
        }

        nxf_gs_upload() {
            local name=\$1
            local gspath=\${2%/}
            local move=\${3:-0}
            if [[ "\$name" == '-' ]]; then
                if command -v "\$NXF_GCLOUD" >/dev/null 2>&1 && \$NXF_GCLOUD storage cp - "\$gspath"; then
                    :
                elif command -v "\$NXF_GSUTIL" >/dev/null 2>&1 && \$NXF_GSUTIL cp - "\$gspath"; then
                    :
                else
                    return 1
                fi
            else
                case "\${NXF_STAGE_OUT_COPY_TRANSPORT:-posix}" in
                gsutil)
                    if nxf_gs_upload_try_gsutil "\$name" "\$gspath"; then :;
                    elif nxf_gs_upload_try_gcloud "\$name" "\$gspath"; then :;
                    elif nxf_gs_upload_try_mount "\$name" "\$gspath"; then :;
                    else
                        >&2 echo "Unable to upload path: \$name"
                        return 1
                    fi
                    ;;
                gcloud)
                    if nxf_gs_upload_try_gcloud "\$name" "\$gspath"; then :;
                    elif nxf_gs_upload_try_gsutil "\$name" "\$gspath"; then :;
                    elif nxf_gs_upload_try_mount "\$name" "\$gspath"; then :;
                    else
                        >&2 echo "Unable to upload path: \$name"
                        return 1
                    fi
                    ;;
                posix|*)
                    if nxf_gs_upload_try_mount "\$name" "\$gspath"; then :;
                    else
                        >&2 echo "Unable to upload path: \$name"
                        return 1
                    fi
                    ;;
                esac
            fi
            if [[ "\$move" == "1" && "\$name" != '-' ]]; then
                rm -rf "\$name"
            fi
        }
        """.stripIndent(true)
    }

    protected String transportExports() {
        final cfg = batchConfig
        if( !cfg )
            return ''
        final inTr = cfg.stageInCopyTransport ?: BatchConfig.COPY_TRANSPORT_POSIX
        final outTr = cfg.stageOutCopyTransport ?: BatchConfig.COPY_TRANSPORT_POSIX
        """
        export NXF_STAGE_IN_COPY_TRANSPORT=${Escape.path(inTr)}
        export NXF_STAGE_OUT_COPY_TRANSPORT=${Escape.path(outTr)}
        """.stripIndent(true)
    }

    @Override
    String render() {
        super.render() + transportExports() + gsLib()
    }

    @Memoized
    static String script(BatchConfig config) {
        new GoogleBatchBashLib()
                .includeCoreFun(true)
                .withMaxParallelTransfers(config.maxParallelTransfers)
                .withMaxTransferAttempts(config.maxTransferAttempts)
                .withDelayBetweenAttempts(config.delayBetweenAttempts)
                .withGcloudCli(config.gcloudCli)
                .withGsutilCli(config.gsutilCli)
                .withBatchConfig(config)
                .render()
    }
}

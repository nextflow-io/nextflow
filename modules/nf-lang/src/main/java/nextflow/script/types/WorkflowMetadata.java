/*
 * Copyright 2024-2025, Seqera Labs
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
package nextflow.script.types;

import java.nio.file.Path;
import java.time.OffsetDateTime;
import java.util.List;
import java.util.Map;
import java.util.UUID;

import groovy.lang.Closure;
import nextflow.script.dsl.Constant;
import nextflow.script.dsl.Description;

public interface WorkflowMetadata {

    @Constant("commandLine")
    @Description("""
        Command line as entered by the user to launch the workflow execution.
    """)
    String getCommandLine();

    @Constant("commitId")
    @Description("""
        Git commit ID of the executed workflow repository.
    """)
    String getCommitId();

    @Constant("complete")
    @Description("""
        Timestamp of workflow when execution is completed.
    """)
    OffsetDateTime getComplete();

    @Constant("configFiles")
    @Description("""
        Configuration files used for the workflow execution.
    """)
    List<Path> getConfigFiles();

    @Constant("container")
    @Description("""
        Docker image used to run workflow tasks, or a map of process names to process containers when multiple images are used.
    """)
    /* String | Map<String,String> */
    Map<String,String> getContainer();

    @Constant("containerEngine")
    @Description("""
        Returns the name of the container engine (e.g. docker or singularity) or null if no container engine is enabled.
    """)
    String getContainerEngine();

    @Constant("duration")
    @Description("""
        Time elapsed to complete workflow execution.
    """)
    Duration getDuration();

    @Constant("errorMessage")
    @Description("""
        Error message of the task that caused the workflow execution to fail.
    """)
    String getErrorMessage();

    @Constant("errorReport")
    @Description("""
        Detailed error of the task that caused the workflow execution to fail.
    """)
    String getErrorReport();

    @Constant("exitStatus")
    @Description("""
        Exit status of the task that caused the workflow execution to fail.
    """)
    int getExitStatus();

    @Constant("failOnIgnore")
    @Description("""
        Whether the `workflow.failOnIgnore` config option was enabled.
    """)
    boolean isFailOnIgnore();

    @Constant("fusion")
    @Description("""
        Map of Fusion runtime information.
    """)
    FusionMetadata getFusion();

    @Constant("homeDir")
    @Description("""
        User system home directory.
    """)
    Path getHomeDir();

    @Constant("launchDir")
    @Description("""
        Directory where the workflow was launched.
    """)
    Path getLaunchDir();

    @Constant("manifest")
    @Description("""
        Map of properties corresponding to the {ref}`config-manifest` config scope.
    """)
    Manifest getManifest();

    @Constant("outputDir")
    @Description("""
        Workflow output directory.
    """)
    Path getOutputDir();

    @Constant("preview")
    @Description("""
        Whether the current workflow run is a preview run.
    """)
    boolean isPreview();

    @Constant("profile")
    @Description("""
        Comma-separated list of active configuration profiles.
    """)
    String getProfile();

    @Constant("projectDir")
    @Description("""
        Directory where the workflow project is located.
    """)
    Path getProjectDir();

    @Constant("repository")
    @Description("""
        Project repository Git remote URL.
    """)
    String getRepository();

    @Constant("resume")
    @Description("""
        Whether the current instance is resumed from a previous execution.
    """)
    boolean isResume();

    @Constant("revision")
    @Description("""
        Git branch/tag of the executed workflow repository.
    """)
    String getRevision();

    @Constant("runName")
    @Description("""
        Mnemonic name assigned to this execution instance.
    """)
    String getRunName();

    @Constant("scriptFile")
    @Description("""
        Project main script file path.
    """)
    Path getScriptFile();

    @Constant("scriptId")
    @Description("""
        Project main script unique hash ID.
    """)
    String getScriptId();

    @Constant("scriptName")
    @Description("""
        Project main script file name.
    """)
    String getScriptName();

    @Constant("sessionId")
    @Description("""
        Unique identifier (UUID) associated to current execution.
    """)
    UUID getSessionId();

    @Constant("start")
    @Description("""
        Timestamp of workflow at execution start.
    """)
    OffsetDateTime getStart();

    @Constant("stubRun")
    @Description("""
        Whether the current instance is a stub-run execution.
    """)
    boolean isStubRun();

    @Constant("success")
    @Description("""
        Whether the execution completed successfully.
    """)
    boolean isSuccess();

    @Constant("userName")
    @Description("""
        User system account name.
    """)
    String getUserName();

    @Constant("wave")
    @Description("""
        Map of Wave runtime information.
    """)
    WaveMetadata getWave();

    @Constant("workDir")
    @Description("""
        The directory where task temporary files are stored.
    """)
    Path getWorkDir();

    @Description("""
        Define an action to take when the workflow completes (whether successful or not).
    """)
    void onComplete(Closure action);

    @Description("""
        Define an action to take if the workflow is terminated due to a runtime error or task failure.
    """)
    void onError(Closure action);

    interface FusionMetadata {

        @Constant("enabled")
        @Description("""
            Whether Fusion is enabled.
        """)
        boolean isEnabled();

        @Constant("version")
        @Description("""
            The Fusion version being used.
        """)
        String getVersion();
    }

    interface WaveMetadata {

        @Constant("enabled")
        @Description("""
            Whether Wave is enabled.
        """)
        boolean isEnabled();
    }

}

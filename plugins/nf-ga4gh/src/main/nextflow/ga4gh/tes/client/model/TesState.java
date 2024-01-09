/*
 * Copyright 2013-2023, Seqera Labs
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

/*
 * task_execution.proto
 * No description provided (generated by Swagger Codegen https://github.com/swagger-api/swagger-codegen)
 *
 * OpenAPI spec version: version not set
 *
 *
 * NOTE: This class is auto generated by the swagger code generator program.
 * https://github.com/swagger-api/swagger-codegen.git
 * Do not edit the class manually.
 */


package nextflow.ga4gh.tes.client.model;

import java.io.IOException;
import com.google.gson.TypeAdapter;
import com.google.gson.annotations.JsonAdapter;
import com.google.gson.stream.JsonReader;
import com.google.gson.stream.JsonWriter;

/**
 * Task states.   - UNKNOWN: The state of the task is unknown.  This provides a safe default for messages where this field is missing, for example, so that a missing field does not accidentally imply that the state is QUEUED.  - QUEUED: The task is queued.  - INITIALIZING: The task has been assigned to a worker and is currently preparing to run. For example, the worker may be turning on, downloading input files, etc.  - RUNNING: The task is running. Input files are downloaded and the first Executor has been started.  - PAUSED: The task is paused.  An implementation may have the ability to pause a task, but this is not required.  - COMPLETE: The task has completed running. Executors have exited without error and output files have been successfully uploaded.  - EXECUTOR_ERROR: The task encountered an error in one of the Executor processes. Generally, this means that an Executor exited with a non-zero exit code.  - SYSTEM_ERROR: The task was stopped due to a system error, but not from an Executor, for example an upload failed due to network issues, the worker&#39;s ran out of disk space, etc.  - CANCELED: The task was canceled by the user.
 */
@JsonAdapter(TesState.Adapter.class)
public enum TesState {

  UNKNOWN("UNKNOWN"),

  QUEUED("QUEUED"),

  INITIALIZING("INITIALIZING"),

  RUNNING("RUNNING"),

  PAUSED("PAUSED"),

  COMPLETE("COMPLETE"),

  EXECUTOR_ERROR("EXECUTOR_ERROR"),

  SYSTEM_ERROR("SYSTEM_ERROR"),

  CANCELED("CANCELED");

  private String value;

  TesState(String value) {
    this.value = value;
  }

  public String getValue() {
    return value;
  }

  @Override
  public String toString() {
    return String.valueOf(value);
  }

  public static TesState fromValue(String text) {
    for (TesState b : TesState.values()) {
      if (String.valueOf(b.value).equals(text)) {
        return b;
      }
    }
    return null;
  }

  public static class Adapter extends TypeAdapter<TesState> {
    @Override
    public void write(final JsonWriter jsonWriter, final TesState enumeration) throws IOException {
      jsonWriter.value(enumeration.getValue());
    }

    @Override
    public TesState read(final JsonReader jsonReader) throws IOException {
      String value = jsonReader.nextString();
      return TesState.fromValue(String.valueOf(value));
    }
  }
}


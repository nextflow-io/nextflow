package nextflow.ga4gh.tes.server.model;

import java.util.Objects;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import nextflow.ga4gh.tes.server.model.TesExecutor;
import nextflow.ga4gh.tes.server.model.TesInput;
import nextflow.ga4gh.tes.server.model.TesOutput;
import nextflow.ga4gh.tes.server.model.TesResources;
import nextflow.ga4gh.tes.server.model.TesState;
import nextflow.ga4gh.tes.server.model.TesTaskLog;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Task describes an instance of a task.
 **/
@JsonInclude(JsonInclude.Include.NON_NULL) 
public class TesTask   {
  
  private String id = null;
  private TesState state = null;
  private String name = null;
  private String description = null;
  private List<TesInput> inputs = new ArrayList<TesInput>();
  private List<TesOutput> outputs = new ArrayList<TesOutput>();
  private TesResources resources = null;
  private List<TesExecutor> executors = new ArrayList<TesExecutor>();
  private List<String> volumes = new ArrayList<String>();
  private Map<String, String> tags = new HashMap<String, String>();
  private List<TesTaskLog> logs = new ArrayList<TesTaskLog>();
  private String creationTime = null;

  public TesTask () {

  }

  public TesTask (String id, TesState state, String name, String description, List<TesInput> inputs, List<TesOutput> outputs, TesResources resources, List<TesExecutor> executors, List<String> volumes, Map<String, String> tags, List<TesTaskLog> logs, String creationTime) {
    this.id = id;
    this.state = state;
    this.name = name;
    this.description = description;
    this.inputs = inputs;
    this.outputs = outputs;
    this.resources = resources;
    this.executors = executors;
    this.volumes = volumes;
    this.tags = tags;
    this.logs = logs;
    this.creationTime = creationTime;
  }

    
  @JsonProperty("id")
  public String getId() {
    return id;
  }
  public void setId(String id) {
    this.id = id;
  }

    
  @JsonProperty("state")
  public TesState getState() {
    return state;
  }
  public void setState(TesState state) {
    this.state = state;
  }

    
  @JsonProperty("name")
  public String getName() {
    return name;
  }
  public void setName(String name) {
    this.name = name;
  }

    
  @JsonProperty("description")
  public String getDescription() {
    return description;
  }
  public void setDescription(String description) {
    this.description = description;
  }

    
  @JsonProperty("inputs")
  public List<TesInput> getInputs() {
    return inputs;
  }
  public void setInputs(List<TesInput> inputs) {
    this.inputs = inputs;
  }

    
  @JsonProperty("outputs")
  public List<TesOutput> getOutputs() {
    return outputs;
  }
  public void setOutputs(List<TesOutput> outputs) {
    this.outputs = outputs;
  }

    
  @JsonProperty("resources")
  public TesResources getResources() {
    return resources;
  }
  public void setResources(TesResources resources) {
    this.resources = resources;
  }

    
  @JsonProperty("executors")
  public List<TesExecutor> getExecutors() {
    return executors;
  }
  public void setExecutors(List<TesExecutor> executors) {
    this.executors = executors;
  }

    
  @JsonProperty("volumes")
  public List<String> getVolumes() {
    return volumes;
  }
  public void setVolumes(List<String> volumes) {
    this.volumes = volumes;
  }

    
  @JsonProperty("tags")
  public Map<String, String> getTags() {
    return tags;
  }
  public void setTags(Map<String, String> tags) {
    this.tags = tags;
  }

    
  @JsonProperty("logs")
  public List<TesTaskLog> getLogs() {
    return logs;
  }
  public void setLogs(List<TesTaskLog> logs) {
    this.logs = logs;
  }

    
  @JsonProperty("creation_time")
  public String getCreationTime() {
    return creationTime;
  }
  public void setCreationTime(String creationTime) {
    this.creationTime = creationTime;
  }


  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    TesTask tesTask = (TesTask) o;
    return Objects.equals(id, tesTask.id) &&
        Objects.equals(state, tesTask.state) &&
        Objects.equals(name, tesTask.name) &&
        Objects.equals(description, tesTask.description) &&
        Objects.equals(inputs, tesTask.inputs) &&
        Objects.equals(outputs, tesTask.outputs) &&
        Objects.equals(resources, tesTask.resources) &&
        Objects.equals(executors, tesTask.executors) &&
        Objects.equals(volumes, tesTask.volumes) &&
        Objects.equals(tags, tesTask.tags) &&
        Objects.equals(logs, tesTask.logs) &&
        Objects.equals(creationTime, tesTask.creationTime);
  }

  @Override
  public int hashCode() {
    return Objects.hash(id, state, name, description, inputs, outputs, resources, executors, volumes, tags, logs, creationTime);
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("class TesTask {\n");
    
    sb.append("    id: ").append(toIndentedString(id)).append("\n");
    sb.append("    state: ").append(toIndentedString(state)).append("\n");
    sb.append("    name: ").append(toIndentedString(name)).append("\n");
    sb.append("    description: ").append(toIndentedString(description)).append("\n");
    sb.append("    inputs: ").append(toIndentedString(inputs)).append("\n");
    sb.append("    outputs: ").append(toIndentedString(outputs)).append("\n");
    sb.append("    resources: ").append(toIndentedString(resources)).append("\n");
    sb.append("    executors: ").append(toIndentedString(executors)).append("\n");
    sb.append("    volumes: ").append(toIndentedString(volumes)).append("\n");
    sb.append("    tags: ").append(toIndentedString(tags)).append("\n");
    sb.append("    logs: ").append(toIndentedString(logs)).append("\n");
    sb.append("    creationTime: ").append(toIndentedString(creationTime)).append("\n");
    sb.append("}");
    return sb.toString();
  }

  /**
   * Convert the given object to string with each line indented by 4 spaces
   * (except the first line).
   */
  private String toIndentedString(Object o) {
    if (o == null) {
      return "null";
    }
    return o.toString().replace("\n", "\n    ");
  }
}

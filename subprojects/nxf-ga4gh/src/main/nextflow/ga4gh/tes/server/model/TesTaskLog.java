package nextflow.ga4gh.tes.server.model;

import java.util.Objects;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import nextflow.ga4gh.tes.server.model.TesExecutorLog;
import nextflow.ga4gh.tes.server.model.TesOutputFileLog;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * TaskLog describes logging information related to a Task.
 **/
@JsonInclude(JsonInclude.Include.NON_NULL) 
public class TesTaskLog   {
  
  private List<TesExecutorLog> logs = new ArrayList<TesExecutorLog>();
  private Map<String, String> metadata = new HashMap<String, String>();
  private String startTime = null;
  private String endTime = null;
  private List<TesOutputFileLog> outputs = new ArrayList<TesOutputFileLog>();
  private List<String> systemLogs = new ArrayList<String>();

  public TesTaskLog () {

  }

  public TesTaskLog (List<TesExecutorLog> logs, Map<String, String> metadata, String startTime, String endTime, List<TesOutputFileLog> outputs, List<String> systemLogs) {
    this.logs = logs;
    this.metadata = metadata;
    this.startTime = startTime;
    this.endTime = endTime;
    this.outputs = outputs;
    this.systemLogs = systemLogs;
  }

    
  @JsonProperty("logs")
  public List<TesExecutorLog> getLogs() {
    return logs;
  }
  public void setLogs(List<TesExecutorLog> logs) {
    this.logs = logs;
  }

    
  @JsonProperty("metadata")
  public Map<String, String> getMetadata() {
    return metadata;
  }
  public void setMetadata(Map<String, String> metadata) {
    this.metadata = metadata;
  }

    
  @JsonProperty("start_time")
  public String getStartTime() {
    return startTime;
  }
  public void setStartTime(String startTime) {
    this.startTime = startTime;
  }

    
  @JsonProperty("end_time")
  public String getEndTime() {
    return endTime;
  }
  public void setEndTime(String endTime) {
    this.endTime = endTime;
  }

    
  @JsonProperty("outputs")
  public List<TesOutputFileLog> getOutputs() {
    return outputs;
  }
  public void setOutputs(List<TesOutputFileLog> outputs) {
    this.outputs = outputs;
  }

    
  @JsonProperty("system_logs")
  public List<String> getSystemLogs() {
    return systemLogs;
  }
  public void setSystemLogs(List<String> systemLogs) {
    this.systemLogs = systemLogs;
  }


  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    TesTaskLog tesTaskLog = (TesTaskLog) o;
    return Objects.equals(logs, tesTaskLog.logs) &&
        Objects.equals(metadata, tesTaskLog.metadata) &&
        Objects.equals(startTime, tesTaskLog.startTime) &&
        Objects.equals(endTime, tesTaskLog.endTime) &&
        Objects.equals(outputs, tesTaskLog.outputs) &&
        Objects.equals(systemLogs, tesTaskLog.systemLogs);
  }

  @Override
  public int hashCode() {
    return Objects.hash(logs, metadata, startTime, endTime, outputs, systemLogs);
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("class TesTaskLog {\n");
    
    sb.append("    logs: ").append(toIndentedString(logs)).append("\n");
    sb.append("    metadata: ").append(toIndentedString(metadata)).append("\n");
    sb.append("    startTime: ").append(toIndentedString(startTime)).append("\n");
    sb.append("    endTime: ").append(toIndentedString(endTime)).append("\n");
    sb.append("    outputs: ").append(toIndentedString(outputs)).append("\n");
    sb.append("    systemLogs: ").append(toIndentedString(systemLogs)).append("\n");
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

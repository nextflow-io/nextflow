package nextflow.ga4gh.tes.server.model;

import java.util.Objects;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;

/**
 * ExecutorLog describes logging information related to an Executor.
 **/
@JsonInclude(JsonInclude.Include.NON_NULL) 
public class TesExecutorLog   {
  
  private String startTime = null;
  private String endTime = null;
  private String stdout = null;
  private String stderr = null;
  private Integer exitCode = null;

  public TesExecutorLog () {

  }

  public TesExecutorLog (String startTime, String endTime, String stdout, String stderr, Integer exitCode) {
    this.startTime = startTime;
    this.endTime = endTime;
    this.stdout = stdout;
    this.stderr = stderr;
    this.exitCode = exitCode;
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

    
  @JsonProperty("stdout")
  public String getStdout() {
    return stdout;
  }
  public void setStdout(String stdout) {
    this.stdout = stdout;
  }

    
  @JsonProperty("stderr")
  public String getStderr() {
    return stderr;
  }
  public void setStderr(String stderr) {
    this.stderr = stderr;
  }

    
  @JsonProperty("exit_code")
  public Integer getExitCode() {
    return exitCode;
  }
  public void setExitCode(Integer exitCode) {
    this.exitCode = exitCode;
  }


  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    TesExecutorLog tesExecutorLog = (TesExecutorLog) o;
    return Objects.equals(startTime, tesExecutorLog.startTime) &&
        Objects.equals(endTime, tesExecutorLog.endTime) &&
        Objects.equals(stdout, tesExecutorLog.stdout) &&
        Objects.equals(stderr, tesExecutorLog.stderr) &&
        Objects.equals(exitCode, tesExecutorLog.exitCode);
  }

  @Override
  public int hashCode() {
    return Objects.hash(startTime, endTime, stdout, stderr, exitCode);
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("class TesExecutorLog {\n");
    
    sb.append("    startTime: ").append(toIndentedString(startTime)).append("\n");
    sb.append("    endTime: ").append(toIndentedString(endTime)).append("\n");
    sb.append("    stdout: ").append(toIndentedString(stdout)).append("\n");
    sb.append("    stderr: ").append(toIndentedString(stderr)).append("\n");
    sb.append("    exitCode: ").append(toIndentedString(exitCode)).append("\n");
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

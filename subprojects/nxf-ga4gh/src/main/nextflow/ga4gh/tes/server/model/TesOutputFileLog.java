package nextflow.ga4gh.tes.server.model;

import java.util.Objects;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;

/**
 * OutputFileLog describes a single output file. This describes file details after the task has completed successfully, for logging purposes.
 **/
@JsonInclude(JsonInclude.Include.NON_NULL) 
public class TesOutputFileLog   {
  
  private String url = null;
  private String path = null;
  private String sizeBytes = null;

  public TesOutputFileLog () {

  }

  public TesOutputFileLog (String url, String path, String sizeBytes) {
    this.url = url;
    this.path = path;
    this.sizeBytes = sizeBytes;
  }

    
  @JsonProperty("url")
  public String getUrl() {
    return url;
  }
  public void setUrl(String url) {
    this.url = url;
  }

    
  @JsonProperty("path")
  public String getPath() {
    return path;
  }
  public void setPath(String path) {
    this.path = path;
  }

    
  @JsonProperty("size_bytes")
  public String getSizeBytes() {
    return sizeBytes;
  }
  public void setSizeBytes(String sizeBytes) {
    this.sizeBytes = sizeBytes;
  }


  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    TesOutputFileLog tesOutputFileLog = (TesOutputFileLog) o;
    return Objects.equals(url, tesOutputFileLog.url) &&
        Objects.equals(path, tesOutputFileLog.path) &&
        Objects.equals(sizeBytes, tesOutputFileLog.sizeBytes);
  }

  @Override
  public int hashCode() {
    return Objects.hash(url, path, sizeBytes);
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("class TesOutputFileLog {\n");
    
    sb.append("    url: ").append(toIndentedString(url)).append("\n");
    sb.append("    path: ").append(toIndentedString(path)).append("\n");
    sb.append("    sizeBytes: ").append(toIndentedString(sizeBytes)).append("\n");
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

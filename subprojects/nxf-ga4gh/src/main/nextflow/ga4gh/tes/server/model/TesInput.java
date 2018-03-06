package nextflow.ga4gh.tes.server.model;

import java.util.Objects;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import nextflow.ga4gh.tes.server.model.TesFileType;

/**
 * Input describes Task input files.
 **/
@JsonInclude(JsonInclude.Include.NON_NULL) 
public class TesInput   {
  
  private String name = null;
  private String description = null;
  private String url = null;
  private String path = null;
  private TesFileType type = null;
  private String content = null;

  public TesInput () {

  }

  public TesInput (String name, String description, String url, String path, TesFileType type, String content) {
    this.name = name;
    this.description = description;
    this.url = url;
    this.path = path;
    this.type = type;
    this.content = content;
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

    
  @JsonProperty("type")
  public TesFileType getType() {
    return type;
  }
  public void setType(TesFileType type) {
    this.type = type;
  }

    
  @JsonProperty("content")
  public String getContent() {
    return content;
  }
  public void setContent(String content) {
    this.content = content;
  }


  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    TesInput tesInput = (TesInput) o;
    return Objects.equals(name, tesInput.name) &&
        Objects.equals(description, tesInput.description) &&
        Objects.equals(url, tesInput.url) &&
        Objects.equals(path, tesInput.path) &&
        Objects.equals(type, tesInput.type) &&
        Objects.equals(content, tesInput.content);
  }

  @Override
  public int hashCode() {
    return Objects.hash(name, description, url, path, type, content);
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("class TesInput {\n");
    
    sb.append("    name: ").append(toIndentedString(name)).append("\n");
    sb.append("    description: ").append(toIndentedString(description)).append("\n");
    sb.append("    url: ").append(toIndentedString(url)).append("\n");
    sb.append("    path: ").append(toIndentedString(path)).append("\n");
    sb.append("    type: ").append(toIndentedString(type)).append("\n");
    sb.append("    content: ").append(toIndentedString(content)).append("\n");
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

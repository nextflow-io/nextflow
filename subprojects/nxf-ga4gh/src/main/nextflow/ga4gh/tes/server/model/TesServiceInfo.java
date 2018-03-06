package nextflow.ga4gh.tes.server.model;

import java.util.Objects;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import java.util.ArrayList;
import java.util.List;

/**
 * ServiceInfo describes information about the service, such as storage details, resource availability, and other documentation.
 **/
@JsonInclude(JsonInclude.Include.NON_NULL) 
public class TesServiceInfo   {
  
  private String name = null;
  private String doc = null;
  private List<String> storage = new ArrayList<String>();

  public TesServiceInfo () {

  }

  public TesServiceInfo (String name, String doc, List<String> storage) {
    this.name = name;
    this.doc = doc;
    this.storage = storage;
  }

    
  @JsonProperty("name")
  public String getName() {
    return name;
  }
  public void setName(String name) {
    this.name = name;
  }

    
  @JsonProperty("doc")
  public String getDoc() {
    return doc;
  }
  public void setDoc(String doc) {
    this.doc = doc;
  }

    
  @JsonProperty("storage")
  public List<String> getStorage() {
    return storage;
  }
  public void setStorage(List<String> storage) {
    this.storage = storage;
  }


  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    TesServiceInfo tesServiceInfo = (TesServiceInfo) o;
    return Objects.equals(name, tesServiceInfo.name) &&
        Objects.equals(doc, tesServiceInfo.doc) &&
        Objects.equals(storage, tesServiceInfo.storage);
  }

  @Override
  public int hashCode() {
    return Objects.hash(name, doc, storage);
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("class TesServiceInfo {\n");
    
    sb.append("    name: ").append(toIndentedString(name)).append("\n");
    sb.append("    doc: ").append(toIndentedString(doc)).append("\n");
    sb.append("    storage: ").append(toIndentedString(storage)).append("\n");
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

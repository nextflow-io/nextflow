package nextflow.ga4gh.tes.server.model;

import java.util.Objects;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import java.util.ArrayList;
import java.util.List;

/**
 * Resources describes the resources requested by a task.
 **/
@JsonInclude(JsonInclude.Include.NON_NULL) 
public class TesResources   {
  
  private Long cpuCores = null;
  private Boolean preemptible = null;
  private Double ramGb = null;
  private Double diskGb = null;
  private List<String> zones = new ArrayList<String>();

  public TesResources () {

  }

  public TesResources (Long cpuCores, Boolean preemptible, Double ramGb, Double diskGb, List<String> zones) {
    this.cpuCores = cpuCores;
    this.preemptible = preemptible;
    this.ramGb = ramGb;
    this.diskGb = diskGb;
    this.zones = zones;
  }

    
  @JsonProperty("cpu_cores")
  public Long getCpuCores() {
    return cpuCores;
  }
  public void setCpuCores(Long cpuCores) {
    this.cpuCores = cpuCores;
  }

    
  @JsonProperty("preemptible")
  public Boolean getPreemptible() {
    return preemptible;
  }
  public void setPreemptible(Boolean preemptible) {
    this.preemptible = preemptible;
  }

    
  @JsonProperty("ram_gb")
  public Double getRamGb() {
    return ramGb;
  }
  public void setRamGb(Double ramGb) {
    this.ramGb = ramGb;
  }

    
  @JsonProperty("disk_gb")
  public Double getDiskGb() {
    return diskGb;
  }
  public void setDiskGb(Double diskGb) {
    this.diskGb = diskGb;
  }

    
  @JsonProperty("zones")
  public List<String> getZones() {
    return zones;
  }
  public void setZones(List<String> zones) {
    this.zones = zones;
  }


  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    TesResources tesResources = (TesResources) o;
    return Objects.equals(cpuCores, tesResources.cpuCores) &&
        Objects.equals(preemptible, tesResources.preemptible) &&
        Objects.equals(ramGb, tesResources.ramGb) &&
        Objects.equals(diskGb, tesResources.diskGb) &&
        Objects.equals(zones, tesResources.zones);
  }

  @Override
  public int hashCode() {
    return Objects.hash(cpuCores, preemptible, ramGb, diskGb, zones);
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("class TesResources {\n");
    
    sb.append("    cpuCores: ").append(toIndentedString(cpuCores)).append("\n");
    sb.append("    preemptible: ").append(toIndentedString(preemptible)).append("\n");
    sb.append("    ramGb: ").append(toIndentedString(ramGb)).append("\n");
    sb.append("    diskGb: ").append(toIndentedString(diskGb)).append("\n");
    sb.append("    zones: ").append(toIndentedString(zones)).append("\n");
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

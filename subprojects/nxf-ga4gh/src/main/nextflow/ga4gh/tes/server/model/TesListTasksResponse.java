package nextflow.ga4gh.tes.server.model;

import java.util.Objects;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import nextflow.ga4gh.tes.server.model.TesTask;
import java.util.ArrayList;
import java.util.List;

/**
 * ListTasksResponse describes a response from the ListTasks endpoint.
 **/
@JsonInclude(JsonInclude.Include.NON_NULL) 
public class TesListTasksResponse   {
  
  private List<TesTask> tasks = new ArrayList<TesTask>();
  private String nextPageToken = null;

  public TesListTasksResponse () {

  }

  public TesListTasksResponse (List<TesTask> tasks, String nextPageToken) {
    this.tasks = tasks;
    this.nextPageToken = nextPageToken;
  }

    
  @JsonProperty("tasks")
  public List<TesTask> getTasks() {
    return tasks;
  }
  public void setTasks(List<TesTask> tasks) {
    this.tasks = tasks;
  }

    
  @JsonProperty("next_page_token")
  public String getNextPageToken() {
    return nextPageToken;
  }
  public void setNextPageToken(String nextPageToken) {
    this.nextPageToken = nextPageToken;
  }


  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    TesListTasksResponse tesListTasksResponse = (TesListTasksResponse) o;
    return Objects.equals(tasks, tesListTasksResponse.tasks) &&
        Objects.equals(nextPageToken, tesListTasksResponse.nextPageToken);
  }

  @Override
  public int hashCode() {
    return Objects.hash(tasks, nextPageToken);
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("class TesListTasksResponse {\n");
    
    sb.append("    tasks: ").append(toIndentedString(tasks)).append("\n");
    sb.append("    nextPageToken: ").append(toIndentedString(nextPageToken)).append("\n");
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

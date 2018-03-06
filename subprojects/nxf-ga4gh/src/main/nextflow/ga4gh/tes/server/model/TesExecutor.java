package nextflow.ga4gh.tes.server.model;

import java.util.Objects;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Executor describes a command to be executed, and its environment.
 **/
@JsonInclude(JsonInclude.Include.NON_NULL) 
public class TesExecutor   {
  
  private String image = null;
  private List<String> command = new ArrayList<String>();
  private String workdir = null;
  private String stdin = null;
  private String stdout = null;
  private String stderr = null;
  private Map<String, String> env = new HashMap<String, String>();

  public TesExecutor () {

  }

  public TesExecutor (String image, List<String> command, String workdir, String stdin, String stdout, String stderr, Map<String, String> env) {
    this.image = image;
    this.command = command;
    this.workdir = workdir;
    this.stdin = stdin;
    this.stdout = stdout;
    this.stderr = stderr;
    this.env = env;
  }

    
  @JsonProperty("image")
  public String getImage() {
    return image;
  }
  public void setImage(String image) {
    this.image = image;
  }

    
  @JsonProperty("command")
  public List<String> getCommand() {
    return command;
  }
  public void setCommand(List<String> command) {
    this.command = command;
  }

    
  @JsonProperty("workdir")
  public String getWorkdir() {
    return workdir;
  }
  public void setWorkdir(String workdir) {
    this.workdir = workdir;
  }

    
  @JsonProperty("stdin")
  public String getStdin() {
    return stdin;
  }
  public void setStdin(String stdin) {
    this.stdin = stdin;
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

    
  @JsonProperty("env")
  public Map<String, String> getEnv() {
    return env;
  }
  public void setEnv(Map<String, String> env) {
    this.env = env;
  }


  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    TesExecutor tesExecutor = (TesExecutor) o;
    return Objects.equals(image, tesExecutor.image) &&
        Objects.equals(command, tesExecutor.command) &&
        Objects.equals(workdir, tesExecutor.workdir) &&
        Objects.equals(stdin, tesExecutor.stdin) &&
        Objects.equals(stdout, tesExecutor.stdout) &&
        Objects.equals(stderr, tesExecutor.stderr) &&
        Objects.equals(env, tesExecutor.env);
  }

  @Override
  public int hashCode() {
    return Objects.hash(image, command, workdir, stdin, stdout, stderr, env);
  }

  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("class TesExecutor {\n");
    
    sb.append("    image: ").append(toIndentedString(image)).append("\n");
    sb.append("    command: ").append(toIndentedString(command)).append("\n");
    sb.append("    workdir: ").append(toIndentedString(workdir)).append("\n");
    sb.append("    stdin: ").append(toIndentedString(stdin)).append("\n");
    sb.append("    stdout: ").append(toIndentedString(stdout)).append("\n");
    sb.append("    stderr: ").append(toIndentedString(stderr)).append("\n");
    sb.append("    env: ").append(toIndentedString(env)).append("\n");
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

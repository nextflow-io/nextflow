/*
 * Copyright 2013-2024, Seqera Labs
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

import java.util.Objects;

import com.google.gson.annotations.SerializedName;
import io.swagger.annotations.ApiModel;
import io.swagger.annotations.ApiModelProperty;

import java.util.ArrayList;
import java.util.List;

/**
 * ServiceInfo describes information about the service, such as storage details, resource availability, and other documentation.
 */
@ApiModel(description = "ServiceInfo describes information about the service, such as storage details, resource availability, and other documentation.")
@javax.annotation.Generated(value = "io.swagger.codegen.languages.JavaClientCodegen", date = "2018-02-01T15:43:49.638Z")
public class TesServiceInfo {
  @SerializedName("name")
  private String name = null;

  @SerializedName("doc")
  private String doc = null;

  @SerializedName("storage")
  private List<String> storage = null;

  public TesServiceInfo name(String name) {
    this.name = name;
    return this;
  }

   /**
   * Returns the name of the service, e.g. \&quot;ohsu-compbio-funnel\&quot;.
   * @return name
  **/
  @ApiModelProperty(value = "Returns the name of the service, e.g. \"ohsu-compbio-funnel\".")
  public String getName() {
    return name;
  }

  public void setName(String name) {
    this.name = name;
  }

  public TesServiceInfo doc(String doc) {
    this.doc = doc;
    return this;
  }

   /**
   * Returns a documentation string, e.g. \&quot;Hey, we&#39;re OHSU Comp. Bio!\&quot;.
   * @return doc
  **/
  @ApiModelProperty(value = "Returns a documentation string, e.g. \"Hey, we're OHSU Comp. Bio!\".")
  public String getDoc() {
    return doc;
  }

  public void setDoc(String doc) {
    this.doc = doc;
  }

  public TesServiceInfo storage(List<String> storage) {
    this.storage = storage;
    return this;
  }

  public TesServiceInfo addStorageItem(String storageItem) {
    if (this.storage == null) {
      this.storage = new ArrayList<String>();
    }
    this.storage.add(storageItem);
    return this;
  }

   /**
   * Lists some, but not necessarily all, storage locations supported by the service.  Must be in a valid URL format. e.g.  file:///path/to/local/funnel-storage s3://ohsu-compbio-funnel/storage etc.
   * @return storage
  **/
  @ApiModelProperty(value = "Lists some, but not necessarily all, storage locations supported by the service.  Must be in a valid URL format. e.g.  file:///path/to/local/funnel-storage s3://ohsu-compbio-funnel/storage etc.")
  public List<String> getStorage() {
    return storage;
  }

  public void setStorage(List<String> storage) {
    this.storage = storage;
  }


  @Override
  public boolean equals(java.lang.Object o) {
    if (this == o) {
      return true;
    }
    if (o == null || getClass() != o.getClass()) {
      return false;
    }
    TesServiceInfo tesServiceInfo = (TesServiceInfo) o;
    return Objects.equals(this.name, tesServiceInfo.name) &&
        Objects.equals(this.doc, tesServiceInfo.doc) &&
        Objects.equals(this.storage, tesServiceInfo.storage);
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
  private String toIndentedString(java.lang.Object o) {
    if (o == null) {
      return "null";
    }
    return o.toString().replace("\n", "\n    ");
  }
  
}


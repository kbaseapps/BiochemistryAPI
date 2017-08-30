
package us.kbase.biochemistryapi;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.annotation.Generated;
import com.fasterxml.jackson.annotation.JsonAnyGetter;
import com.fasterxml.jackson.annotation.JsonAnySetter;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonPropertyOrder;


/**
 * <p>Original spec-file type: get_compounds_params</p>
 * <pre>
 * Input parameters for the "get_compounds" function.
 *     list<compound_id> compounds - a list of the compound IDs for the compounds to be returned (a required argument)
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "compounds"
})
public class GetCompoundsParams {

    @JsonProperty("compounds")
    private List<String> compounds;
    private Map<java.lang.String, Object> additionalProperties = new HashMap<java.lang.String, Object>();

    @JsonProperty("compounds")
    public List<String> getCompounds() {
        return compounds;
    }

    @JsonProperty("compounds")
    public void setCompounds(List<String> compounds) {
        this.compounds = compounds;
    }

    public GetCompoundsParams withCompounds(List<String> compounds) {
        this.compounds = compounds;
        return this;
    }

    @JsonAnyGetter
    public Map<java.lang.String, Object> getAdditionalProperties() {
        return this.additionalProperties;
    }

    @JsonAnySetter
    public void setAdditionalProperties(java.lang.String name, Object value) {
        this.additionalProperties.put(name, value);
    }

    @Override
    public java.lang.String toString() {
        return ((((("GetCompoundsParams"+" [compounds=")+ compounds)+", additionalProperties=")+ additionalProperties)+"]");
    }

}

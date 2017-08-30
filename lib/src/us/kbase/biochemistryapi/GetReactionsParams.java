
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
 * <p>Original spec-file type: get_reactions_params</p>
 * <pre>
 * Input parameters for the "get_reactions" function.
 *                 list<reaction_id> reactions - a list of the reaction IDs for the reactions to be returned (a required argument)
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "reactions"
})
public class GetReactionsParams {

    @JsonProperty("reactions")
    private List<String> reactions;
    private Map<java.lang.String, Object> additionalProperties = new HashMap<java.lang.String, Object>();

    @JsonProperty("reactions")
    public List<String> getReactions() {
        return reactions;
    }

    @JsonProperty("reactions")
    public void setReactions(List<String> reactions) {
        this.reactions = reactions;
    }

    public GetReactionsParams withReactions(List<String> reactions) {
        this.reactions = reactions;
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
        return ((((("GetReactionsParams"+" [reactions=")+ reactions)+", additionalProperties=")+ additionalProperties)+"]");
    }

}

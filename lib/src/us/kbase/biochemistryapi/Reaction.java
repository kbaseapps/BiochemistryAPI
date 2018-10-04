
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
 * <p>Original spec-file type: Reaction</p>
 * <pre>
 * Data structures for reactions
 *                 reaction_id id - ID of reaction
 *                 string name - primary name of reaction
 *                 string abbrev - abbreviated name of reaction
 *                 list<string> enzymes - list of EC numbers for reaction
 *                 string direction - directionality of reaction
 *                 string reversibility - reversibility of reaction
 *                 float deltaG - estimated delta G of reaction
 *                 float deltaGErr - uncertainty in estimated delta G of reaction
 *                 string equation - reaction equation in terms of compound IDs
 *                 string definition - reaction equation in terms of compound names
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "id",
    "name",
    "abbrev",
    "enzymes",
    "direction",
    "reversibility",
    "deltaG",
    "deltaGErr",
    "equation",
    "definition"
})
public class Reaction {

    @JsonProperty("id")
    private java.lang.String id;
    @JsonProperty("name")
    private java.lang.String name;
    @JsonProperty("abbrev")
    private java.lang.String abbrev;
    @JsonProperty("enzymes")
    private List<String> enzymes;
    @JsonProperty("direction")
    private java.lang.String direction;
    @JsonProperty("reversibility")
    private java.lang.String reversibility;
    @JsonProperty("deltaG")
    private Double deltaG;
    @JsonProperty("deltaGErr")
    private Double deltaGErr;
    @JsonProperty("equation")
    private java.lang.String equation;
    @JsonProperty("definition")
    private java.lang.String definition;
    private Map<java.lang.String, Object> additionalProperties = new HashMap<java.lang.String, Object>();

    @JsonProperty("id")
    public java.lang.String getId() {
        return id;
    }

    @JsonProperty("id")
    public void setId(java.lang.String id) {
        this.id = id;
    }

    public Reaction withId(java.lang.String id) {
        this.id = id;
        return this;
    }

    @JsonProperty("name")
    public java.lang.String getName() {
        return name;
    }

    @JsonProperty("name")
    public void setName(java.lang.String name) {
        this.name = name;
    }

    public Reaction withName(java.lang.String name) {
        this.name = name;
        return this;
    }

    @JsonProperty("abbrev")
    public java.lang.String getAbbrev() {
        return abbrev;
    }

    @JsonProperty("abbrev")
    public void setAbbrev(java.lang.String abbrev) {
        this.abbrev = abbrev;
    }

    public Reaction withAbbrev(java.lang.String abbrev) {
        this.abbrev = abbrev;
        return this;
    }

    @JsonProperty("enzymes")
    public List<String> getEnzymes() {
        return enzymes;
    }

    @JsonProperty("enzymes")
    public void setEnzymes(List<String> enzymes) {
        this.enzymes = enzymes;
    }

    public Reaction withEnzymes(List<String> enzymes) {
        this.enzymes = enzymes;
        return this;
    }

    @JsonProperty("direction")
    public java.lang.String getDirection() {
        return direction;
    }

    @JsonProperty("direction")
    public void setDirection(java.lang.String direction) {
        this.direction = direction;
    }

    public Reaction withDirection(java.lang.String direction) {
        this.direction = direction;
        return this;
    }

    @JsonProperty("reversibility")
    public java.lang.String getReversibility() {
        return reversibility;
    }

    @JsonProperty("reversibility")
    public void setReversibility(java.lang.String reversibility) {
        this.reversibility = reversibility;
    }

    public Reaction withReversibility(java.lang.String reversibility) {
        this.reversibility = reversibility;
        return this;
    }

    @JsonProperty("deltaG")
    public Double getDeltaG() {
        return deltaG;
    }

    @JsonProperty("deltaG")
    public void setDeltaG(Double deltaG) {
        this.deltaG = deltaG;
    }

    public Reaction withDeltaG(Double deltaG) {
        this.deltaG = deltaG;
        return this;
    }

    @JsonProperty("deltaGErr")
    public Double getDeltaGErr() {
        return deltaGErr;
    }

    @JsonProperty("deltaGErr")
    public void setDeltaGErr(Double deltaGErr) {
        this.deltaGErr = deltaGErr;
    }

    public Reaction withDeltaGErr(Double deltaGErr) {
        this.deltaGErr = deltaGErr;
        return this;
    }

    @JsonProperty("equation")
    public java.lang.String getEquation() {
        return equation;
    }

    @JsonProperty("equation")
    public void setEquation(java.lang.String equation) {
        this.equation = equation;
    }

    public Reaction withEquation(java.lang.String equation) {
        this.equation = equation;
        return this;
    }

    @JsonProperty("definition")
    public java.lang.String getDefinition() {
        return definition;
    }

    @JsonProperty("definition")
    public void setDefinition(java.lang.String definition) {
        this.definition = definition;
    }

    public Reaction withDefinition(java.lang.String definition) {
        this.definition = definition;
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
        return ((((((((((((((((((((((("Reaction"+" [id=")+ id)+", name=")+ name)+", abbrev=")+ abbrev)+", enzymes=")+ enzymes)+", direction=")+ direction)+", reversibility=")+ reversibility)+", deltaG=")+ deltaG)+", deltaGErr=")+ deltaGErr)+", equation=")+ equation)+", definition=")+ definition)+", additionalProperties=")+ additionalProperties)+"]");
    }

}

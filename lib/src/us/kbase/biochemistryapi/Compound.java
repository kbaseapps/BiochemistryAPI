
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
 * <p>Original spec-file type: Compound</p>
 * <pre>
 * Data structures for media formulation
 *                 compound_id id - ID of compound
 *                 string abbrev - abbreviated name of compound
 *                 string name - primary name of compound
 *                 list<string> aliases - list of aliases for compound
 *                 float charge - molecular charge of compound
 *                 float deltaG - estimated compound delta G
 *                 float deltaGErr - uncertainty in estimated compound delta G
 *                 string formula - molecular formula of compound
 * </pre>
 * 
 */
@JsonInclude(JsonInclude.Include.NON_NULL)
@Generated("com.googlecode.jsonschema2pojo")
@JsonPropertyOrder({
    "id",
    "abbrev",
    "name",
    "aliases",
    "charge",
    "deltaG",
    "deltaGErr",
    "formula"
})
public class Compound {

    @JsonProperty("id")
    private java.lang.String id;
    @JsonProperty("abbrev")
    private java.lang.String abbrev;
    @JsonProperty("name")
    private java.lang.String name;
    @JsonProperty("aliases")
    private List<String> aliases;
    @JsonProperty("charge")
    private Double charge;
    @JsonProperty("deltaG")
    private Double deltaG;
    @JsonProperty("deltaGErr")
    private Double deltaGErr;
    @JsonProperty("formula")
    private java.lang.String formula;
    private Map<java.lang.String, Object> additionalProperties = new HashMap<java.lang.String, Object>();

    @JsonProperty("id")
    public java.lang.String getId() {
        return id;
    }

    @JsonProperty("id")
    public void setId(java.lang.String id) {
        this.id = id;
    }

    public Compound withId(java.lang.String id) {
        this.id = id;
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

    public Compound withAbbrev(java.lang.String abbrev) {
        this.abbrev = abbrev;
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

    public Compound withName(java.lang.String name) {
        this.name = name;
        return this;
    }

    @JsonProperty("aliases")
    public List<String> getAliases() {
        return aliases;
    }

    @JsonProperty("aliases")
    public void setAliases(List<String> aliases) {
        this.aliases = aliases;
    }

    public Compound withAliases(List<String> aliases) {
        this.aliases = aliases;
        return this;
    }

    @JsonProperty("charge")
    public Double getCharge() {
        return charge;
    }

    @JsonProperty("charge")
    public void setCharge(Double charge) {
        this.charge = charge;
    }

    public Compound withCharge(Double charge) {
        this.charge = charge;
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

    public Compound withDeltaG(Double deltaG) {
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

    public Compound withDeltaGErr(Double deltaGErr) {
        this.deltaGErr = deltaGErr;
        return this;
    }

    @JsonProperty("formula")
    public java.lang.String getFormula() {
        return formula;
    }

    @JsonProperty("formula")
    public void setFormula(java.lang.String formula) {
        this.formula = formula;
    }

    public Compound withFormula(java.lang.String formula) {
        this.formula = formula;
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
        return ((((((((((((((((((("Compound"+" [id=")+ id)+", abbrev=")+ abbrev)+", name=")+ name)+", aliases=")+ aliases)+", charge=")+ charge)+", deltaG=")+ deltaG)+", deltaGErr=")+ deltaGErr)+", formula=")+ formula)+", additionalProperties=")+ additionalProperties)+"]");
    }

}

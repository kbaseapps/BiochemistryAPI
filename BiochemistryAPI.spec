/*
A KBase module: BiochemistryAPI
*/

module BiochemistryAPI {
    /* An identifier for compounds in the KBase biochemistry database. e.g. cpd00001 */
	typedef string compound_id;

	/* A string identifier used for a reaction in a KBase biochemistry. */
    typedef string reaction_id;

    /* Data structures for media formulation

		compound_id id - ID of compound
		string abbrev - abbreviated name of compound
		string name - primary name of compound
		list<string> aliases - list of aliases for compound
		float charge - molecular charge of compound
		float deltaG - estimated compound delta G
		float deltaGErr - uncertainty in estimated compound delta G
		string formula - molecular formula of compound

	*/
    typedef structure {
		compound_id id;
		string abbrev;
		string name;
		list<string> aliases;
		float charge;
		float deltaG;
		float deltaGErr;
		string formula;
    } Compound;

    /* Data structures for media formulation

		reaction_id id - ID of reaction
		string name - primary name of reaction
		string abbrev - abbreviated name of reaction
		list<string> enzymes - list of EC numbers for reaction
		string direction - directionality of reaction
		string reversibility - reversibility of reaction
		float deltaG - estimated delta G of reaction
		float deltaGErr - uncertainty in estimated delta G of reaction
		string equation - reaction equation in terms of compound IDs
		string definition - reaction equation in terms of compound names

	*/
    typedef structure {
		reaction_id id;
		string name;
		string abbrev;
		list<string> enzymes;
		string direction;
		string reversibility;
		float deltaG;
		float deltaGErr;
		string equation;
		string definition;
    } Reaction;
    /*
        This module serves biochemistry content
    */
    /* Input parameters for the "get_reactions" function.

		list<reaction_id> reactions - a list of the reaction IDs for the reactions to be returned (a required argument)
	*/

    typedef structure {
		list<reaction_id> reactions;
    } get_reactions_params;
    /*
    	Returns data for the requested reactions
    */
    funcdef get_reactions(get_reactions_params params) returns (list<Reaction> out_reactions) authentication required;

	/*
	    Input parameters for the "get_compounds" function.
		list<compound_id> compounds - a list of the compound IDs for the compounds to be returned (a required argument)
	*/
    typedef structure {
		list<compound_id> compounds;
    } get_compounds_params;
    /*
    	Returns data for the requested compounds
    */
    funcdef get_compounds(get_compounds_params params) returns (list<Compound> out_compounds) authentication required;

    typedef structure {
		list<string> compound_structures;
    } depict_compounds_params;
    /*
    	Returns a list of depictions for the compound_structures in SVG format
    */
    funcdef depict_compounds(depict_compounds_params params) returns (list<string> depictions) authentication required;


};

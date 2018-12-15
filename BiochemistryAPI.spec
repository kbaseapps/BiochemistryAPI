/*
A KBase module: BiochemistryAPI
*/

module BiochemistryAPI {
    /* @range(0,1)*/
    typedef int bool;

    /* An identifier for compounds in the KBase biochemistry database. e.g. cpd00001 */
	typedef string compound_id;

	/* A string identifier used for a reaction in a KBase biochemistry. */
    typedef string reaction_id;

    /* Data structures for compounds

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

    /* Data structures for reactions

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
        A molecule structure in InChI or SMILES format
    */

    typedef string mol_structure;

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
    funcdef get_compounds(get_compounds_params params) returns (list<Compound> out_compounds);

    /*
	    Input parameters for the "search_compounds" function.
		string query - a query string to match against names & aliases
		int limit - maximum number of results to return, defaults to 10
	*/
    typedef structure {
		string query;
		int limit;
    } search_compounds_params;
    /*
    	Returns compounds which match a string
    */
    funcdef search_compounds(search_compounds_params params) returns (list<Compound> out_compounds);

    /*
	    Input parameters for the "search_reactions" function.
		string query - a query string to match against names & aliases
		int limit - maximum number of results to return, defaults to 10
	*/
    typedef structure {
		string query;
		int limit;
    } search_reactions_params;
    /*
    	Returns reactions which match a string
    */
    funcdef search_reactions(search_reactions_params params) returns (list<Reaction> out_reactions);

    typedef structure {
		mol_structure query;
    } substructure_search_params;
    /*
    	Returns compound ids for compounds that contain the query substructure
    */
    funcdef substructure_search(substructure_search_params params) returns (list<compound_id> matching_ids);

    /*
        mol_structure query: Either InChI or SMILES string
		string fp_type: Either MACCS or Morgan fingerprints
		float min_similarity: In range 0-1
    */

    typedef structure {
		mol_structure query;
		string fp_type;
		float min_similarity;
    } similarity_search_params;
    /*
    	Returns compound ids for compounds that have greater fingerprint similarity than the min_similarity threshold
    */
    funcdef similarity_search(similarity_search_params params) returns (list<compound_id> matching_ids);


    typedef structure {
		list<mol_structure> structures;
    } depict_compounds_params;
    /*
    	Returns a list of depictions for the compound_structures in SVG format
    */
    funcdef depict_compounds(depict_compounds_params params) returns (list<string> depictions);

    typedef structure {
		list<mol_structure> structures;
		bool optimize;
		string output;
    } calculate_3D_coords_params;
    /*
    	Returns molecules with 3D coordinates.

    	list<mol_structure> compound_structures: compounds in InChI or SMILES
		bool optimize: should forcefeild optimization be run?
		string output: The outpuf format, one of 'mol' or 'pdb'
    */
    funcdef calculate_3D_coords(calculate_3D_coords_params params) returns (list<string> output);
};

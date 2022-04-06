# -*- coding: utf-8 -*-
# BEGIN_HEADER
import logging

from BiochemistryAPI import utils

logging.basicConfig(level=logging.INFO)
# END_HEADER


class BiochemistryAPI:
    """
    Module Name:
    BiochemistryAPI

    Module Description:
    A KBase module: BiochemistryAPI
    """

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.3.0"
    GIT_URL = "https://github.com/kbaseapps/BiochemistryAPI.git"
    GIT_COMMIT_HASH = "410765239da8490da9476411476beb7f2d60140d"

    # BEGIN_CLASS_HEADER

    # END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        # BEGIN_CONSTRUCTOR
        self.config = config
        self.scratch = config["scratch"]
        data_dir = "/kb/module/data/"
        self.compounds, self.comp_search_dict = utils.dict_from_file(
            data_dir + "compounds.tsv"
        )

        self.reactions, self.rxn_search_dict = utils.dict_from_file(
            data_dir + "reactions.tsv"
        )
        self.comp_aliases = utils.alias_dict_from_file(
            data_dir + "Compounds_Aliases.tsv"
        )
        self.rxn_aliases = utils.alias_dict_from_file(
            data_dir + "Reactions_Aliases.tsv"
        )
        self.ec_classes = utils.alias_dict_from_file(
            data_dir + "Enzyme_Class_Reactions_Aliases.tsv"
        )

        logging.info(
            "Loaded {} compounds and {} reactions".format(
                len(self.compounds), len(self.reactions)
            )
        )
        self.structures = utils.make_mol_tuples(self.compounds.values())
        logging.info("Cached compound structures")
        # END_CONSTRUCTOR
        pass

    def get_reactions(self, ctx, params):
        """
        Returns data for the requested reactions
        :param params: instance of type "get_reactions_params" (Input
           parameters for the "get_reactions" function. list<reaction_id>
           reactions - a list of the reaction IDs for the reactions to be
           returned (a required argument)) -> structure: parameter
           "reactions" of list of type "reaction_id" (A string identifier
           used for a reaction in a KBase biochemistry.)
        :returns: instance of list of type "Reaction" (Data structures for
           reactions reaction_id id - ID of reaction string name - primary
           name of reaction string abbrev - abbreviated name of reaction
           list<string> enzymes - list of EC numbers for reaction string
           direction - directionality of reaction string reversibility -
           reversibility of reaction float deltaG - estimated delta G of
           reaction float deltaGErr - uncertainty in estimated delta G of
           reaction string equation - reaction equation in terms of compound
           IDs string definition - reaction equation in terms of compound
           names) -> structure: parameter "id" of type "reaction_id" (A
           string identifier used for a reaction in a KBase biochemistry.),
           parameter "name" of String, parameter "abbrev" of String,
           parameter "enzymes" of list of String, parameter "direction" of
           String, parameter "reversibility" of String, parameter "deltaG" of
           Double, parameter "deltaGErr" of Double, parameter "equation" of
           String, parameter "definition" of String
        """
        # ctx is the context object
        # return variables are: out_reactions
        # BEGIN get_reactions
        logging.info("Starting get_reactions")
        logging.info("Params: {}".format(params))
        utils.check_param(params, ["reactions"])
        out_reactions = []
        for x in params["reactions"]:
            id = x.split("/")[-1]
            rxn = self.reactions.get(id, None)
            if rxn:
                rxn["aliases"] = self.rxn_aliases.get(id, "")
                rxn["enzymes"] = self.ec_classes.get(id, "")
            out_reactions.append(rxn)
        # END get_reactions

        # At some point might do deeper type checking...
        if not isinstance(out_reactions, list):
            raise ValueError(
                "Method get_reactions return value "
                + "out_reactions is not type list as required."
            )
        # return the results
        return [out_reactions]

    def get_compounds(self, ctx, params):
        """
        Returns data for the requested compounds
        :param params: instance of type "get_compounds_params" (Input
           parameters for the "get_compounds" function. list<compound_id>
           compounds - a list of the compound IDs for the compounds to be
           returned (a required argument)) -> structure: parameter
           "compounds" of list of type "compound_id" (An identifier for
           compounds in the KBase biochemistry database. e.g. cpd00001)
        :returns: instance of list of type "Compound" (Data structures for
           compounds compound_id id - ID of compound string abbrev -
           abbreviated name of compound string name - primary name of
           compound list<string> aliases - list of aliases for compound float
           charge - molecular charge of compound float deltaG - estimated
           compound delta G float deltaGErr - uncertainty in estimated
           compound delta G string formula - molecular formula of compound)
           -> structure: parameter "id" of type "compound_id" (An identifier
           for compounds in the KBase biochemistry database. e.g. cpd00001),
           parameter "abbrev" of String, parameter "name" of String,
           parameter "aliases" of list of String, parameter "charge" of
           Double, parameter "deltaG" of Double, parameter "deltaGErr" of
           Double, parameter "formula" of String
        """
        # ctx is the context object
        # return variables are: out_compounds
        # BEGIN get_compounds
        logging.info("Starting get_compounds")
        logging.info("Params: {}".format(params))
        utils.check_param(params, ["compounds"])
        out_compounds = []
        for x in params["compounds"]:
            id = x.split("/")[-1]
            comp = self.compounds.get(id, None)
            if comp:
                comp["aliases"] = self.comp_aliases.get(id, "")
            out_compounds.append(comp)
        # END get_compounds

        # At some point might do deeper type checking...
        if not isinstance(out_compounds, list):
            raise ValueError(
                "Method get_compounds return value "
                + "out_compounds is not type list as required."
            )
        # return the results
        return [out_compounds]

    def search_compounds(self, ctx, params):
        """
        Returns compounds which match a string
        :param params: instance of type "search_compounds_params" (Input
           parameters for the "search_compounds" function. string query - a
           query string to match against names & aliases int limit - maximum
           number of results to return, defaults to 10) -> structure:
           parameter "query" of String, parameter "limit" of Long
        :returns: instance of list of type "Compound" (Data structures for
           compounds compound_id id - ID of compound string abbrev -
           abbreviated name of compound string name - primary name of
           compound list<string> aliases - list of aliases for compound float
           charge - molecular charge of compound float deltaG - estimated
           compound delta G float deltaGErr - uncertainty in estimated
           compound delta G string formula - molecular formula of compound)
           -> structure: parameter "id" of type "compound_id" (An identifier
           for compounds in the KBase biochemistry database. e.g. cpd00001),
           parameter "abbrev" of String, parameter "name" of String,
           parameter "aliases" of list of String, parameter "charge" of
           Double, parameter "deltaG" of Double, parameter "deltaGErr" of
           Double, parameter "formula" of String
        """
        # ctx is the context object
        # return variables are: out_compounds
        # BEGIN search_compounds
        logging.info("Starting search_compounds")
        logging.info("Params: {}".format(params))
        utils.check_param(params, ["query"], ["limit"])
        limit = params.get("limit", 10)
        normed_query = params["query"].strip().lower().translate(utils.ttable)
        matching_ids = self.comp_search_dict.get(normed_query, [])
        if len(matching_ids) > limit:
            matching_ids = matching_ids[:limit]
        out_compounds = [self.compounds[cid] for cid in matching_ids]
        # END search_compounds

        # At some point might do deeper type checking...
        if not isinstance(out_compounds, list):
            raise ValueError(
                "Method search_compounds return value "
                + "out_compounds is not type list as required."
            )
        # return the results
        return [out_compounds]

    def search_reactions(self, ctx, params):
        """
        Returns reactions which match a string
        :param params: instance of type "search_reactions_params" (Input
           parameters for the "search_reactions" function. string query - a
           query string to match against names & aliases int limit - maximum
           number of results to return, defaults to 10) -> structure:
           parameter "query" of String, parameter "limit" of Long
        :returns: instance of list of type "Reaction" (Data structures for
           reactions reaction_id id - ID of reaction string name - primary
           name of reaction string abbrev - abbreviated name of reaction
           list<string> enzymes - list of EC numbers for reaction string
           direction - directionality of reaction string reversibility -
           reversibility of reaction float deltaG - estimated delta G of
           reaction float deltaGErr - uncertainty in estimated delta G of
           reaction string equation - reaction equation in terms of compound
           IDs string definition - reaction equation in terms of compound
           names) -> structure: parameter "id" of type "reaction_id" (A
           string identifier used for a reaction in a KBase biochemistry.),
           parameter "name" of String, parameter "abbrev" of String,
           parameter "enzymes" of list of String, parameter "direction" of
           String, parameter "reversibility" of String, parameter "deltaG" of
           Double, parameter "deltaGErr" of Double, parameter "equation" of
           String, parameter "definition" of String
        """
        # ctx is the context object
        # return variables are: out_reactions
        # BEGIN search_reactions
        logging.info("Starting search_reactions")
        logging.info("Params: {}".format(params))
        utils.check_param(params, ["query"], ["limit"])
        limit = params.get("limit", 10)
        normed_query = params["query"].strip().lower().translate(utils.ttable)
        matching_ids = self.rxn_search_dict.get(normed_query, [])
        if len(matching_ids) > limit:
            matching_ids = matching_ids[:limit]
        out_reactions = [self.reactions[cid] for cid in matching_ids]
        # END search_reactions

        # At some point might do deeper type checking...
        if not isinstance(out_reactions, list):
            raise ValueError(
                "Method search_reactions return value "
                + "out_reactions is not type list as required."
            )
        # return the results
        return [out_reactions]

    def substructure_search(self, ctx, params):
        """
        Returns compound ids for compounds that contain the query substructure
        :param params: instance of type "substructure_search_params" ->
           structure: parameter "query" of type "mol_structure" (A molecule
           structure in InChI or SMILES format)
        :returns: instance of list of type "compound_id" (An identifier for
           compounds in the KBase biochemistry database. e.g. cpd00001)
        """
        # ctx is the context object
        # return variables are: matching_ids
        # BEGIN substructure_search
        logging.info("Starting substructure_search")
        logging.info("Params: {}".format(params))
        utils.check_param(params, ["query"])
        matching_ids = utils.substructure_search(params["query"], self.structures)
        # END substructure_search

        # At some point might do deeper type checking...
        if not isinstance(matching_ids, list):
            raise ValueError(
                "Method substructure_search return value "
                + "matching_ids is not type list as required."
            )
        # return the results
        return [matching_ids]

    def similarity_search(self, ctx, params):
        """
        Returns compound ids for compounds that have greater fingerprint similarity than the min_similarity threshold
        :param params: instance of type "similarity_search_params"
           (mol_structure query: Either InChI or SMILES string string
           fp_type: Either MACCS or Morgan fingerprints float min_similarity:
           In range 0-1) -> structure: parameter "query" of type
           "mol_structure" (A molecule structure in InChI or SMILES format),
           parameter "fp_type" of String, parameter "min_similarity" of Double
        :returns: instance of list of type "compound_id" (An identifier for
           compounds in the KBase biochemistry database. e.g. cpd00001)
        """
        # ctx is the context object
        # return variables are: matching_ids
        # BEGIN similarity_search
        logging.info("Starting similarity_search")
        logging.info("Params: {}".format(params))
        utils.check_param(params, ["query"], ["fp_type", "min_similarity"])
        params["structures"] = self.structures
        matching_ids = utils.similarity_search(**params)
        # END similarity_search

        # At some point might do deeper type checking...
        if not isinstance(matching_ids, list):
            raise ValueError(
                "Method similarity_search return value "
                + "matching_ids is not type list as required."
            )
        # return the results
        return [matching_ids]

    def depict_compounds(self, ctx, params):
        """
        Returns a list of depictions for the compound_structures in SVG format
        :param params: instance of type "depict_compounds_params" ->
           structure: parameter "structures" of list of type "mol_structure"
           (A molecule structure in InChI or SMILES format)
        :returns: instance of list of String
        """
        # ctx is the context object
        # return variables are: depictions
        # BEGIN depict_compounds
        logging.info("Starting depict_compounds")
        logging.info("Params: {}".format(params))
        utils.check_param(params, ["structures"])
        depictions = [utils.depict_compound(struct) for struct in params["structures"]]
        # END depict_compounds

        # At some point might do deeper type checking...
        if not isinstance(depictions, list):
            raise ValueError(
                "Method depict_compounds return value "
                + "depictions is not type list as required."
            )
        # return the results
        return [depictions]

    def calculate_3D_coords(self, ctx, params):
        """
        Returns molecules with 3D coordinates.
        list<mol_structure> compound_structures: compounds in InChI or SMILES
            bool optimize: should forcefeild optimization be run?
            string output: The outpuf format, one of 'mol' or 'pdb'
        :param params: instance of type "calculate_3D_coords_params" ->
           structure: parameter "structures" of list of type "mol_structure"
           (A molecule structure in InChI or SMILES format), parameter
           "optimize" of type "bool" (@range(0,1)), parameter "output" of
           String
        :returns: instance of list of String
        """
        # ctx is the context object
        # return variables are: output
        # BEGIN calculate_3D_coords
        logging.info("Starting calculate_3D_coords")
        logging.info("Params: {}".format(params))
        utils.check_param(params, ["structures"])
        output = [
            utils.get_3d_mol(
                struct,
                params.get("output", "mol"),
                params.get("optimize"),
            )
            for struct in params["structures"]
        ]
        # END calculate_3D_coords

        # At some point might do deeper type checking...
        if not isinstance(output, list):
            raise ValueError(
                "Method calculate_3D_coords return value "
                + "output is not type list as required."
            )
        # return the results
        return [output]

    def status(self, ctx):
        # BEGIN_STATUS
        returnVal = {
            "state": "OK",
            "message": "",
            "version": self.VERSION,
            "git_url": self.GIT_URL,
            "git_commit_hash": self.GIT_COMMIT_HASH,
        }
        # END_STATUS
        return [returnVal]

# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
import csv
#END_HEADER


class BiochemistryAPI:
    '''
    Module Name:
    BiochemistryAPI

    Module Description:
    A KBase module: BiochemistryAPI
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    @staticmethod
    def _check_param(in_params, req_param, opt_param=list()):
        """
        Check if each of the params in the list are in the input params
        """
        for param in req_param:
            if param not in in_params:
                raise ValueError('{} parameter is required'.format(param))
        defined_param = set(req_param + opt_param)
        for param in in_params:
            if param not in defined_param:
                print(
                    "WARNING: received unexpected parameter {}".format(param))

    def dict_from_file(self, path, key='id', dialect='excel-tab'):
        """
        Build a dictionary from an object array in a file
        :param path: local path to object
        :param key: what field should be used as the key
        :param dialect: excel-tab for TSV or excel for CSV
        :return:
        """
        if not os.path.exists(path):
            raise ValueError("File not found: {}".format(path))
        reader = csv.DictReader(open(path), dialect=dialect)
        return dict([(x[key], x) for x in reader])

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.scratch = config['scratch']
        self.compounds = self.dict_from_file("/kb/module/data/compounds.tsv")
        self.reactions = self.dict_from_file("/kb/module/data/reactions.tsv")
        print("Loaded {} compounds and {} reactions".format(
            len(self.compounds), len(self.reactions)))
        #END_CONSTRUCTOR
        pass


    def get_reactions(self, ctx, input):
        """
        Returns data for the requested reactions
        :param input: instance of type "get_reactions_params" (Input
           parameters for the "get_reactions" function. list<reaction_id>
           reactions - a list of the reaction IDs for the reactions to be
           returned (a required argument)) -> structure: parameter
           "reactions" of list of type "reaction_id" (A string identifier
           used for a reaction in a KBase biochemistry.)
        :returns: instance of list of type "Reaction" (Data structures for
           media formulation reaction_id id - ID of reaction string name -
           primary name of reaction string abbrev - abbreviated name of
           reaction list<string> enzymes - list of EC numbers for reaction
           string direction - directionality of reaction string reversibility
           - reversibility of reaction float deltaG - estimated delta G of
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
        #BEGIN get_reactions
        self._check_param(input, ['reactions'])
        out_reactions = [self.reactions.get(x.split('/')[-1], None) for x in
                         input['reactions']]
        #END get_reactions

        # At some point might do deeper type checking...
        if not isinstance(out_reactions, list):
            raise ValueError('Method get_reactions return value ' +
                             'out_reactions is not type list as required.')
        # return the results
        return [out_reactions]

    def get_compounds(self, ctx, input):
        """
        Returns data for the requested compounds
        :param input: instance of type "get_compounds_params" (Input
           parameters for the "get_compounds" function. list<compound_id>
           compounds - a list of the compound IDs for the compounds to be
           returned (a required argument)) -> structure: parameter
           "compounds" of list of type "compound_id" (An identifier for
           compounds in the KBase biochemistry database. e.g. cpd00001)
        :returns: instance of list of type "Compound" (Data structures for
           media formulation compound_id id - ID of compound string abbrev -
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
        #BEGIN get_compounds
        self._check_param(input, ['compounds'])
        out_compounds = [self.compounds.get(x.split('/')[-1]) for x in
                         input['compounds']]
        #END get_compounds

        # At some point might do deeper type checking...
        if not isinstance(out_compounds, list):
            raise ValueError('Method get_compounds return value ' +
                             'out_compounds is not type list as required.')
        # return the results
        return [out_compounds]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]

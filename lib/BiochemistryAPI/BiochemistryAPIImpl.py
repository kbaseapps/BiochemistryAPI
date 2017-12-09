# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os
import csv
from BiochemistryAPI.utils import depict_compound
from collections import defaultdict
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
    VERSION = "0.1.1"
    GIT_URL = "git@github.com:kbaseapps/BiochemistryAPI.git"
    GIT_COMMIT_HASH = "b82efb604f55537e62cdf2bc8960369a8c640d27"

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

    @staticmethod
    def dict_from_file(path, key='id', dialect='excel-tab'):
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

    @staticmethod
    def alias_dict_from_file(path, dialect='excel-tab'):
        """
        Build a dictionary from an object array in a file
        :param path: local path to object
        :param dialect: excel-tab for TSV or excel for CSV
        :return:
        """
        alias_mappings = defaultdict(list)
        with open(path) as infile:
            r = csv.DictReader(infile, dialect=dialect)
            for line in r:
                for seed_id in line['MS ID'].split('|'):
                    if line['Source'] == 'Enzyme Class':
                        alias_mappings[seed_id].append(line['External ID'])
                    else:
                        alias_mappings[seed_id].append('%s:%s' % (
                            line['Source'].strip(), line['External ID']))
        return alias_mappings

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.scratch = config['scratch']
        data_dir = '/kb/module/data/'
        self.compounds = self.dict_from_file(data_dir + "compounds.tsv")
        self.reactions = self.dict_from_file(data_dir + "reactions.tsv")
        self.comp_aliases = self.alias_dict_from_file(data_dir +
                                                      "Compounds_Aliases.tsv")
        self.rxn_aliases = self.alias_dict_from_file(data_dir +
                                                      "Reactions_Aliases.tsv")
        self.ec_classes = self.alias_dict_from_file(data_dir + 'Enzyme_Class_Reactions_Aliases.tsv')

        print("Loaded {} compounds and {} reactions".format(
            len(self.compounds), len(self.reactions)))
        #END_CONSTRUCTOR
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
        self._check_param(params, ['reactions'])
        out_reactions = []
        for x in params['reactions']:
            id = x.split('/')[-1]
            rxn = self.reactions.get(id, None)
            if rxn:
                rxn['aliases'] = self.rxn_aliases.get(id, '')
                rxn['enzymes'] = self.ec_classes.get(id, '')
            out_reactions.append(rxn)
        #END get_reactions

        # At some point might do deeper type checking...
        if not isinstance(out_reactions, list):
            raise ValueError('Method get_reactions return value ' +
                             'out_reactions is not type list as required.')
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
        self._check_param(params, ['compounds'])
        out_compounds = []
        for x in params['compounds']:
            id = x.split('/')[-1]
            comp = self.compounds.get(id, None)
            if comp:
                comp['aliases'] = self.comp_aliases.get(id, '')
            out_compounds.append(comp)
        #END get_compounds

        # At some point might do deeper type checking...
        if not isinstance(out_compounds, list):
            raise ValueError('Method get_compounds return value ' +
                             'out_compounds is not type list as required.')
        # return the results
        return [out_compounds]

    def depict_compounds(self, ctx, params):
        """
        Returns a list of depictions for the compound_structures in SVG format
        :param params: instance of type "depict_compounds_params" ->
           structure: parameter "compound_structures" of list of String
        :returns: instance of list of String
        """
        # ctx is the context object
        # return variables are: depictions
        #BEGIN depict_compounds
        self._check_param(params, ['structures'])
        depictions = [depict_compound(struct)
                      for struct in params['structures']]
        #END depict_compounds

        # At some point might do deeper type checking...
        if not isinstance(depictions, list):
            raise ValueError('Method depict_compounds return value ' +
                             'depictions is not type list as required.')
        # return the results
        return [depictions]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]

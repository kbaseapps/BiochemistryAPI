# -*- coding: utf-8 -*-

try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

import os
import time
import unittest

from biokbase.workspace.client import Workspace as workspaceService
from BiochemistryAPI.BiochemistryAPIImpl import BiochemistryAPI
from BiochemistryAPI.BiochemistryAPIServer import MethodContext
from BiochemistryAPI.authclient import KBaseAuth as _KBaseAuth


class BiochemistryAPITest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('BiochemistryAPI'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'BiochemistryAPI',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL)
        cls.serviceImpl = BiochemistryAPI(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_BiochemistryAPI_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})  # noqa
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    def test_bad_input(self):
        with self.assertRaisesRegexp(ValueError, 'parameter is required'):
            self.getImpl().get_compounds(self.ctx, {})
        with self.assertRaisesRegexp(ValueError, 'parameter is required'):
            self.getImpl().get_reactions(self.ctx, {})
        with self.assertRaisesRegexp(ValueError, 'parameter is required'):
            self.getImpl().substructure_search(self.ctx, {})
        with self.assertRaisesRegexp(ValueError, 'parameter is required'):
            self.getImpl().similarity_search(self.ctx, {})
        with self.assertRaisesRegexp(ValueError, 'parameter is required'):
            self.getImpl().depict_compounds(self.ctx, {})
        with self.assertRaisesRegexp(ValueError, 'Invalid fingerprint type'):
            self.getImpl().similarity_search(self.ctx, {'query': 'C(=O)O', 'fp_type': 'foo'})
        with self.assertRaisesRegexp(ValueError, 'Invalid min_similarity'):
            self.getImpl().similarity_search(self.ctx, {'query': 'C', 'min_similarity': 1.1})
        with self.assertRaisesRegexp(ValueError, 'Invalid min_similarity'):
            self.getImpl().similarity_search(self.ctx, {'query': 'C', 'min_similarity': '0.8'})

    def test_get_compounds(self):
        cpds = self.getImpl().get_compounds(self.ctx, {"compounds":
                                            ["48/1/1/compounds/id/cpd00011",
                                             'cpd00002', "cpd00007"]})[0]
        assert len(cpds) == 3
        assert cpds[0]['id'] == 'cpd00011'
        assert cpds[1]['id'] == 'cpd00002'
        assert cpds[0]['aliases'] == ['BiGG1:co2', 'BiGG1:dco2',
                                      'KEGG:C00011', 'MetaCyc:CARBON-DIOXIDE',
                                      'PlantCyc:CARBON-DIOXIDE', 'BiGG:co2']
        missing_col = {'name', 'formula', 'charge', 'deltaG', 'deltaGErr',
                       'abbreviation', 'aliases'} - set(cpds[0].keys())
        if missing_col:
            raise AssertionError("Missing Columns:", missing_col)

    def test_get_reactions(self):
        rxns = self.getImpl().get_reactions(self.ctx, {"reactions":
                                            ["48/1/1/reactions/id/rxn00011",
                                             'rxn00002', "rxn00007"]})[0]
        assert len(rxns) == 3
        assert rxns[0]['id'] == 'rxn00011'
        assert rxns[1]['id'] == 'rxn00002'
        assert rxns[0]['aliases'] == ['KEGG:R00014', 'MetaCyc:RXN-12583.c',
                                      'PlantCyc:RXN-12583.c']
        assert rxns[0]['enzymes'] == ['1.2.4.1', '2.2.1.6', '4.1.1.1']
        missing_col = {'name', 'direction'} - set(rxns[0].keys())
        if missing_col:
            raise AssertionError("Missing Columns:", missing_col)

    def test_search_compounds(self):
        res = self.getImpl().search_compounds(self.ctx, {'query': 'cpd00007'})[0]
        self.assertEqual(len(res), 1)
        self.assertEqual(res[0]['id'], 'cpd00007')
        res = self.getImpl().search_compounds(self.ctx, {'query': 'O2'})[0]
        self.assertEqual(len(res), 3)
        self.assertEqual(res[0]['id'], 'cpd00007')
        res = self.getImpl().search_compounds(self.ctx, {'query': 'OXYGEN-MOLECULE'})[0]
        self.assertEqual(len(res), 1)
        self.assertEqual(res[0]['id'], 'cpd00007')
        res = self.getImpl().search_compounds(self.ctx, {'query': 'Pyruvate',
                                                         'limit': 5})[0]
        self.assertEqual(len(res), 5)
        self.assertEqual(res[0]['id'], 'cpd00020')

    def test_search_reactions(self):
        res = self.getImpl().search_reactions(self.ctx, {'query': 'rxn00001'})[0]
        self.assertEqual(len(res), 1)
        self.assertEqual(res[0]['id'], 'rxn00001')
        res = self.getImpl().search_reactions(self.ctx, {'query': 'R00004'})[0]
        self.assertEqual(len(res), 1)
        self.assertEqual(res[0]['id'], 'rxn00001')
        res = self.getImpl().search_reactions(self.ctx, {'query': 'diphosphate phosphohydrolase'})[0]
        self.assertEqual(len(res), 1)
        self.assertEqual(res[0]['id'], 'rxn00001')
        res = self.getImpl().search_reactions(self.ctx, {'query': 'phosphohydrolase'})[0]
        self.assertEqual(len(res), 10)
        self.assertEqual(res[0]['id'], 'rxn00001')
        res = self.getImpl().search_reactions(self.ctx, {'query': 'phosphohydrolase',
                                                         'limit': 500})[0]
        self.assertEqual(len(res), 385)
        self.assertEqual(res[0]['id'], 'rxn00001')

    def test_substructure_search(self):
        ids = self.getImpl().substructure_search(self.ctx, {'query': 'C(=O)O'})[0]
        self.assertEqual(len(ids), 7262)
        self.assertEqual(ids[0], 'cpd00017')

    def test_similarity_search(self):
        ids = self.getImpl().similarity_search(
            self.ctx, {'query': 'InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)'})[0]
        self.assertEqual(len(ids), 4)
        self.assertEqual(ids[0], 'cpd00029')

        ids = self.getImpl().similarity_search(self.ctx, {'query': 'C(=O)O', "fp_type": "RDKit",
                                                          "min_similarity": 0.5})[0]
        self.assertEqual(len(ids), 6)
        self.assertEqual(ids[0], 'cpd00047')

    def test_depict_compounds(self):
        svgs = self.getImpl().depict_compounds(
            self.ctx, {'structures': ['C(=O)O', 'InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)', 'foo']})[0]
        self.assertEqual(len(svgs), 3)
        self.assertIn('<svg', svgs[0])
        self.assertIn('<svg', svgs[1])
        self.assertEqual(svgs[2], '')

    def test_calculate_3D_coords(self):
        mols = self.getImpl().calculate_3D_coords(
            self.ctx, {'structures': ['C(=O)O', 'InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)', 'foo']})[0]
        self.assertEqual(len(mols), 3)
        self.assertIn('  -1.1054   -0.2402   -0.0000 O', mols[0])
        self.assertIn('  -0.3889    1.2404    0.0153 O', mols[1])
        self.assertEqual(mols[2], '')

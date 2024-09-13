from __future__ import absolute_import

# import non-local packages
from pandas import DataFrame, read_csv, concat, set_option
from cobra.io import write_sbml_model, read_sbml_model
from configparser import ConfigParser
from numpy import load, nan, ndarray
from zipfile import ZipFile
from typing import Iterable
from math import isclose
from pathlib import Path
import json, time, copy
import sys, os, re
import hashlib
import sigfig
import logging
import uuid

# define a local pointer to import package versions
path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
path2 = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
kbmodules_path = "/scratch/shared/code/KBBaseModules"
modelseed_path = "/scratch/shared/ModelSEEDpy_APF"
modelseed_path = "/scratch/shared/ModelSEEDpy_APF"
mscommunity_path = "/scratch/shared/code/MSCommunity"
utilsModule_path = "/scratch/shared/code/chenry_utility_module/lib"
msrecon_path = "/scratch/shared/code/KB-ModelSEEDReconstruction/lib"
commscores_path = "/scratch/shared/code/CommScores"
for p in [path, path2, kbmodules_path, modelseed_path, mscommunity_path, utilsModule_path, msrecon_path, commscores_path]:
    sys.path.insert(0, p)
# print(sys.path) 

## cobrakbase imports
from cobrakbase.core.kbasefba import FBAModel

## modelseedpy imports
from modelseedpy.core.fbahelper import FBAHelper
# from modelseedpy import MSPackageManager
# from modelseedpy import AnnotationOntology, MSPackageManager, MSMedia, MSModelUtil, MSBuilder, MSATPCorrection, MSGapfill, MSGrowthPhenotype, MSGrowthPhenotypes, ModelSEEDBiochem
# from modelseedpy.core.mstemplate import MSTemplateBuilder
# from modelseedpy.core.msgenome import normalize_role
from modelseedpy.core import MSMinimalMedia
# from modelseedpy.core.msensemble import MSEnsemble
from modelseedpy.helpers import get_template

## mscommunity imports
from mscommunity.commhelper import build_from_species_models

# from installed_clients.WorkspaceClient import Workspace

# # define the standard files
# codebase = os.environ.get("CODE_BASE","/scratch/shared/code")
# config_file = os.environ.get("KBDEVUTIL_CONFIG", codebase+"/sharedconfig.cfg")
# config_parse = ConfigParser()
# config_parse.read(config_file)
# config = {}
# for nameval in config_parse.items("DevEnv"):
#     config[nameval[0]] = nameval[1]
# paths = config["syspaths"].split(";")
# for i,filepath in enumerate(paths):
#     if filepath[0:1] != "/":
#         paths[i] = codebase+"/"+filepath
# sys.path = paths + sys.path


# logger = logging.getLogger(__name__)
# logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
#                             level=logging.INFO)


class BadModels(Exception):
    def __init__(self, message):
        print(message)


# class BaseModule:
#     def __init__(self,name,config,module_dir="/kb/module",working_dir=None,token=None,clients={},callback=None):
#         #Initializing flexible container for client libraries which will be lazy loaded as needed
#         self.version = "0.1.1.bm"
#         self.clients = {}
#         self.callback_url = callback
#         for item in clients:
#             self.clients[item] = clients[item]
#         #Initializing config, name, and token
#         self.config = config
#         self.validate_args(self.config,[],{
#             "max_retry":3,
#             "workspace-url":"https://kbase.us/services/ws",
#         })
#         self.token = token
#         self.name = name
#         self.module_dir = module_dir
#         #Initializing working directory if specified, otherwise using config scratch
#         if working_dir:
#             self.working_dir = working_dir
#         else:
#             self.working_dir = config['scratch']
#         self.reset_attributes()
    
#     #########METHOD CALL INITIALIZATION FUNCTIONS#######################
#     def initialize_call(self,method,params,print_params=False,no_print=[],no_prov_params=[]):
#         if not self.initialized:
#             self.obj_created = []
#             self.input_objects = []
#             self.method = method
#             self.params = copy.deepcopy(params)
#             for item in no_prov_params:
#                 if item in self.params:
#                     del self.params[item]
#             self.initialized = True
#             if "workspace" in params:
#                 self.set_ws(params["workspace"])
#             elif "output_workspace" in params:
#                 self.set_ws(params["output_workspace"])
#             if print_params:
#                 saved_data = {}
#                 for item in no_print:
#                     if item in params:
#                         saved_data[item] = params[item]
#                         del params[item]
#                 logger.info(method+":"+json.dumps(params,indent=4))
#                 for item in saved_data:
#                     params[item] = saved_data[item]
                
#     def reset_attributes(self):
#         #Initializing stores tracking objects created and input objects
#         self.obj_created = []
#         self.input_objects = []
#         #Initializing attributes tracking method data to support provencance and context
#         self.method = None
#         self.params = {}
#         self.initialized = False
#         self.ws_id = None
#         self.ws_name = None
#         #Computing timestamp
#         ts = time.gmtime()
#         self.timestamp = time.strftime("%Y-%m-%d %H:%M:%S", ts)
    
#     #########CLIENT RETRIEVAL AND INITIALIZATION FUNCTIONS#######################
#     def ws_client(self):
#         if "Workspace" not in self.clients:
#             if "devenv" in self.config and self.config["devenv"] == "1":
#                 from kbbasemodules import Workspace
#             else:
#                 from installed_clients.WorkspaceClient import Workspace
#             self.clients["Workspace"] = Workspace(self.config["workspace-url"], token=self.token)
#         return self.clients["Workspace"]
    
#     def report_client(self):
#         if "KBaseReport" not in self.clients:
#             if "devenv" in self.config and self.config["devenv"] == "1":
#                 from kbbasemodules import KBaseReport
#             else:
#                 from installed_clients.KBaseReportClient import KBaseReport
#             self.clients["KBaseReport"] = KBaseReport(self.callback_url,token=self.token)
#         return self.clients["KBaseReport"]
    
#     def dfu_client(self):
#         if "DataFileUtil" not in self.clients:
#             if "devenv" in self.config and self.config["devenv"] == "1":
#                 from kbbasemodules import DataFileUtil
#             else:
#                 from installed_clients.DataFileUtilClient import DataFileUtil
#             self.clients["DataFileUtil"] = DataFileUtil(self.callback_url,token=self.token)
#         return self.clients["DataFileUtil"]
    
#     def gfu_client(self):
#         if "GenomeFileUtil" not in self.clients:
#             if "devenv" in self.config and self.config["devenv"] == "1":
#                 from kbbasemodules import GenomeFileUtil
#             else:
#                 from installed_clients.GenomeFileUtilClient import GenomeFileUtil
#             self.clients["GenomeFileUtil"] = GenomeFileUtil(self.callback_url,token=self.token)
#         return self.clients["GenomeFileUtil"]
    
#     def afu_client(self):
#         if "AssemblyUtil" not in self.clients:
#             if "devenv" in self.config and self.config["devenv"] == "1":
#                 from kbbasemodules import AssemblyUtil
#             else:
#                 from installed_clients.AssemblyUtilClient import AssemblyUtil
#             self.clients["AssemblyUtil"] = AssemblyUtil(self.callback_url,token=self.token)
#         return self.clients["AssemblyUtil"]
    
#     def rast_client(self):
#         if "RAST_SDK" not in self.clients:
#             if "devenv" in self.config and self.config["devenv"] == "1":
#                 from kbbasemodules import RAST_SDK
#             else:
#                 from installed_clients.RAST_SDKClient import RAST_SDK
#             self.clients["RAST_SDK"] = RAST_SDK(self.callback_url,token=self.token)
#         return self.clients["RAST_SDK"]
    
#     def anno_client(self,native_python_api=False):
#         if "cb_annotation_ontology_api" not in self.clients:
#             if native_python_api:
#                 from cb_annotation_ontology_api.annotation_ontology_api import AnnotationOntologyModule
#                 self.clients["cb_annotation_ontology_api"] = AnnotationOntologyModule("cb_annotation_ontology_api",{"data" :"/data/"},module_dir=self.module_dir+"/../cb_annotation_ontology_api",working_dir=self.working_dir,token=self.token,clients=self.clients,callback=self.callback_url)
#             else:
#                 if "devenv" in self.config and self.config["devenv"] == "1":
#                     from kbbasemodules import cb_annotation_ontology_api
#                 else:
#                     from installed_clients.cb_annotation_ontology_apiClient import cb_annotation_ontology_api
#                 self.clients["cb_annotation_ontology_api"] = cb_annotation_ontology_api(self.callback_url,token=self.token)
#         return self.clients["cb_annotation_ontology_api"]
    
#     #########GENERAL UTILITY FUNCTIONS#######################
#     def process_genome_list(self,input_references,workspace=None):
#         ws_identities = []
#         for ref in input_references:
#             ws_identities.append(self.process_ws_ids(ref,workspace))
#         output = self.ws_client().get_object_info(ws_identities,1)
#         output_references = []
#         for info in output:
#             if info[2].startswith("KBaseSearch.GenomeSet"):
#                 genomeset = self.get_object(self.wsinfo_to_ref(info))["data"]
#                 for label in genomeset["elements"]:
#                     output_references.append(genomeset["elements"][label]["ref"])
#             else:
#                 output_references.append(self.wsinfo_to_ref(info))
#         return output_references
    
#     def kb_version(self):
#         wsclient = self.ws_client()
#         if "appdev.kbase.us" in wsclient._client.url:
#             return "dev"
#         elif "/kbase.us" in wsclient._client.url:
#             return "prod"
#         elif "ci.kbase.us" in wsclient._client.url:
#             return "ci"
#         elif "next.kbase.us" in wsclient._client.url:
#             return "next"
#         else:
#             return "unknown"
    
#     def validate_args(self,params,required,defaults):
#         for item in required:
#             if item not in params:
#                 raise ValueError('Required argument '+item+' is missing!')
#         for key in defaults:
#             if key not in params:
#                 params[key] = defaults[key]
#         return params
    
#     def transfer_outputs(self,output,api_output,key_list):
#         for key in key_list:
#             if key in api_output:
#                 output[key] = api_output[key]
                
#     def print_json_debug_file(self,filename,data):
#         with open(self.working_dir+'/'+filename, 'w') as f:
#             json.dump(data, f,indent=4,skipkeys=True)
    
#     #########WORKSPACE RELATED FUNCTIONS#######################
#     def provenance(self):
#         return [{
#             'description': self.method,
#             'input_ws_objects': self.input_objects,
#             'method': self.method,
#             'script_command_line': "",
#             'method_params': [self.params],
#             'service': self.name,
#             'service_ver': self.version
#         }]
    
#     def set_ws(self,workspace):
#         if self.ws_id == workspace or self.ws_name == workspace:
#             return 
#         if not isinstance(workspace, str) or re.search('^\\d+$',workspace) != None:
#             if isinstance(workspace, str):
#                 workspace = int(workspace)
#             self.ws_id = workspace
#             info = self.ws_client().get_workspace_info({"id":workspace})
#             self.ws_name = info[1]
#         else:
#             self.ws_name = workspace
#             info = self.ws_client().get_workspace_info({"workspace":workspace})
#             self.ws_id = info[0]
    
#     def process_ws_ids(self,id_or_ref,workspace=None,no_ref=False):
#         """
#         IDs should always be processed through this function so we can interchangeably use
#         refs, IDs, and names for workspaces and objects
#         """
#         objspec = {}
#         if len(id_or_ref.split("/")) > 1:
#             if no_ref:
#                 array = id_or_ref.split("/")
#                 workspace = array[0]
#                 id_or_ref = array[1]
#             else:
#                 objspec["ref"] = id_or_ref
                
#         if "ref" not in objspec:
#             if isinstance(workspace, int):
#                 objspec['wsid'] = workspace
#             else:
#                 objspec['workspace'] = workspace
#             if isinstance(id_or_ref, int):
#                 objspec['objid'] = id_or_ref
#             else:
#                 objspec['name'] = id_or_ref
#         return objspec
          
#     def ws_get_objects(self, args):
#         """
#         All functions calling get_objects2 should call this function to ensure they get the retry
#         code because workspace periodically times out
#         :param args:
#         :return:
#         """
#         tries = 0
#         while tries < self.config["max_retry"]:
#             try:
#                 return self.ws_client().get_objects2(args)
#             except:
#                 logger.warning("Workspace get_objects2 call failed [%s:%s - %s]. Trying again!")
#                 tries += 1
#                 time.sleep(500)  # Give half second
#         logger.warning("get_objects2 failed after multiple tries: %s", sys.exc_info()[0])
#         raise

#     def get_object_info(self, id_or_ref, ws=None):
#         ws_identities = [self.process_ws_ids(id_or_ref, ws)]
#         return self.ws_client().get_object_info(ws_identities,1)[0]

#     def get_object(self, id_or_ref, ws=None):
#         res = self.ws_get_objects({"objects": [self.process_ws_ids(id_or_ref, ws)]})
#         if res is None:
#             return None
#         return res["data"][0]
    
#     def save_genome_or_metagenome(self,objid,workspace,obj_json):
#         self.set_ws(workspace)
#         save_output = self.gfu_client().save_one_genome({
#             "name" : objid,
#             "data" : obj_json,
#             "upgrade" : 1,
#             "provenance" : self.provenance(),
#             "hidden" : 0,
#             "workspace" : self.ws_name
#         });
#         self.obj_created.append({"ref":self.create_ref(objid,self.ws_name),"description":""})
#         return save_output["info"]
    
#     def save_ws_object(self,objid,workspace,obj_json,obj_type):
#         self.set_ws(workspace)
#         params = {
#             'id':self.ws_id,
#             'objects': [{
#                 'data': obj_json,
#                 'name': objid,
#                 'type': obj_type,
#                 'meta': {},
#                 'provenance': self.provenance()
#             }]
#         }
#         self.obj_created.append({"ref":self.create_ref(objid,self.ws_name),"description":""})
#         return self.ws_client().save_objects(params)
    
#     def wsinfo_to_ref(self,info):
#         return str(info[6])+"/"+str(info[0])+"/"+str(info[4])
    
#     def create_ref(self,id_or_ref,ws=None):
#         if isinstance(id_or_ref, int):
#             id_or_ref=str(id_or_ref)
#         if len(id_or_ref.split("/")) > 1:
#             return id_or_ref
#         if isinstance(ws, int):
#             ws=str(ws)
#         return ws+"/"+id_or_ref
    
#     #########REPORT RELATED FUNCTIONS#######################
#     def save_report_to_kbase(self,height=700,message="",warnings=[],file_links=[],summary_height=None):
#         rootDir = self.working_dir+"/html/"
#         files = [{'path': "/kb/module/work/tmp/html/",'name': "index.html",'description': 'HTML report'}]
#         for dirName, subdirList, fileList in os.walk(rootDir):
#             for fname in fileList:
#                 if fname != "index.html":
#                     files.append({'path': dirName.replace(rootDir,"/kb/module/work/tmp/html/"),'name': fname,'description': 'Files related to HTML report'})
#         report_name = self.method+"-"+str(uuid.uuid4())
#         output = self.report_client().create_extended_report({
#             'message': message,
#             'warnings': warnings,            
#             'html_links': files,
#             'file_links': file_links,
#             'direct_html_link_index': 0,
#             'html_window_height': height,
#             'objects_created': self.obj_created,
#             'workspace_name': self.ws_name,
#             'report_object_name': report_name,
#             'summary_window_height': summary_height
#         })
#         return {"report_name":report_name,"report_ref":output["ref"],'workspace_name':self.ws_name}
    
#     def build_dataframe_report(self,table,column_list):        
#         #Convert columns to this format:
#         columns = []
#         for item in column_list:
#             columns.append({"data":item})
#         #for index, row in table.iterrows():
#         #    pass
#         json_str = table.to_json(orient='records')
#         #columns=column_list
#         html_data = """
#     <html>
#     <header>
#         <link href="https://cdn.datatables.net/1.11.5/css/jquery.dataTables.min.css" rel="stylesheet">
#     </header>
#     <body>
#     <script src="https://code.jquery.com/jquery-3.6.0.slim.min.js" integrity="sha256-u7e5khyithlIdTpu22PHhENmPcRdFiHRjhAuHcs05RI=" crossorigin="anonymous"></script>
#     <script type="text/javascript" src="https://cdn.datatables.net/1.11.5/js/jquery.dataTables.min.js"></script>
#     <script>
#         $(document).ready(function() {
#             $('#example').DataTable( {
#                 "ajax": {
#                     "url": "data.json"
#                 },
#                 "columns": """+json.dumps(columns,indent=4)+"""
#             } );
#         } );
#     </script>
#     </body>
#     </html>
#     """
#         os.makedirs(self.working_dir+"/html", exist_ok=True)
#         with open(self.working_dir+"/html/index.html", 'w') as f:
#             f.write(html_data)
#         with open(self.working_dir+"/html/data.json", 'w') as f:
#             f.write(json_str)


# class KBDevUtils(BaseModule):
#     def __init__(self,study_name,token_file=str(Path.home()) + '/.kbase/token',ws_version="prod",sdkhome=None,output_root=None):
#         wsurl = None
#         if ws_version == "prod":
#             wsurl = "https://kbase.us/services/ws"
#         elif ws_version == "appdev":
#             wsurl = "https://appdev.kbase.us/services/ws"
#         elif ws_version == "ci":
#             wsurl = "https://ci.kbase.us/services/ws"    
#         token = None
#         with open(token_file, 'r') as fh:
#             token = fh.read().strip()
#         self.root_path = config.get("callback_path","")    
#         if output_root:
#             self.output_root = output_root
#         else:
#             self.output_root = config.get("output_root","KBaseAnalysis/")
#         if self.output_root[0] != "/":
#             self.output_root = str(Path.home())+"/"+self.output_root
#         if not sdkhome:
#             if ws_version == "prod":
#                 sdkhome = config.get("prod_sdk_home","prod")
#             elif ws_version == "appdev":
#                 sdkhome = config.get("appdev_sdk_home","appdev")
#         self.callback_file = self.root_path+"/"+sdkhome+"/kb_sdk_home/run_local/workdir/CallBack.txt"
#         self.working_directory = self.root_path+"/"+sdkhome+"/kb_sdk_home/run_local/workdir/tmp/"
#         callback = None
#         if exists(self.callback_file):
#             with open(self.callback_file, 'r') as fh:
#                 callback = fh.read()     
#         BaseModule.__init__(self,"KBDevUtils."+study_name,config,codebase+"/chenry_utility_module/",str(Path.home()) + "/scratch/" + study_name,token,{"Workspace":Workspace(wsurl, token=token)},callback)
#         self.version = "0.1.1.kbdu"
#         self.study_name = study_name
#         print("Output files printed to:"+self.out_dir()+" when using KBDevUtils.out_dir()")
#         self.msrecon = None
    
#     def msseedrecon(self):
#         if self.msrecon == None:
#             from ModelSEEDReconstruction.modelseedrecon import ModelSEEDRecon
#             self.msrecon = ModelSEEDRecon(self.config,codebase+"/KB-ModelSEEDReconstruction/",self.working_dir,self.token,self.clients,self.callback_url)
#         return self.msrecon
    
#     def devutil_client(self):
#         if "KBDevUtils" not in self.clients:
#             from installed_clients.chenry_utility_moduleClient import chenry_utility_module
#             self.clients["KBDevUtils"] = chenry_utility_module(self.callback_url,token=self.token)
#         return self.clients["KBDevUtils"]
        
#     def clear_sdk_dir(self):
#         return self.devutil_client().run_command({"command":"clear"})
        
#     def sdk_dir_perms(self):
#         return self.devutil_client().run_command({"command":"perms"})
    
#     def out_dir(self,create=True):
#         output_path = self.output_root+"/"+self.study_name+"/"
#         if create:
#             if not exists(output_path):
#                 os.makedirs(output_path, exist_ok=True)
#         return output_path
        

        
        
class CommScoresUtil():
    @staticmethod
    def _categorize_mets(metIDs):
        # load compound categories
        package_dir = os.path.abspath(os.path.dirname(__file__))
        categories_dir = os.path.join(package_dir, "data", "categories")
        sugars, aminoacids = (
            load(os.path.join(categories_dir, "sugars.npy")),
            load(os.path.join(categories_dir, "aminoAcids.npy")),
        )
        vitamins, minerals = (
            load(os.path.join(categories_dir, "vitamins.npy")),
            load(os.path.join(categories_dir, "minerals.npy")),
        )
        energy_compounds = load(os.path.join(categories_dir, "energy_compounds.npy"))
        energy_compounds_dic = dict(zip(energy_compounds[:, 0], energy_compounds[:, 1]))
        sugars_dic, aminoacids_dic = (
            dict(zip(sugars[:, 0], sugars[:, 1])),
            dict(zip(aminoacids[:, 0], aminoacids[:, 1])),
        )
        vitamins_dic, minerals_dic = (
            dict(zip(vitamins[:, 0], vitamins[:, 1])),
            dict(zip(minerals[:, 0], minerals[:, 1])),
        )
        
        # partition the compounds into their respective directories
        met_sugars, met_aminoAcids, met_vitamins, met_minerals, met_energy, met_other = ([], [], [], [], [], [])
        for metID in metIDs:
            if metID in sugars[:, 0]:
                met_sugars.append(f"{sugars_dic[metID]} ({metID})")
            elif metID in aminoacids[:, 0]:
                met_aminoAcids.append(f"{aminoacids_dic[metID]} ({metID})")
            elif metID in vitamins[:, 0]:
                met_vitamins.append(f"{vitamins_dic[metID]} ({metID})")
            elif metID in minerals[:, 0]:
                met_minerals.append(f"{minerals_dic[metID]} ({metID})")
            elif metID in energy_compounds[:, 0]:
                met_energy.append(f"{energy_compounds_dic[metID]} ({metID})")
            else:
                met_other.append(metID)
        return met_sugars, met_aminoAcids, met_vitamins, met_minerals, met_energy, met_other


    @staticmethod
    def remove_metadata(element):
        rm_costless = re.compile("(\s\(.+\))")
        try:
            element = float(rm_costless.sub("", str(element)).replace("%", ""))
        except:
            pass
        return element

    @staticmethod
    def convert_to_int(element):
        try:
            element = int(element)
        except:
            pass
        return element

    @staticmethod
    def _process_mets(metIDs):
        return [", ".join(lst) for lst in CommScoresUtil._categorize_mets(metIDs)]

    @staticmethod
    def _compatibilize(member_models: Iterable, printing=False):
        return member_models

    @staticmethod
    def _load_models(
        member_models: Iterable, com_model=None, compatibilize=True, commID=None, printing=False
    ):
        # ic(member_models, com_model, compatibilize)
        if not com_model and member_models:
            return member_models, build_from_species_models(member_models, commID, "CommScores_community")  # (model, names=names, abundances=abundances)
        # elif com_model and not member_models:
        #     return com_model.members, com_model  # TODO the individual models of a community model can be parsed
        if compatibilize:
            return (
                CommScoresUtil._compatibilize(member_models, printing),
                CommScoresUtil._compatibilize([com_model], printing)[0],
            )
        return member_models, com_model


    @staticmethod
    def _get_media(
        media=None,
        com_model=None,
        model_s_=None,
        min_growth=0.1,
        environment=None,
        interacting=True,
        printing=False,
        minimization_method="minFlux",
    ):
        assert com_model is not None or model_s_ is not None, "com_model or model_s_ must be parameterized."
        if com_model is True and isinstance(model_s_, (set, list, tuple)):
                com_model = build_from_species_models(model_s_)
        if media is not None:
            if model_s_ is not None and not isinstance(model_s_, (list, set, tuple)):
                return media["members"][model_s_.id]["media"]
            elif com_model is not None:
                return media["community_media"]
            return media
        com_media = None
        if com_model is not None:
            minGrowth = min_growth
            while com_media is None or com_media == {}:
                com_media, media_sol = MSMinimalMedia.determine_min_media(
                    com_model, minimization_method, minGrowth, environment, interacting, printing)
                minGrowth *= 1.1
            if minGrowth != min_growth:   print(f"{com_model.id} needed {minGrowth} v {min_growth}\n")
            if model_s_ is None:
                return com_media, media_sol
        if model_s_ is not None:
            min_media = None
            if not isinstance(model_s_, (list, set, tuple, ndarray)):
                minGrowth = min_growth
                while min_media is None or min_media == {}:
                    min_media, media_sol = MSMinimalMedia.determine_min_media(
                        model_s_, minimization_method, minGrowth, environment, interacting, printing)
                    minGrowth *= 1.1
                return min_media, media_sol
            members_media = {}
            for model in model_s_:
                minGrowth = min_growth
                while min_media is None or min_media == {}:
                    min_media, media_sol = MSMinimalMedia.determine_min_media(
                        model, minimization_method, minGrowth, environment, interacting, printing)
                    minGrowth *= 1.1
                members_media[model.id] = {"media": min_media, "solution": media_sol}
                min_media = None
                if minGrowth != min_growth:  print(f"{model.id} needed {minGrowth} v {min_growth}")
            if com_model is None:
                return members_media
            return {"community_media": com_media, "members": members_media}
        raise BadModels(f"The parameterized community model of type {type(com_model)} and member models {model_s_} are not properly captured.")

    @staticmethod
    def _sigfig_check(value, sigfigs, default):
        if str(value) in ["inf", "nan"]:
            value = ""
        if FBAHelper.isnumber(value):
            return sigfig.round(value, sigfigs)
        else:
            return default

    @staticmethod
    def _calculate_jaccard_score(set1, set2):
        if set1 == set2:
            print(f"The sets are identical, with a length of {len(set1)}.")
        if len(set1.union(set2)) == 0:
            return (None, None)
        return (
            set1.intersection(set2),
            len(set1.intersection(set2)) / len(set1.union(set2)),
        )


    @staticmethod
    def _check_model(model_util, media, model_str=None):
        default_media = model_util.model.medium
        # print("test")
        if media is not None:
            model_util.add_medium(media)
        obj_val = model_util.model.slim_optimize()
        model_str = model_str or model_util.model.id
        # print(model_str, obj_val)
        if isclose(obj_val, 0, abs_tol=1e-6) or not FBAHelper.isnumber(obj_val):
            print(f"The {model_str} model is not operational on the provided media")
            model_util.add_medium(default_media)
        return model_util.model


    @staticmethod
    def _determine_growths(modelUtils, environ, sigfig=5):
        obj_vals = []
        for util in modelUtils:
            util.add_medium(environ)
            val = util.model.slim_optimize()
            sigfiggedVal = CommScoresUtil._sigfig_check(val, sigfig, "")
            obj_vals.append(sigfiggedVal)
        return obj_vals

    @staticmethod
    def nanFilter(value, string=True):
        if isinstance(value, str) or value is None:
            if string:
                return value
            else:
                return nan
        if any([value < 0, value > 1e5]):
            return "" if string else nan
        return value

    
util = CommScoresUtil()
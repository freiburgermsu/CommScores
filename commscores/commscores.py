# from icecream import ic
import os
import re
import warnings
from collections import Counter
from itertools import chain, combinations, permutations
from math import isclose
from pprint import pprint
from typing import Iterable

from deepdiff import DeepDiff  # (old, new)
from modelseedpy.community.commhelper import build_from_species_models
from modelseedpy.core.exceptions import ObjectiveError, ParameterError
from modelseedpy.core.fbahelper import FBAHelper
from modelseedpy.core.msminimalmedia import MSMinimalMedia
from modelseedpy.core.msmodelutl import MSModelUtil

# silence deprecation warnings from DeepDiff parsing the syntrophy
from numpy import array

from .logger import logger
from .scores import antiSMASH, bss, cip, fs, gyd, mip, mp, mro, mu, pc, sc, smetana

warnings.simplefilter("ignore", category=DeprecationWarning)

package_dir = os.path.abspath(os.path.dirname(__file__))

rm_comp = FBAHelper.remove_compartment

from .utils import (
    _categorize_mets,
    _compatibilize,
    _get_media,
    _load_models,
    _process_mets,
    _sigfig_check,
    convert_to_int,
    nanFilter,
    remove_metadata,
)


class CommScores:
    def __init__(
        self,
        member_models,
        min_growth=0.1,
        n_solutions=100,
        environment=None,
        abstol=1e-3,
        media_dict=None,
        printing=True,
        raw_content=False,
        antismash_json_path: str = None,
        antismash_zip_path: str = None,
        minimal_media_method="minFlux",
    ):
        self.min_growth = min_growth
        self.abstol = abstol
        self.n_solutions = n_solutions
        self.printing = printing
        self.raw_content = raw_content
        self.antismash_json_path = antismash_json_path
        self.antismash_zip_path = antismash_zip_path

        # process the models
        self.models = _compatibilize(member_models)
        self.community = MSModelUtil(build_from_species_models(self.models))
        ## define the environment
        if environment:
            if hasattr(environment, "get_media_constraints"):
                ### standardize modelseed media into COBRApy media
                environment = {
                    "EX_" + exID: -bound[0]
                    for exID, bound in environment.get_media_constraints().items()
                }
            self.community.add_medium(environment)
        self.environment = environment
        ## test growth
        for model in self.models:
            if model.slim_optimize() == 0:
                raise ObjectiveError(
                    f"The model {model.id} possesses an objective value of 0 in complete media, "
                    "which is incompatible with minimal media computations and hence SMETANA."
                )
        if self.community.model.slim_optimize() == 0:
            raise ObjectiveError(
                f"The community model {self.community.model.id} possesses an objective "
                "value of 0 in complete media, which is incompatible with minimal "
                "media computations and hence SMETANA."
            )
        ## determine the minimal media for each model, including the community
        self.media = (
            media_dict
            if media_dict
            else MSMinimalMedia.comm_media_est(
                member_models,
                self.community.model,
                minimal_media_method,
                min_growth,
                self.environment,
                True,
                n_solutions,
                printing,
            )
        )

    def all_scores(
        self,
        mp_score=True,
        kbase_obj=None,
        token_string: str = None,
        kbase_token_path: str = None,
        annotated_genomes: dict = None,
    ):
        mro = self.mro_score()
        mip = self.mip_score(interacting_media=self.media)
        mp = None if not mp_score else self.mp_score()
        mu = None  # self.mu_score()
        sc = None  # self.sc_score()
        smetana = None  # self.smetana_score()
        gyd = self.gyd_score()
        fs = (
            self.fs_score()
            if any(
                [
                    kbase_obj is not None,
                    annotated_genomes != [],
                    token_string is not None and kbase_token_path is not None,
                ]
            )
            else None
        )
        return {
            "mro": mro,
            "mip": mip,
            "mp": mp,
            "mu": mu,
            "sc": sc,
            "smetana": smetana,
            "gyd": gyd,
            "fs": fs,
        }

    def mro_score(self):
        self.mro_val = mro(
            self.models,
            self.media["members"],
            self.min_growth,
            self.media,
            self.raw_content,
            self.environment,
            self.printing,
            True,
        )
        if not self.printing:
            return self.mro_val
        if self.raw_content:
            for pair, (interaction, media) in self.mro_val.items():
                newcomer, established = pair.split("---")
                print(
                    f"\n(MRO) The {newcomer} media {media} possesses {interaction} shared "
                    f"requirements with the {established} established member."
                )
                return self.mro_val
        for pair, mro in self.mro_val.items():
            newcomer, established = pair.split("---")
            print(
                f"\nThe {newcomer} on {established} MRO score: {mro[0]} ({mro[0]*100:.2f}%). "
                f"This is the percent of nutritional requirements in {newcomer} "
                f"that overlap with {established} ({mro[1]}/{mro[2]})."
            )
        return self.mro_val

    def mip_score(
        self, interacting_media: dict = None, noninteracting_media: dict = None
    ):
        interacting_media = interacting_media or self.media or None
        diff, self.mip_val = mip(
            self.models,
            self.community.model,
            self.min_growth,
            interacting_media,
            noninteracting_media,
            self.environment,
            self.printing,
            True,
        )
        if not self.printing:
            return self.mip_val
        print(
            f"\nMIP score: {self.mip_val}\t\t\t{self.mip_val} required compound(s) can be sourced via syntrophy:"
        )
        if self.raw_content:
            pprint(diff)
        return self.mip_val

    def gyd_score(self, coculture_growth=False):
        self.gyd_val = gyd(
            self.models, environment=self.environment, coculture_growth=coculture_growth
        )
        if not self.printing:
            return self.gyd
        growth_type = "monocultural" if not coculture_growth else "cocultural"
        for pair, score in self.gyd_val.items():
            print(
                f"\nGYD score: The {growth_type} growth difference between the {pair} member models"
                f" is {score} times greater than the growth of the slower member."
            )
        return self.gyd

    def fs_score(
        self,
        kbase_obj=None,
        token_string: str = None,
        kbase_token_path: str = None,
        annotated_genomes: dict = None,
    ):
        self.fs_val = fs(
            self.models, kbase_obj, token_string, kbase_token_path, annotated_genomes
        )
        if not self.printing:
            return fs
        for pair, score in self.fs_val.items():
            print(
                f"\nFS Score: The similarity of RAST functional SSO ontology "
                f"terms between the {pair} members is {score}."
            )
        return fs

    def mp_score(self):
        print("executing MP")
        self.mp_val = mp(
            self.models,
            self.environment,
            self.community.model,
            None,
            self.abstol,
            self.printing,
        )
        if not self.printing:
            return self.mp_val
        if self.raw_content:
            print(
                "\n(MP) The possible contributions of each member in the member media include:\n"
            )
            pprint(self.mp_val)
        else:
            print(
                "\nMP score:\t\t\tEach member can possibly contribute the following to the community:\n"
            )
            for member, contributions in self.mp_val.items():
                print(member, "\t", len(contributions))
        return self.mp_val

    def mu_score(self):
        member_excreta = self.mp_score() if not hasattr(self, "mp_val") else self.mp_val
        self.mu_val = mu(
            self.models,
            self.environment,
            member_excreta,
            self.n_solutions,
            self.abstol,
            True,
            self.printing,
        )
        if not self.printing:
            return self.mu_val
        print(
            "\nMU score:\t\t\tThe fraction of solutions in which each member is the "
            "syntrophic receiver that contain a respective metabolite:\n"
        )
        pprint(self.mu_val)
        return self.mu_val

    def sc_score(self):
        self.sc_val = sc(
            self.models,
            self.community.model,
            self.min_growth,
            self.n_solutions,
            self.abstol,
            True,
            self.printing,
        )
        if not self.printing:
            return self.sc_val
        print(
            "\nSC score:\t\t\tThe fraction of community members who syntrophically contribute to each species:\n"
        )
        pprint(self.sc_val)
        return self.sc_val

    def smetana_score(self):
        if not hasattr(self, "sc_val"):
            self.sc_val = self.sc_score()
        sc_coupling = all(array(list(sc.values())) is not None)
        if not hasattr(self, "mu_val"):
            self.mu_val = self.mu_score()
        if not hasattr(self, "mp_val"):
            self.mp_val = self.mp_score()

        self.smetana = smetana(
            self.models,
            self.community.model,
            self.min_growth,
            self.n_solutions,
            self.abstol,
            (self.sc_val, self.mu_val, self.mp_val),
            True,
            sc_coupling,
            self.printing,
        )
        if self.printing:
            print("\nsmetana score:\n")
            pprint(self.smetana)
        return self.smetana

    def antiSMASH_scores(self, antismash_json_path=None):
        self.antismash = antiSMASH(antismash_json_path or self.antismash_json_path)
        if not self.printing:
            return self.antismash
        if self.raw_content:
            print(
                "\n(antismash) The biosynthetic_areas, BGCs, protein_annotations, clusterBlast, and "
                "num_clusterBlast from the provided antiSMASH results:\n"
            )
            print(
                "The 'areas' that antiSMASH determines produce biosynthetic products:"
            )
            pprint(self.antismash[0])
            print("The set of biosynthetic gene clusters:")
            pprint(self.antismash[1])
            print("The set of clusterblast protein annotations:")
            pprint(self.antismash[2])
            print("Resistance information from clusterblast")
            pprint(self.antismash[3])
            print("The number of proteins associated with resistance")
            pprint(self.antismash[4])
            return self.antismash
        print("\nantiSMASH scores:\n")
        print(
            "The community exhibited:"
            f"- {len(self.antismash[0])}'areas' that antiSMASH determines produce biosynthetic products."
            f"- {len(self.antismash[1])} biosynthetic gene clusters."
            f"- {len(self.antismash[2])} clusterblast protein annotations."
            f"- {len(self.antismash[3])} parcels of resistance information from clusterblast."
            f"- {self.antismash[4]} proteins associated with resistance."
        )
        return list(map(len, self.antismash[:4])) + [self.antismash[4]]

    ###### STATIC METHODS OF THE SMETANA SCORES, WHICH ARE APPLIED IN THE ABOVE CLASS OBJECT ######

    @staticmethod
    def mqs():
        pass

    @staticmethod
    def get_all_genomes_from_ws(
        ws_id,
        kbase_object=None,
        cobrakbase_repo_path: str = None,
        kbase_token_path: str = None,
    ):
        def get_genome(genome_name):
            return kbase_object.ws_client.get_objects2(
                {"objects": [{"ref": f"{ws_id}/{genome_name}"}]}
            )["data"][0]["data"]

        # load the kbase client instance
        if not kbase_object:
            import os

            os.environ["HOME"] = cobrakbase_repo_path
            import cobrakbase

            with open(kbase_token_path) as token_file:
                kbase_object = cobrakbase.KBaseAPI(token_file.readline())

        # calculate the complementarity
        genome_list = kbase_object.ws_client.list_objects(
            {
                "ids": [ws_id],
                "type": "KBaseGenomes.Genome",
                "minObjectID": 0,
                "maxObjectID": 10000,
            }
        )
        genome_names = [g[1] for g in genome_list if g[1].endswith("RAST")]
        return {
            genome_name: set(
                [
                    sso
                    for j in get_genome(genome_name)["cdss"]
                    for sso in j["ontology_terms"]["SSO"].keys()
                ]
            )
            for genome_name in genome_names
        }

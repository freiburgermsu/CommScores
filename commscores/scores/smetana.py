from collections import Counter
from itertools import chain, combinations, permutations
from typing import Iterable

from deepdiff import DeepDiff
from modelseedpy.community.commhelper import build_from_species_models
from modelseedpy.core.msminimalmedia import MSMinimalMedia
from modelseedpy.core.msmodelutl import MSModelUtil
from optlang import Constraint, Objective, Variable

from commscores.utils import _compatibilize, _load_models


def contributions(org_possible_contributions, scores, model_util, abstol):
    # identify and log excreta from the solution
    model_util.add_objective(
        sum(ex_rxn.flux_expression for ex_rxn in org_possible_contributions)
    )
    sol = model_util.model.optimize()
    if sol.status != "optimal":
        # exit the while loop by returning the original possible_contributions,
        ## hence DeepDiff == {} and the while loop terminates
        return scores, org_possible_contributions
    # identify and log excreta from the solution
    possible_contributions = org_possible_contributions[:]
    for ex in org_possible_contributions:
        if ex.id in sol.fluxes.keys() and sol.fluxes[ex.id] >= abstol:
            possible_contributions.remove(ex)
            scores[model_util.model.id].update([met.id for met in ex.metabolites])
    return scores, possible_contributions


def mu(
    member_models: Iterable,
    environment=None,
    member_excreta=None,
    n_solutions=100,
    abstol=1e-3,
    compatibilized=False,
    printing=True,
):
    """the fractional frequency of each received metabolite amongst all possible alternative syntrophic solutions"""
    # member_solutions = member_solutions if member_solutions else {model.id: model.optimize() for model in member_models}
    scores = {}
    member_models = (
        member_models if compatibilized else _compatibilize(member_models, printing)
    )
    if member_excreta:
        missing_members = [
            model for model in member_models if model.id not in member_excreta
        ]
        if missing_members:
            print(
                f"The {','.join(missing_members)} members are missing from the defined "
                f"excreta list and will therefore be determined through an additional MP simulation."
            )
            member_excreta.update(mp(missing_members, environment))
    else:
        member_excreta = mp(member_models, environment, None, abstol, printing)
    for org_model in member_models:
        other_excreta = set(
            chain.from_iterable(
                [
                    excreta
                    for model, excreta in member_excreta.items()
                    if model != org_model.id
                ]
            )
        )
        print(f"\n{org_model.id}\tOther Excreta", other_excreta)
        model_util = MSModelUtil(org_model, True)
        if environment:
            model_util.add_medium(environment)
        ex_rxns = {
            ex_rxn: list(ex_rxn.metabolites)[0] for ex_rxn in model_util.exchange_list()
        }
        print(f"\n{org_model.id}\tExtracellular reactions", ex_rxns)
        variables = {
            ex_rxn.id: Variable(
                "___".join([model_util.model.id, ex_rxn.id]),
                lb=0,
                ub=1,
                type="binary",
            )
            for ex_rxn in ex_rxns
        }
        model_util.add_cons_vars(list(variables.values()))
        media, solutions = [], []
        sol = model_util.model.optimize()
        while sol.status == "optimal" and len(solutions) < n_solutions:
            solutions.append(sol)
            medium = set(
                [
                    ex
                    for ex in ex_rxns
                    if sol.fluxes[ex.id] < -abstol and ex in other_excreta
                ]
            )
            model_util.create_constraint(
                Constraint(
                    sum([variables[ex.id] for ex in medium]),
                    ub=len(medium) - 1,
                    name=f"iteration_{len(solutions)}",
                )
            )
            media.append(medium)
            sol = model_util.model.optimize()
        counter = Counter(chain(*media))
        scores[model_util.model.id] = {
            met.id: counter[ex] / len(media)
            for ex, met in ex_rxns.items()
            if counter[ex] > 0
        }
    return scores


def mp(
    member_models: Iterable,
    environment,
    com_model=None,
    minimal_media=None,
    abstol=1e-3,
    printing=False,
):
    """Discover the metabolites that each species can contribute to a community"""
    community = (
        _compatibilize(com_model)
        if com_model
        else build_from_species_models(member_models, standardize=True)
    )
    community.medium = minimal_media or MSMinimalMedia.minimize_flux(community)
    scores = {}
    for org_model in (
        member_models
    ):  # TODO support parsing the individual members through the MSCommunity object
        model_util = MSModelUtil(org_model)
        model_util.compatibilize(printing=printing)
        if environment:
            model_util.add_medium(environment)
        scores[model_util.model.id] = set()
        # determines possible member contributions in the community environment, where the excretion of media compounds is irrelevant
        org_possible_contr = [
            ex_rxn
            for ex_rxn in model_util.exchange_list()
            if (ex_rxn.id not in community.medium and ex_rxn.upper_bound > 0)
        ]
        # ic(org_possible_contributions, len(model_util.exchange_list()), len(community.medium))
        scores, possible_contr = contributions(
            org_possible_contr, scores, model_util, abstol
        )
        while DeepDiff(org_possible_contr, possible_contr):
            print("remaining possible_contributions", len(possible_contr), end="\r")
            ## optimize the sum of the remaining exchanges that have not surpassed the abstol
            org_possible_contr = possible_contr[:]
            scores, possible_contr = contributions(
                org_possible_contr, scores, model_util, abstol
            )

        ## individually checks the remaining possible contributions
        for ex_rxn in possible_contr:
            model_util.model.objective = Objective(ex_rxn.flux_expression)
            sol = model_util.model.optimize()
            if sol.status == "optimal" or sol.objective_value > abstol:
                for met in ex_rxn.metabolites:
                    if met.id in scores[model_util.model.id]:
                        scores[model_util.model.id].remove(met.id)
                        print("removing", met.id)
    return scores


def sc(
    member_models: Iterable = None,
    com_model=None,
    min_growth=0.1,
    n_solutions=100,
    abstol=1e-6,
    compatibilized=True,
    printing=False,
):
    """Calculate the frequency of interspecies dependency in a community"""
    member_models, community = _load_models(
        member_models, com_model, not compatibilized, printing=printing
    )
    for rxn in com_model.reactions:
        rxn.lower_bound = 0 if "bio" in rxn.id else rxn.lower_bound

    # c_{rxn.id}_lb: rxn < 1000*y_{species_id}
    # c_{rxn.id}_ub: rxn > -1000*y_{species_id}
    variables = {}
    constraints = []
    # TODO this can be converted to an MSCommunity object by looping through each index
    # leverage CommKinetics
    for org_model in member_models:
        model_util = MSModelUtil(org_model, True)
        variables[model_util.model.id] = Variable(
            name=f"y_{model_util.model.id}", lb=0, ub=1, type="binary"
        )
        model_util.add_cons_vars([variables[model_util.model.id]])
        for rxn in model_util.model.reactions:
            if "bio" not in rxn.id:
                # print(rxn.flux_expression)
                lb = Constraint(
                    rxn.flux_expression + 1000 * variables[model_util.model.id],
                    name="_".join(["c", model_util.model.id, rxn.id, "lb"]),
                    lb=0,
                )
                ub = Constraint(
                    rxn.flux_expression - 1000 * variables[model_util.model.id],
                    name="_".join(["c", model_util.model.id, rxn.id, "ub"]),
                    ub=0,
                )
                constraints.extend([lb, ub])

    # calculate the SCS
    scores = {}
    for model in member_models:
        com_model_util = MSModelUtil(com_model)
        com_model_util.add_cons_vars(constraints, sloppy=True)
        # model growth is guaranteed while minimizing the growing members of the community
        ## SMETANA_Biomass: {biomass_reactions} > {min_growth}
        com_model_util.create_constraint(
            Constraint(
                sum(rxn.flux_expression for rxn in model.reactions if "bio" in rxn.id),
                name="SMETANA_Biomass",
                lb=min_growth,
            )
        )  # sloppy = True)
        other_members = [other for other in member_models if other.id != model.id]
        com_model_util.add_objective(
            sum([variables[other.id] for other in other_members]), "min"
        )
        previous_constraints, donors_list = [], []
        for i in range(n_solutions):
            sol = com_model.optimize()  # FIXME The solution is not optimal
            if sol.status != "optimal":
                scores[model.id] = None
                break
            donors = [
                o
                for o in other_members
                if com_model.solver.primal_values[f"y_{o.id}"] > abstol
            ]
            donors_list.append(donors)
            previous_con = f"iteration_{i}"
            previous_constraints.append(previous_con)
            com_model_util.add_cons_vars(
                [
                    Constraint(
                        sum(variables[o.id] for o in donors),
                        name=previous_con,
                        ub=len(previous_constraints) - 1,
                    )
                ],
                sloppy=True,
            )
        if i != 0:
            donors_counter = Counter(chain(*donors_list))
            scores[model.id] = {
                o.id: donors_counter[o] / len(donors_list) for o in other_members
            }
    return scores


def smetana(
    member_models: Iterable,
    environment,
    com_model=None,
    min_growth=0.1,
    n_solutions=100,
    abstol=1e-6,
    prior_values=None,
    compatibilized=False,
    sc_coupling=False,
    printing=False,
):
    """Quantifies the extent of syntrophy as the sum of all exchanges in a given nutritional environment"""
    member_models, community = _load_models(
        member_models, com_model, compatibilized == False, printing=printing
    )
    sc = None
    if not prior_values:
        mp = mp(member_models, environment, com_model, abstol)
        mu = mu(member_models, environment, mp, n_solutions, abstol, compatibilized)
        if sc_coupling:
            sc = sc(
                member_models,
                com_model,
                min_growth,
                n_solutions,
                abstol,
                compatibilized,
            )
    elif len(prior_values) == 3:
        sc, mu, mp = prior_values
    else:
        mu, mp = prior_values

    smetana_scores = {}
    for pairs in combinations(member_models, 2):
        for model1, model2 in permutations(pairs):
            if model1.id not in smetana_scores:
                smetana_scores[model1.id] = {}
            if not any([not mu[model1.id], not mp[model1.id]]):
                sc_score = 1 if not sc_coupling else sc[model1.id][model2.id]
                models_mets = list(model1.metabolites) + list(model2.metabolites)
                unique_mets = set([met.id for met in models_mets])
                smetana_scores[model1.id][model2.id] = 0
                for met in models_mets:
                    if met.id in unique_mets:
                        mp_score = 0 if met.id not in mp[model1.id] else 1
                        smetana_scores[model1.id][model2.id] += (
                            mu[model1.id].get(met.id, 0) * sc_score * mp_score
                        )
    return smetana_scores

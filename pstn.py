"""
"""

from __future__ import annotations

#################################################################################

from typing import Callable, Dict, FrozenSet, Iterable, List, NamedTuple, Optional, Set, Tuple

import mip
import pycryptosat
import statistics

from common import Node, Weight
from stnu import STNU, check_stnu_dc

#################################################################################

class ProbDistr(NamedTuple):
    """
    Represents a probability distribution with a finite support / finite (truncated) definition domain.
    
    The representation uses the cumulative distribution function (CDF) instead
    of the probability distribution function (PDF).
    """

    cumul_distr_func: Callable[[float], float]

    definition_domain: Tuple[float, float]

    discr_points: Tuple[float,...]

    @classmethod
    def uniform(cls,
        definition_domain: Tuple[float, float],
    ):
        
        if definition_domain[0] > definition_domain[1]:
            raise ValueError("Lower bound of distribution domain is greater than upper bound")

        def cumul_distr_func(x):
            if x < definition_domain[0]:
                return 0
            elif x >= definition_domain[1]:
                return 1
            else:
                assert definition_domain[1] > definition_domain[0]
                return (x - definition_domain[0]) / (definition_domain[1] - definition_domain[0])

        return ProbDistr(cumul_distr_func,
                         definition_domain,
                         definition_domain)

    @classmethod
    def truncated_normal(cls,
        mu: float,
        sigma: float,
        definition_domain: Tuple[float, float],
        num_discr_points: int
    ):
        
        if definition_domain[0] > definition_domain[1]:
            raise ValueError("Lower bound of distribution domain is greater than upper bound")

        normal = statistics.NormalDist(0, 1)
        alpha = (definition_domain[0] - mu) / sigma
        beta = (definition_domain[1] - mu) / sigma

        def cumul_distr_func(x):
            return (normal.cdf((x - mu) / sigma) - normal.cdf(alpha)) / (normal.cdf(beta) - normal.cdf(alpha))

        discr_step = (definition_domain[1]-definition_domain[0]) / num_discr_points
        discr_points = tuple(definition_domain[0] + i*discr_step for i in range(num_discr_points))

        return ProbDistr(cumul_distr_func,
                         definition_domain,
                         discr_points)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class PSTN(NamedTuple):

    reverse_graph: Dict[Node, Set[Node]]
    """
    Contains the (reverse) graph representation of the PSTN.
    The key of the dictionary corresponds to the target node, and
    its value is the set of source nodes from which an activities
    or requirements targeting the target node are defined. 
    """

    requirements: Dict[Tuple[Node, Node], Tuple[Weight, Weight]]
    """
    Contains the requirements links of the PSTN, which have a slightly
    different meaning/interpretation/semantics than that of activity links.
    """

    activities: Dict[Tuple[Node, Node], Tuple[Weight, Weight] | ProbDistr]
    """
    Contains both the probabilistic and non probabilistic activity links of the PSTN.
    """

    @classmethod
    def from_requirements_and_activities(cls,
        requirements: Iterable[Tuple[Node, Tuple[Weight, Weight], Node]],
        activities: Iterable[Tuple[Node, Tuple[Weight, Weight] | ProbDistr, Node]],
    ):

        pstn = PSTN({}, {}, {})

        for (src_node, (l, u), tgt_node) in requirements:
            
            if ((src_node, tgt_node) in pstn.requirements
                or (tgt_node, src_node) in pstn.requirements
            ):
                raise ValueError("Multiple requirements defined between the same two nodes.")

            if u < 0:
                raise ValueError("The upper bound of the requirement must be positive.")

            if l > u:
                raise ValueError("The lower bound of the requirement is larger than its upper bound.")

            pstn.requirements[(src_node, tgt_node)] = (l, u)
            pstn.reverse_graph.setdefault(tgt_node, set()).add(src_node)

        for (src_node, l_u_or_prob_distr, tgt_node) in activities:

            if ((src_node, tgt_node) in pstn.activities
                or (tgt_node, src_node) in pstn.activities
            ):
                raise ValueError("Multiple activities defined between the same two nodes.")

            if isinstance(l_u_or_prob_distr, ProbDistr):
                prob_distr = l_u_or_prob_distr

                # TODO check that the lower bound of the activity's duration should not be <(=?) 0.
                # Or maybe that should be done ProbDistrInfo ?
                if prob_distr.definition_domain[0] < 0:     # FIXME: < or <= ?
                    raise ValueError(("The lower bound of the activity's duration cannot be 0 or lower ",
                                      "(the lower bound of the probability distribution's support is negative)."))

                pstn.activities[(src_node, tgt_node)] = prob_distr

            else:
                l, u = l_u_or_prob_distr

                if l > u:
                    raise ValueError("The lower bound of the activity's duration is larger than its upper bound.")

                if l < 0:   # FIXME <= 0 ? and use an epsilon instead of 0 for minimum duration ?
                    raise ValueError("The lower bound of the activity's duration cannot be 0 or lower.")

                pstn.activities[(src_node, tgt_node)] = (l, u)

            pstn.reverse_graph.setdefault(tgt_node, set()).add(src_node)

        return pstn
    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def determine_relevant_probabilistic_activities(
    pstn: PSTN,
    requirement: Tuple[Node, Node],
) -> Tuple[Tuple[Node, Node],...]:

    # To deal with multiple (local) chance constraints, we need to be able to determine
    # the probabilistic activities that are relevant to (i.e. influence the satisfaction of)
    # "requirements" links in the PSTN, on which chance constraints are defined.
    # Insights on that are given in chapter 8.2. of A.J. Wang's thesis.
    #
    # The idea is that we need to collect the probabilistic activities that may
    # influence the satisfaction of the requirement. This corresponds to
    # "non-shared" probabilistic durations that may precede the execution of
    # the source and target nodes of the requirement. By "non-shared" we mean
    # probabilistic activities that do NOT appear both in the threads leading
    # to source and target nodes of the requirement. Indeed, the satisfaction
    # of the requirement corresponds to the difference of execution "dates" of
    # the target and source node fitting into the "window" defined by the requirement.
    # This means that probabilistic activities that influence the execution instant
    # of both the source and target node of the requirement "cancel out".
    # As such, we do not consider as relevant the probabilistic activities
    # that happen before "fork points" common to threads leading to the
    # source and target nodes of the requirement.
    # 
    # This intuition follows from the insights and example given by Wang in chapter 8.2.
    #
    # Thus, we propose this simple (albeit brutal) algorithm:
    # 
    # - Collect the probabilistic activities encountered while performing
    #   "backwards" DFS along activities, starting at the source node of the requirement.
    # 
    # - Do the same starting at the target node of the requirement.
    #
    # - Take the symmetric difference of the collected probabilistic activities.

    src_node, tgt_node = requirement

    if (src_node, tgt_node) not in pstn.requirements:
        raise ValueError("Requirement is not defined.")

    src_node_dfs_stack: List[Node] = [src_node]
    src_node_dfs_visited_nodes: Set[Node] = set()

    probabilistic_activities_collected_from_src_node: Set[Tuple[Node, Node]] = set()

    tgt_node_dfs_stack: List[Node] = [tgt_node]
    tgt_node_dfs_visited_nodes: Set[Node] = set()

    probabilistic_activities_collected_from_tgt_node: Set[Tuple[Node, Node]] = set()

    while src_node_dfs_stack:

        cur_node = src_node_dfs_stack.pop()

        if cur_node == tgt_node:
            continue

        if cur_node in src_node_dfs_visited_nodes:
            continue

        src_node_dfs_visited_nodes.add(cur_node)

        if cur_node not in pstn.reverse_graph:
            continue

        for node in pstn.reverse_graph[cur_node]:

            if node in src_node_dfs_visited_nodes:
                continue

            if (node, cur_node) in pstn.activities:
                src_node_dfs_stack.append(node)

                if isinstance(pstn.activities[(node, cur_node)], ProbDistr):
                    probabilistic_activities_collected_from_src_node.add((node, cur_node))

    while tgt_node_dfs_stack:

        cur_node = tgt_node_dfs_stack.pop()

        if cur_node == src_node:
            continue

        if cur_node in tgt_node_dfs_visited_nodes:
            continue

        tgt_node_dfs_visited_nodes.add(cur_node)

        if cur_node not in pstn.reverse_graph:
            continue

        for node in pstn.reverse_graph[cur_node]:

            if node in tgt_node_dfs_visited_nodes:
                continue

            if (node, cur_node) in pstn.activities:
                tgt_node_dfs_stack.append(node)

                if isinstance(pstn.activities[(node, cur_node)], ProbDistr):
                    probabilistic_activities_collected_from_tgt_node.add((node, cur_node))

    return tuple(probabilistic_activities_collected_from_tgt_node \
                 .symmetric_difference(probabilistic_activities_collected_from_src_node))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def reformulate_chance_constraints(
    pstn: PSTN,
    chance_constraints: Iterable[ChanceConstraint],
) -> Tuple[ReformulatedChanceConstraint,...]:

    rccs_dict: Dict[Tuple[RAVar,...], float] = {}

    for requirements, risk in chance_constraints:

        ra_vars_set: Set[RAVar] = set()

        for req_src_node, req_tgt_node in requirements:

            rpa = determine_relevant_probabilistic_activities(pstn, (req_src_node, req_tgt_node))
            
            for pa_src_node, pa_tgt_node in rpa:
                ra_vars_set.add(((pa_src_node, pa_tgt_node), False))
                ra_vars_set.add(((pa_src_node, pa_tgt_node), True))
        
        ra_vars_tuple = tuple(ra_vars_set)

        if (ra_vars_tuple not in rccs_dict
            or rccs_dict[ra_vars_tuple] > risk
        ):
            rccs_dict[ra_vars_tuple] = risk

    return tuple(ReformulatedChanceConstraint(ra_vars, risk) for ra_vars, risk in rccs_dict.items())

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class ChanceConstraint(NamedTuple):
    """
    A chance constraint is defined by a set of requirements of the PSTN
    and a "risk" value. "1 - this risk value" describes the minimal probability
    with which we allow the requirements of to hold during execution.
    """
    requirements: Tuple[Tuple[Node, Node],...]
    risk: float

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

RAVar = Tuple[Tuple[Node, Node], bool]
"""
We encode / identify a risk allocation variable using the activity from which it originates
and its kind, i.e. "upper bound" (True) or "lower bound" (False) risk allocation variable.
"""

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class ReformulatedChanceConstraint(NamedTuple):
    """
    A reformulated chance constraint is obtained after... well, reformulating a
    list/set of chance constraints. It consists of a set of probabilistic activities of
    the PSTN (not requirements, like for "raw" chance constraints!), and a risk value.

    The reformulated chance constraint expresses that
    `sum(cdf_i(l_i) + (1 - cdf_i(u_i))) must be <= risk_j`, where `l_i / u_i` are the
    risk allocation variables for the probabilistic activity indexed by `i`, and
    `cdf_i` is its cumulative probability distribution function.

    The `cdf_i` functions and the constraint expression itself are not represented here,
    but are fetched right before building the MIP constraints / program.
    """
    ra_vars: Tuple[RAVar,...]
    risk: float
    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class RALinearConstraint(NamedTuple):
    terms: FrozenSet[Tuple[RAVar, float]]
    constant: float
    sense: bool
    """
    True: <=
    False: >=
    """

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class RAMIPCache(NamedTuple):
    mip_model: mip.Model
    ra_vars_col_idx: Dict[RAVar, int]
    ra_vars_col_idx_reverse: Dict[int, RAVar]
    ra_vars_cdf: Dict[RAVar, mip.LinExpr]
    num_constraints_related_to_rccs: int

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def check_pstn_cc_dc(
    pstn: PSTN,
    reformulated_chance_constraints: Iterable[ReformulatedChanceConstraint],
) -> Tuple[bool,
           Optional[Dict[RAVar, float]],
           Optional[STNU],
           Tuple[RALinearConstraint,...]]:
    """
    Returns:
        - bool: Whether a risk allocation was found (succesfully).

        - Optional[Dict[RAVar, float]]: In case of success, the assignments     \
            to risk allocation variables. None in case of failure.

        - Optional[STNU]: In case of success, the "projection STNU", implied    \
            by the risk allocation, i.e. resulting from the projection of risk  \
            allocation variables' assignments to PSTN probabilistic activities  \
            durations.The assignments to risk allocation variables, in case of  \
            success. None in case of failure.
        
        - Tuple[RALinearConstraint,...]: In case of success, a set of conflict      \
            resolution constraints, as proof of the risk allocation's feasibility.  \
            In case of failure, a set of SRNC conflict resolution. (Thus, two       \
            different meanings of success and failure cases !!)
    """

    # We map individual linear constraints (which we may also
    # call "linear disjuncts") to so-called boolean "search variables".
    # In other words, each search variable is mapped to an individual linear constraint,
    # and disjunctive linear constraints are encoded as tuples of search variables.
    sch_vars_to_linear_disjuncts_mapping: Dict[int, RALinearConstraint] = {}
    same_terms_linear_disjuncts_sch_vars: Dict[FrozenSet[Tuple[RAVar, float]], List[int]] = {}
    sch_lits_clauses: List[Tuple[int,...]] = []

    sat_solver = pycryptosat.Solver()

    cached_projection_stnu: List[STNU] = []     # Will contain at most 1 element !!
    cached_ra_mip_info: List[RAMIPCache] = []   # Will contain at most 1 element !! 
#    _n = 0
    while True:
#        if _n > 30:
#            return (False,
#                    None,
#                    None,
#                    tuple(sch_vars_to_linear_disjuncts_mapping[abs(sch_lit)]
#                          for clause in sch_lits_clauses
#                          for sch_lit in clause))
        result = sat_solver.solve()

        exists: bool = result[0]
        sch_vars_assignments: Tuple[Optional[bool],...] = result[1] 
        # NOTE: the 1st element (index: 0) does not correspond to a search variable.
        # This is because "0" is special in pycryptosat.

        print(sch_vars_assignments)

        if not exists:
            return (False,
                    None,
                    None,
                    tuple(sch_vars_to_linear_disjuncts_mapping[abs(sch_lit)]
                          for clause in sch_lits_clauses
                          for sch_lit in clause))

        selected_linear_constraints: Set[RALinearConstraint] = { sch_vars_to_linear_disjuncts_mapping[sch_var]
                                                                 for sch_var, sch_var_val in enumerate(sch_vars_assignments)
                                                                 if sch_var_val is True }

        non_disjunctive_srnc_conflict_resolution_constraints: Set[RALinearConstraint] = set()
#        _n += 1
        while True:
#            _n += 1
            current_linear_constraints = \
                selected_linear_constraints \
                    .union(non_disjunctive_srnc_conflict_resolution_constraints)
#                selected_linear_constraints.union({ sch_vars_to_linear_disjuncts_mapping[abs(sch_lit)]
#                                                    for sch_lit in unary_conflict_resolution_conflicts})
            print("--")
            print(len(selected_linear_constraints))
            print(len(non_disjunctive_srnc_conflict_resolution_constraints))

            (feasible,
            ra_vars_assignments,
            ra_conflict_resolution_constraint) = build_and_solve_ra_mip(pstn,
                                                                        reformulated_chance_constraints,
                                                                        current_linear_constraints,
                                                                        cached_ra_mip_info)

            if not feasible:

                if sch_vars_assignments is not None: 
                    sat_solver.add_clause([-sch_var if sch_var_val else sch_var
                                           for sch_var, sch_var_val in enumerate(sch_vars_assignments)
                                           if sch_var_val is not None])

                assert ra_conflict_resolution_constraint is not None
                learn_constraint(sat_solver,
                                 ra_conflict_resolution_constraint,
                                 sch_vars_to_linear_disjuncts_mapping,
                                 same_terms_linear_disjuncts_sch_vars,
                                 sch_lits_clauses)
                break

            assert ra_vars_assignments is not None

            dc, srnc_conflict_resolution_constraint = \
                project_ra_onto_pstn_and_check_projection_stnu_dc(pstn,
                                                                  ra_vars_assignments,
                                                                  cached_projection_stnu)

            assert len(cached_projection_stnu) > 0

            if dc:
                return (True,
                        ra_vars_assignments,
                        cached_projection_stnu[0],
                        tuple(current_linear_constraints))

            assert srnc_conflict_resolution_constraint is not None

            learn_constraint(sat_solver,
                             srnc_conflict_resolution_constraint,
                             sch_vars_to_linear_disjuncts_mapping,
                             same_terms_linear_disjuncts_sch_vars,
                             sch_lits_clauses)
                
            if len(srnc_conflict_resolution_constraint) == 0:
                assert False

            elif len(srnc_conflict_resolution_constraint) == 1:
            
                sch_lit, = sch_lits_clauses[-1]
                assert sch_lit > 0
                non_disjunctive_srnc_conflict_resolution_constraints \
                    .add(sch_vars_to_linear_disjuncts_mapping[abs(sch_lit)])

            else:
                break

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def learn_constraint(
    sat_solver: pycryptosat.Solver,
    constraint: Iterable[RALinearConstraint],
    sch_vars_to_linear_disjuncts_mapping: Dict[int, RALinearConstraint],
    same_terms_linear_disjuncts_sch_vars: Dict[FrozenSet[Tuple[RAVar, float]], List[int]],
    sch_lits_clauses: List[Tuple[int,...]],
) -> bool:

    clause: List[int] = []

    for lin_disj in constraint:

        if lin_disj.terms not in same_terms_linear_disjuncts_sch_vars:

            sch_var = len(sch_vars_to_linear_disjuncts_mapping)+1
            sch_vars_to_linear_disjuncts_mapping[sch_var] = lin_disj
            same_terms_linear_disjuncts_sch_vars[lin_disj.terms] = [sch_var]
            clause.append(sch_var)

        else:

            for sch_var2 in same_terms_linear_disjuncts_sch_vars[lin_disj.terms].copy():
                lin_disj2 = sch_vars_to_linear_disjuncts_mapping[sch_var2]

                match lin_disj.sense, lin_disj2.sense:

                    case True, True:
                        if lin_disj.constant == lin_disj2.constant:
                            clause.append(sch_var2)
                        else:
                            sch_var = len(sch_vars_to_linear_disjuncts_mapping)+1
                            sch_vars_to_linear_disjuncts_mapping[sch_var] = lin_disj
                            same_terms_linear_disjuncts_sch_vars[lin_disj.terms].append(sch_var)
                            clause.append(sch_var)                                
                            
                            if lin_disj.constant < lin_disj2.constant:
#                                sat_solver.add_xor_clause((-sch_var, sch_var2))
                                sat_solver.add_clause((-sch_var, sch_var2))
                            elif lin_disj.constant > lin_disj2.constant:
#                                sat_solver.add_xor_clause((-sch_var2, sch_var))
                                sat_solver.add_clause((-sch_var2, sch_var))
                            else:
                                assert False

                    case False, False:
                        if lin_disj.constant == lin_disj2.constant:
                            clause.append(sch_var2)
                        else:
                            sch_var = len(sch_vars_to_linear_disjuncts_mapping)+1
                            sch_vars_to_linear_disjuncts_mapping[sch_var] = lin_disj
                            same_terms_linear_disjuncts_sch_vars[lin_disj.terms].append(sch_var)
                            clause.append(sch_var)                                

                            if lin_disj.constant <= lin_disj2.constant:
#                                sat_solver.add_xor_clause((-sch_var2, sch_var))
                                sat_solver.add_clause((-sch_var2, sch_var))
                            elif lin_disj.constant >= lin_disj2.constant:
#                                sat_solver.add_xor_clause((-sch_var, sch_var2))
                                sat_solver.add_clause((-sch_var, sch_var2))
                            else:
                                assert False

                    case True, False:
                        sch_var = len(sch_vars_to_linear_disjuncts_mapping)+1
                        sch_vars_to_linear_disjuncts_mapping[sch_var] = lin_disj
                        same_terms_linear_disjuncts_sch_vars[lin_disj.terms].append(sch_var)
                        clause.append(sch_var)                                

                        if lin_disj.constant >= lin_disj2.constant:
#                            sat_solver.add_xor_clause((sch_var, sch_var2))
                            sat_solver.add_clause((sch_var, sch_var2))
                            sat_solver.add_clause((-sch_var, -sch_var2))
                        else:
                            return False

                    case True, False:
                        sch_var = len(sch_vars_to_linear_disjuncts_mapping)+1
                        sch_vars_to_linear_disjuncts_mapping[sch_var] = lin_disj
                        same_terms_linear_disjuncts_sch_vars[lin_disj.terms].append(sch_var)
                        clause.append(sch_var)                                

                        if lin_disj.constant <= lin_disj2.constant:
#                            sat_solver.add_xor_clause((sch_var, sch_var2))
                            sat_solver.add_clause((sch_var, sch_var2))
                            sat_solver.add_clause((-sch_var, -sch_var2))
                        else:
                            return False

    if len(clause) == 0:
        return False
        
    sat_solver.add_clause(clause)
    sch_lits_clauses.append(tuple(clause))
    return True

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def build_and_solve_ra_mip(
    pstn: PSTN,
    reformulated_chance_constraints: Iterable[ReformulatedChanceConstraint],
    other_constraints: Iterable[RALinearConstraint],
    cached_ra_mip_info: List[RAMIPCache],
) -> Tuple[bool,
           Optional[Dict[RAVar, float]],
           Optional[Set[RALinearConstraint]]]:

    if len(cached_ra_mip_info) == 0:

        mip_base_model = mip.Model()
        mip_base_model.verbose = 0
        mip_base_model.round_int_vars = False
        
        ra_vars_col_idx: Dict[RAVar, int] = {}
        ra_vars_col_idx_reverse: Dict[int, RAVar] = {}
        ra_vars_cdf: Dict[RAVar, mip.LinExpr] = {}

        for rcc in reformulated_chance_constraints:
            # FIXME: avoid duplication when a ra_var can appear in multiple chance constraints
            for ra_var in rcc.ra_vars:

                prob_distr = pstn.activities[ra_var[0]]
                assert isinstance(prob_distr, ProbDistr)
                cdf = prob_distr.cumul_distr_func
                discr_points = prob_distr.discr_points
                num_discr_points = len(discr_points)

                if ra_var not in ra_vars_col_idx:

                    v = mip_base_model.add_var()

                    ws = [mip_base_model.add_var(lb=0, ub=1, var_type=mip.CONTINUOUS) for _ in range(num_discr_points)]      # type: ignore
                    mip_base_model.add_constr(mip.xsum(ws) == 1)
                    mip_base_model.add_constr(v == mip.xsum(discr_points[k] * ws[k] for k in range(num_discr_points)))       # type: ignore
                    mip_base_model.add_sos([(ws[k], discr_points[k]) for k in range(num_discr_points)], 2)                   # type: ignore

                    F_v = mip.xsum(cdf(discr_points[k]) * ws[k] for k in range(num_discr_points))               # type: ignore

                    ra_vars_col_idx[ra_var] = v.idx
                    ra_vars_col_idx_reverse[v.idx] = ra_var
                    ra_vars_cdf[ra_var] = F_v

            for ra_var, idx in ra_vars_col_idx.items():
                if ra_var[1]:
                    mip_base_model.add_constr(mip_base_model.vars[ra_vars_col_idx[(ra_var[0], not ra_var[1])]] <= mip_base_model.vars[idx])    # type: ignore

    #        m.add_constr(rcc.risk >= mip.xsum((1-ra_vars_cdf[ra_var]) if ra_var[1] else ra_vars_cdf[ra_var] # type: ignore
    #                                          for ra_var in rcc.ra_vars))                                   # type: ignore

            expr = mip.xsum((1-ra_vars_cdf[ra_var]) if ra_var[1] else ra_vars_cdf[ra_var] # type: ignore
                            for ra_var in rcc.ra_vars)

            mip_base_model.add_constr(expr <= rcc.risk)  # type: ignore
            mip_base_model.add_constr(0 <= expr)         # type: ignore
        
        num_constraints_related_to_rccs: int = mip_base_model.num_rows
        cached_ra_mip_info.append(RAMIPCache(mip_base_model,
                                             ra_vars_col_idx,
                                             ra_vars_col_idx_reverse,
                                             ra_vars_cdf,
                                             num_constraints_related_to_rccs))
    else:

        (mip_base_model,
         ra_vars_col_idx,
         ra_vars_col_idx_reverse,
         ra_vars_cdf,
         num_constraints_related_to_rccs) = cached_ra_mip_info[0]
    
    mip_model = mip_base_model.copy()

    for terms, constant, sense in other_constraints:

        if sense:
            mip_model.add_constr(mip.xsum([coeff*mip_model.vars[ra_vars_col_idx[ra_var]]    # type: ignore
                                           for ra_var, coeff in terms]) <= constant)   # type: ignore
        else:
            mip_model.add_constr(mip.xsum([coeff*mip_model.vars[ra_vars_col_idx[ra_var]]    # type: ignore
                                           for ra_var, coeff in terms]) >= constant)   # type: ignore
    print(mip_model.num_rows)
    mip_model.optimize()

    if (mip_model.status == mip.OptimizationStatus.FEASIBLE
        or mip_model.status == mip.OptimizationStatus.OPTIMAL
    ):
#        print(expr.x)   # type: ignore
        return (True,
                { ra_var: mip_model.vars[idx].x for ra_var, idx in ra_vars_col_idx.items() },   # type: ignore
                None)
    
    if mip_model.status == mip.OptimizationStatus.INFEASIBLE:
    
        mip_aux_model = mip_base_model.copy()

        mip_aux_model.emphasis = mip.SearchEmphasis.FEASIBILITY
        mip_aux_model.preprocess = 1

        mip_aux_model.objective = 0     # type: ignore

        n = num_constraints_related_to_rccs
        for constr in mip_model.constrs[num_constraints_related_to_rccs:]: # type: ignore

            n+=1
            mip_aux_model.add_constr(constr.expr)
            mip_aux_model.optimize()

            if (mip_aux_model.status == mip.OptimizationStatus.INFEASIBLE
                or mip_aux_model.status == mip.OptimizationStatus.INT_INFEASIBLE
            ):
                break

        for constr in reversed(mip_model.constrs[num_constraints_related_to_rccs:n]):   # type: ignore
            print(constr.expr)
#            mip_aux_model.constrs.remove([constr])
            mip_aux_model.remove(constr)
            mip_aux_model.optimize()

            if (mip_aux_model.status == mip.OptimizationStatus.FEASIBLE
                or mip_aux_model.status == mip.OptimizationStatus.OPTIMAL
            ):
                mip_aux_model.add_constr(constr.expr)     

        iis = mip_aux_model.constrs[num_constraints_related_to_rccs:]
        print(iis)
        return (False, 
                None,
                { RALinearConstraint(frozenset((ra_vars_col_idx_reverse[v.idx], vx)            # type: ignore
                                               for v, vx in cstr.expr.expr.items()),           # type: ignore
                                     cstr.rhs,
#                                     cstr.rhs+0.01 if cstr.expr.sense=="<" else cstr.rhs-0.01,
                                     False if cstr.expr.sense=="<" else True)       # FIXME
                  for cstr in iis })                                                            # type: ignore

    assert False

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def project_ra_onto_pstn_and_check_projection_stnu_dc(
    pstn: PSTN,
    ra_vars_assignments: Dict[RAVar, float],
    cached_projection_stnu: List[STNU],
) -> Tuple[bool, Optional[Set[RALinearConstraint]]]:

    if len(cached_projection_stnu) == 0:

        ordinary_as_dict = { (src_node, tgt_node): l_u_or_prob_distr
                        for (src_node, tgt_node), l_u_or_prob_distr in pstn.activities.items()
                        if not isinstance(l_u_or_prob_distr, ProbDistr) }

        # overwrite (non probabilistic) activities with requirements defined on top of them !
        for (src_node, tgt_node), (l, u) in pstn.requirements.items():
            ordinary_as_dict[(src_node, tgt_node)] = (l, u)

        contingent = [(src_node,
                       (ra_vars_assignments[((src_node, tgt_node), False)],
                        ra_vars_assignments[((src_node, tgt_node), True)]),
                       tgt_node)
                      for (src_node, tgt_node), l_u_or_prob_distr in pstn.activities.items()
                      if isinstance(l_u_or_prob_distr, ProbDistr)]

        cached_projection_stnu.append(STNU.from_links([(src_node, w_or_l_u, tgt_node)
                                                      for (src_node, tgt_node), w_or_l_u in ordinary_as_dict.items()],
                                                      contingent))

    else:

        for (src_node, tgt_node), l_u_or_prob_distr in pstn.activities.items():

            if isinstance(l_u_or_prob_distr, ProbDistr):

                l = ra_vars_assignments[((src_node, tgt_node), False)]
                u = ra_vars_assignments[((src_node, tgt_node), True)]

                lbl_src_tgt, lbl_tgt_src = (cached_projection_stnu[0].labeled_weights[(src_node, tgt_node)][0],
                                            cached_projection_stnu[0].labeled_weights[(tgt_node, src_node)][0])

                cached_projection_stnu[0].labeled_weights[(src_node, tgt_node)] = (lbl_src_tgt, l)
                cached_projection_stnu[0].labeled_weights[(tgt_node, src_node)] = (lbl_tgt_src, -u)

    dc, srnc_conflict = check_stnu_dc(cached_projection_stnu[0])
    # FIXME !!!!!!!!!!! ('E', 'D') edge twice (successively) in one of the paths !!!!!

    # To obtain conflict resolutions constraints (negated linear constraints),
    # replace all non contingent edges with their weights and
    # contingent edges with their risk allocation variables
    # (and change the sense of the inequality) of the linear constraint.
    srnc_conflict_resolution_constraint: Optional[Set[RALinearConstraint]] = None

    for path in srnc_conflict:

        terms: Dict[RAVar, float] = {}
        constant = 0

        for (src_node, tgt_node) in path:

            if ((src_node, tgt_node) in pstn.activities
                and isinstance(pstn.activities[(src_node, tgt_node)], ProbDistr)
            ):
                if ((src_node, tgt_node), False) in terms:
                    terms[((src_node, tgt_node), False)] += 1
                else:
                    terms[((src_node, tgt_node), False)] = 1

            elif ((tgt_node, src_node) in pstn.activities
                and isinstance(pstn.activities[(tgt_node, src_node)], ProbDistr)
            ):
                if ((tgt_node, src_node), True) in terms:
                    terms[((tgt_node, src_node), True)] += -1
                else:
                    terms[((tgt_node, src_node), True)] = -1

            else:
                constant += -cached_projection_stnu[0].labeled_weights[(src_node, tgt_node)][1]

        if terms:
            if srnc_conflict_resolution_constraint is None:
                srnc_conflict_resolution_constraint = set()
            srnc_conflict_resolution_constraint.add(RALinearConstraint(frozenset(terms.items()),
                                                                       constant,
                                                                       False))

    return dc, srnc_conflict_resolution_constraint

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
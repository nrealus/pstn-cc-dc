"""
The functions defined in this module allow us to check dynamic
controllability of a STNU and, if it is uncontrollable, to extract
"conflicts" which provide a reason to why it is dynamically uncontrollable.

The algorithm for dynamic controllability checking and conflicts
extraction is based on [1], which in turn is based on Morris' O(n^3)
dynamic controllability checking algorithm [2].

[1]: N.Bhargava et al. *Faster Conflict Generation for Dynamic Controllability* (2017)
[2]: P.Morris *Dynamic Controllability and Dispatchability Relationships* (2014)
"""

from __future__ import annotations

#################################################################################

from heapdict import heapdict
from typing import Dict, List, NamedTuple, Optional, Set, Tuple
from enum import Enum

from common import Node, Weight

#################################################################################

class Label(NamedTuple):
    """
    A type representing the label of an edge in an STNU.
    Labels with `kind = Kind.NONE` must have `node = None`.
    Similarly, nodes with `kind = Kind.LOWERCASE` or
    `kind = Kind.UPPERCASE` must have `node != None`
    """

    class Kind(Enum):
        NONE = 0
        UPPERCASE = 1
        LOWERCASE = 2

    kind: Kind
    node: Optional[Node]

    @classmethod
    def nocase(cls):
        return Label(Label.Kind.NONE, None)

    @classmethod
    def lowercase(cls,
        node: Node
    ):
        return Label(Label.Kind.LOWERCASE, node)

    @classmethod
    def uppercase(cls,
        node: Node
    ):
        return Label(Label.Kind.UPPERCASE, node)

Edge = Tuple[Node, Node, Label, Weight]
"""
Representation of an STNU edge as a tuple:

(source node, target node, label, weight)
"""

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class STNU(NamedTuple):

    reverse_adjacency: Dict[Node, Set[Node]]
    """
    The reverse adjacency set/list representation.
    
    The key of the dictionary represents the target node of an edge in the STNU.
    The value of the dictionary corresponds to the set of source nodes of STNU
    edges incoming to the target node.
    """

    labeled_weights: Dict[Tuple[Node, Node], Tuple[Label, Weight]]
    """
    Stores the label and weight of an STNU edge (source node, target node).
    """

    contingent_links: Set[Tuple[Node, Node]]
    """
    Contains the contingent links (source node, target node) of the STNU.
    """

    @classmethod
    def from_links(cls,
        ordinary: List[Tuple[Node, Weight | Tuple[Weight, Weight], Node]],
        contingent: List[Tuple[Node, Tuple[Weight, Weight], Node]],
    ):

        stnu = STNU({}, {}, set())
        
        for (src_node, (l, u), tgt_node) in contingent:
            
            if ((tgt_node in stnu.reverse_adjacency and src_node in stnu.reverse_adjacency[tgt_node])
                or (src_node in stnu.reverse_adjacency and tgt_node in stnu.reverse_adjacency[src_node])
            ):
                raise ValueError("Two links defined between the same two nodes.")

            stnu.reverse_adjacency.setdefault(tgt_node, set()).add(src_node)
            stnu.labeled_weights[(src_node, tgt_node)] = (Label.lowercase(tgt_node), l)

            stnu.reverse_adjacency.setdefault(src_node, set()).add(tgt_node)
            stnu.labeled_weights[(tgt_node, src_node)] = (Label.uppercase(tgt_node), -u)

            stnu.contingent_links.add((src_node, tgt_node))

        for (src_node, w_or_l_u, tgt_node) in ordinary:

            if ((tgt_node in stnu.reverse_adjacency and src_node in stnu.reverse_adjacency[tgt_node])
                or (src_node in stnu.reverse_adjacency and tgt_node in stnu.reverse_adjacency[src_node])
            ):
                raise ValueError("Two links defined between the same two nodes.")

            if isinstance(w_or_l_u, tuple):
                l, u = w_or_l_u

                stnu.reverse_adjacency.setdefault(tgt_node, set()).add(src_node)
                stnu.labeled_weights[(src_node, tgt_node)] = (Label.nocase(), u)

                stnu.reverse_adjacency.setdefault(src_node, set()).add(tgt_node)
                stnu.labeled_weights[(tgt_node, src_node)] = (Label.nocase(), -l)

            else:
                w = w_or_l_u

                stnu.reverse_adjacency.setdefault(tgt_node, set()).add(src_node)
                stnu.labeled_weights[(src_node, tgt_node)] = (Label.nocase(), w)

        return stnu

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def convert_stnu_to_normal_form(
    stnu: STNU,
) -> Tuple[STNU, Dict[Node, Tuple[Node, Node]]]:
    """
    Converts the given STNU into a new STNU in "normal form".

    The normal form STNU of an original STNU, is a STNU where the original
    contingent links have been split into two links: one ordinary (or
    requirement) link, and one contingent link. Both of these links
    are connected using a newly added "helper" node. The length of the
    newly created ordinary link is fixed (lower bound = -upper bound),
    and the lower bound of the duration of the new contingent link is 0.

    Example: The contingent link A ==[l,u]==> C is transformed into
    these two links: A --[l,l]--> A'' ==[0,u-l]==> C.

    Returns:
        - A new STNU, corresponding to the "normal form" of the original STNU.
        
        - A mapping from the helper nodes in the new STNU to the
        original STNU contingent links that these nodes helped to split.

    """
    
    normal_form_stnu = STNU({ node : stnu.reverse_adjacency[node].copy()
                              for node in stnu.reverse_adjacency },
                            stnu.labeled_weights.copy(),
                            stnu.contingent_links.copy())

    normal_form_stnu_helper_nodes: Dict[Node, Tuple[Node, Node]] = {}

    helpers_counter: int = 0

    for src_node, tgt_node in stnu.contingent_links:

        weight_src_to_tgt = stnu.labeled_weights[(src_node, tgt_node)][1]
        weight_tgt_to_src = -stnu.labeled_weights[(tgt_node, src_node)][1]

        normal_form_stnu.reverse_adjacency[tgt_node].remove(src_node)
        normal_form_stnu.reverse_adjacency[src_node].remove(tgt_node)

        helpers_counter += 1
        helper_node = Node(src_node+"''"+str(helpers_counter))

        normal_form_stnu.reverse_adjacency[helper_node] = set()

        normal_form_stnu.reverse_adjacency[helper_node].add(src_node)
        normal_form_stnu.labeled_weights[(src_node, helper_node)] = (Label.nocase(), weight_src_to_tgt)

        normal_form_stnu.reverse_adjacency[src_node].add(helper_node)
        normal_form_stnu.labeled_weights[(helper_node, src_node)] = (Label.nocase(), -weight_src_to_tgt)

        normal_form_stnu.reverse_adjacency[tgt_node].add(helper_node)
        normal_form_stnu.labeled_weights[(helper_node, tgt_node)] = (Label.lowercase(tgt_node), Weight(0))

        normal_form_stnu.reverse_adjacency[helper_node].add(tgt_node)
        normal_form_stnu.labeled_weights[(tgt_node, helper_node)] = (Label.uppercase(tgt_node), weight_src_to_tgt
                                                                                                - weight_tgt_to_src)

        normal_form_stnu.contingent_links.add((helper_node, tgt_node))

        normal_form_stnu_helper_nodes[helper_node] = (src_node, tgt_node)

    return (normal_form_stnu, normal_form_stnu_helper_nodes)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def check_stnu_dc(
    stnu: STNU,
) -> Tuple[bool, List[List[Tuple[Node, Node]]]]:
    """
    Checks whether the given STNU is dynamically controllable (DC),
    and if not, provides a reason to why it isn't, in the form of a
    set of "conflicts". These take the form of a semi-reducible
    negative cycle (SRNC) combined with a (possibly empty) set
    of "extension paths" of lowercase edges appearing in the SRNC.

    Implements the top level of Morris' O(n^3) DC-checking algorithm,
    with slight modifications as in Bhargava et al. (Algorithm 1).

    Returns:
        - Whether the STNU is DC.

        - If the STNU is not DC, a set of "conflicts".
    """
    
    # The inverse labeled distance graph of the normal
    # form STNU may be modified by `dc_disjktra`. Newly added edges
    normal_form_stnu, normal_form_stnu_helper_nodes = convert_stnu_to_normal_form(stnu)

    # The booleans in `negative_nodes` indicate whether a negative node
    # (i.e. a node with an incoming edge with a negative weight) has
    # already been (successfully) processed during `dc_dijkstra`,
    # (i.e. no SRNC was found by starting backwards from that node).
    negative_nodes: Dict[Node, bool] = { tgt_node: False
                                         for tgt_node, src_nodes in normal_form_stnu.reverse_adjacency.items()
                                         for src_node in src_nodes
                                         if normal_form_stnu.labeled_weights[(src_node, tgt_node)][1] < 0 }

    # The edges newly added to the (normal form) STNU during runs of `dc_dijkstra`.
    # Note that updated already existing edges won't be recorded.
    novel_edges: Set[Tuple[Node, Node]] = set()

    # `distances_to` stores the shortest distance from a source node
    # (key of the inner dictionary) to a target node (key of the outer 
    # dictionary) as well as the first edge of the shortest path between them.
    distances_to: Dict[Node, Dict[Node, Tuple[Weight, Optional[Edge]]]] = {}

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    for start_node, was_processed in negative_nodes.items():
        
        if was_processed:
            continue

        dc, raw_paths, _ = dc_dijkstra(start_node,
                                       normal_form_stnu.reverse_adjacency,
                                       normal_form_stnu.labeled_weights,
                                       novel_edges,
                                       negative_nodes,
                                       distances_to,
                                       [])

        if not dc:

            # The paths in `raw_paths` form a semi-reducible negative cycle.
            # These paths result from the application of a series of
            # "reduction rules" (see STNU theory) on lowercase edges.
            #
            # The "conflicts" are the full semi-reducible negative cycle
            # itself (or rather the path of edges composing it), as well
            # as these paths / sets of edges described above,
            # without their first edge.

            # But first, before we extract the conflicts, we need to
            # process these raw paths. This is done in two steps, for
            # each raw path:
            #
            # - Decomposing those of their edges that were
            #   not present in the initial STNU in normal form.
            # 
            # - Converting the edges from the initial normal form STNU
            #   to edges of the original, not necessarily normal form STNU.

            processed_paths: List[List[Tuple[Node, Node]]] = []

            for raw_path in raw_paths:
                processed_paths.append(process_raw_path(raw_path,
                                                        novel_edges,
                                                        distances_to,
                                                        normal_form_stnu_helper_nodes))
            
            conflicts: List[List[Tuple[Node, Node]]] = [[]]

            for processed_path in processed_paths:
                conflicts[0].extend(processed_path)

            if len(processed_paths) > 1:
                for processed_path in processed_paths:
                    if len(processed_path) > 1:
                        conflicts.append(processed_path[1:])

#            print("----")
#            print(raw_paths)
#
#            print("--")
#            print(processed_paths)
#
#            print("--")
#            print(conflicts)
#            
#            ws = []
#            for conflict in conflicts:
#                w = 0
#                for src_node, tgt_node in conflict:
#                    w += stnu.labeled_weights[(src_node, tgt_node)][1]
#                ws.append(w)
#            print(ws)
#            
#            print("----")

            return False, conflicts

    return True, []

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def dc_dijkstra(
    start_node: Node,
    normal_form_stnu_reverse_adjacency: Dict[Node, Set[Node]],
    normal_form_stnu_labeled_weights: Dict[Tuple[Node, Node], Tuple[Label, Weight]],
    novel_edges: Set[Tuple[Node, Node]],
    negative_nodes: Dict[Node, bool],
    distances_to: Dict[Node, Dict[Node, Tuple[Weight, Optional[Edge]]]],
    call_stack: List[Node],
) -> Tuple[bool,
           List[List[Tuple[Node, Node]]],
           Optional[Node]]:
    """
    The main function of the DC-checking algorithm.

    Implements Bhargava et al.'s Algorithm 2, which is based on
    Morris' DCbackprop procedure.
    """

    call_stack.append(start_node)

    if start_node in call_stack[:-1]:
        call_stack.pop()
        return False, [], start_node

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    # Items type: Tuple[Node, Label]
    # Priority type: Weight
    priority_queue = heapdict()

    distances_to_start_node_from: Dict[Node, Tuple[Weight, Optional[Edge]]] = \
        { start_node : (Weight(0), None) }

    for e_source_node in normal_form_stnu_reverse_adjacency[start_node]:

        e_label, e_weight = normal_form_stnu_labeled_weights[(e_source_node, start_node)]
        edge: Edge = (e_source_node, start_node, e_label, e_weight)

        if e_weight < 0:
            distances_to_start_node_from[e_source_node] = (e_weight, edge)
            priority_queue[(e_source_node, e_label)] = e_weight

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    distances_to[start_node] = distances_to_start_node_from

    while priority_queue:

        (node, label), weight = priority_queue.popitem()

        node: Node
        label: Label
        weight: Weight

        if weight >= 0:
            if (start_node not in normal_form_stnu_reverse_adjacency
                or node not in normal_form_stnu_reverse_adjacency[start_node]
            ):
                novel_edges.add((node, start_node))
            normal_form_stnu_reverse_adjacency.setdefault(start_node, set()).add(node)
            normal_form_stnu_labeled_weights[(node, start_node)] = (Label.nocase(), weight)
            continue

        if node in negative_nodes:

            srnc_found, raw_paths, end_node = dc_dijkstra(node,
                                                          normal_form_stnu_reverse_adjacency,
                                                          normal_form_stnu_labeled_weights,
                                                          novel_edges,
                                                          negative_nodes,
                                                          distances_to,
                                                          call_stack)
            if not srnc_found:

                if end_node is not None:

                    raw_path: List[Tuple[Node, Node]] = build_raw_path(node,
                                                                       start_node,
                                                                       distances_to)
                    raw_paths.append(raw_path)

                if end_node == start_node:
                    end_node = None

                call_stack.pop()
                return False, raw_paths, end_node

        for e_source_node in normal_form_stnu_reverse_adjacency[node]:

            e_label, e_weight = normal_form_stnu_labeled_weights[(e_source_node, node)]
            edge: Edge = (e_source_node, node, e_label, e_weight)

            if e_weight < 0:
                continue

            if (e_label.kind == Label.Kind.LOWERCASE
                and label.kind != Label.Kind.NONE
                and e_label.node == label.node 
            ):
                continue

            w = weight + e_weight
            
            if (e_source_node not in distances_to_start_node_from
                or w < distances_to_start_node_from[e_source_node][0]
            ):
                distances_to_start_node_from[e_source_node] = (w, edge)
                priority_queue[(e_source_node, label)] = w

    negative_nodes[start_node] = True

    call_stack.pop()
    return True, [], None

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def build_raw_path(
    src_node: Node,
    tgt_node: Node,
    distances_to: Dict[Node, Dict[Node, Tuple[Weight, Optional[Edge]]]],
) -> List[Tuple[Node, Node]]:
    """
    Builds the shortest path from `src_node` to `tgt_node`,
    by "hopping" from node to node using `distances_to`.

    Note that the path may go through edges absent from
    the initial normal form STNU. These edges would
    be added during the run of `dc_dijkstra`.
    """

    path = []
    cur_node = src_node

    while True:
        _, edge = distances_to[tgt_node][cur_node]
        assert edge is not None

        next_node = edge[1]

        path.append((cur_node, next_node))

        if next_node == tgt_node:
            break
        cur_node = next_node

    return path

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def process_raw_path(
    raw_path: List[Tuple[Node, Node]],
    novel_edges: Set[Tuple[Node, Node]],
    distances_to: Dict[Node, Dict[Node, Tuple[Weight, Optional[Edge]]]],
    normal_form_stnu_helper_nodes: Dict[Node, Tuple[Node, Node]],
) -> List[Tuple[Node, Node]]:
    """
    Processes a raw path in two passes:

    1. Decomposing the edges that were not present
    in the initial normal form STNU

    2. Converting the edges (that now all match to edges in
    the normal form STNU) to edges of the original, not
    necessarily normal form STNU. In other words, remove
    helper nodes (however it is slightly more complicated
    than just removing them! see comment below.)
    """

    normal_form_edges_path = []

    for src_node, tgt_node in raw_path:

        stack = [(src_node, tgt_node)]

        while stack:
            src_node, tgt_node = stack.pop()

            if (src_node,tgt_node) not in novel_edges:
                normal_form_edges_path.append((src_node, tgt_node))
            
            else:

                ext = []
                _, e = distances_to[tgt_node][src_node]
                assert e is not None

                while True:
                    ext.append((e[0], e[1]))

                    if e[1] == tgt_node:
                        break

                    _, e = distances_to[tgt_node][e[1]]
                    assert e is not None
                
                stack.extend(reversed(ext))

    # To remove helper nodes in the path, the following rules are applied:
    # 
    # - If the first edge of the path starts from a helper node,
    # replace that helper node with the source node of the original
    # contingent link which was split using that helper node.
    # 
    # - Similarly, if the last edge of the path ends with a helper node,
    # replace it with the target node of the original contingent link
    # which was split using that helper node.
    # 
    # - If any (except the last one) edge of the path ends with a helper node,
    # check whether this edge's source node is the same as the next edge's
    # target node (i.e. whether these two edges are "symmetric"). If they
    # are symmetric, we can "unify" these two edges into one, "skipping"
    # the helper node inbetween. If they are symmetric, we extend them both
    # (symmetrically) up to the target node of the original contingent link.

#    print(normal_form_edges_path)

    original_edges_path = []

    # TODO: The following works but... hic sunt dracones...!
    # should come back to this at some point and try to simplify / clarify this.

    n = len(normal_form_edges_path)
    for i, (src_node, tgt_node) in enumerate(normal_form_edges_path):

        stack = [(src_node, tgt_node)]

        if src_node in normal_form_stnu_helper_nodes:

            if i == 0:
                s, t = normal_form_stnu_helper_nodes[src_node]
                if s != tgt_node:
                    original_edges_path.append((s, tgt_node))
                else:
                    original_edges_path.append((t, tgt_node))

        elif tgt_node in normal_form_stnu_helper_nodes:

            if i == n-1:
                s, _ = normal_form_stnu_helper_nodes[tgt_node]
                original_edges_path.append((src_node, s))

            elif i < n-1:
                
                s, t = normal_form_stnu_helper_nodes[tgt_node]
                _, tt = normal_form_edges_path[i+1]
                
                if src_node == tt:
                    if src_node == s:
                        original_edges_path.append((src_node, t))
                        original_edges_path.append((t, src_node))
                    else:
                        original_edges_path.append((src_node, s))
                        original_edges_path.append((s, src_node))
                else:
                    original_edges_path.append((src_node, tt))

                continue

            else:
                assert False

        else:
            original_edges_path.append((src_node, tgt_node))

    return original_edges_path

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
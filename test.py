from __future__ import annotations

#################################################################################

import unittest

#from pstn_and_stnu.common import Node, Weight
from pstn_and_stnu.stnu import STNU, check_stnu_dc
from pstn_and_stnu.pstn import (PSTN, ProbDistr, ChanceConstraint,
                                build_and_solve_ra_mip,
                                determine_relevant_probabilistic_activities,
                                reformulate_chance_constraints,
                                check_pstn_cc_dc)

#################################################################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class TestSTNU(unittest.TestCase):

    def test01_non_dc(self):
        stnu = STNU.from_links([("C", -2 ,"A"),
                                ("A", 3, "B"),
                                ("B", -2, "C")],
                               [])
        results = check_stnu_dc(stnu)
        self.assertFalse(results[0])

    def test02_dc(self):
        stnu = STNU.from_links([("C", -1 ,"A"),
                                ("A", 3, "B"),
                                ("B", -2, "C")],
                               [])
        results = check_stnu_dc(stnu)
        self.assertTrue(results[0])

    def test03_dc(self):
        stnu = STNU.from_links([("A", 2, "C"),
                                ("C", 2, "B")],
                               [("A", (0, 3), "B")])
        results = check_stnu_dc(stnu)
        self.assertTrue(results[0])

    def test04_non_dc(self):
        stnu = STNU.from_links([("E", (0.5, 1), "H"),
                                ("H", (0, 2), "I"),
                                ("G", (0, 3), "I")],
                               [("D", (6, 11), "E"),
                                ("F", (2, 4), "G")])
        results = check_stnu_dc(stnu)
        self.assertFalse(results[0])

    def test05_dc(self):
        stnu = STNU.from_links([("X", (7, 12), "C0"),
                                ("C1", (1, 11), "C0")],
                               [("A0", (1, 3), "C0"),
                                ("A1", (1, 10), "C1")])
        results = check_stnu_dc(stnu)
        self.assertTrue(results[0])

    def test06_dc(self):
        stnu = STNU.from_links([("X", (29, 48), "C0"),
                                ("C1", (1, 8), "C0"),
                                ("C2", (7, 37), "C0")],
                               [("A0", (1, 3), "C0"),
                                ("A1", (1, 10), "C1"),
                                ("A2", (1, 36), "C2")])
        results = check_stnu_dc(stnu)
        self.assertTrue(results[0])

    def test07_non_dc(self):
        stnu = STNU.from_links([("E1", (0, 50), "E3"),
                                ("E2", (40, 45), "E3")],
                               [("E1", (20, 30), "E2")])
        results = check_stnu_dc(stnu)
        self.assertFalse(results[0])

    def test08_non_dc(self):
        stnu = STNU.from_links([("E2", (1, 2), "E3"),
                                ("E1", (1, 10), "E2")],
                               [("E1", (1, 8), "E3")])
        results = check_stnu_dc(stnu)
        self.assertFalse(results[0])

    def test09_dc(self):
        stnu = STNU.from_links([("10", (399, 650), "1"),
                                ("9", (109, 468), "1"),
                                ("7", (29, 128), "1"),
                                ("5", (7, 34), "1"),
                                ("3", (1, 8), "1")],
                               [("0", (1, 3), "1"),
                                ("2", (1, 10), "3"),
                                ("4", (1, 36), "5"),
                                ("6", (1, 130), "7"),
                                ("8", (1, 470), "9")])
        results = check_stnu_dc(stnu)
        self.assertTrue(results[0])

    def test10_non_dc(self):
        stnu = STNU.from_links([("E1", (1, 100), "E2"),
                                ("E2", (0, 100), "E5"),
                                ("E2", (1, 100), "E3"),
                                ("E3", (1.5, 100), "E4"),
                                ("E1", (0, 3.5), "E4")],
                               [("E1", (0.6294, 18.8554), "E5")])
        results = check_stnu_dc(stnu)
        self.assertFalse(results[0])

    def test11_non_dc(self):
        stnu = STNU.from_links([("E4", (-1, 3), "E2"),
                                ("E5", (2, 4), "E2")],
                               [("E1", (0, 2), "E2"),
                                ("E3", (0, 3), "E4")])
        results = check_stnu_dc(stnu)
        self.assertFalse(results[0])

    def test12_non_dc(self):
        stnu = STNU.from_links([("E2", (0, 0), "E4"),
                                ("E3", (0, 0), "E4")],
                               [("E1", (0, 1), "E2"),
                                ("E1", (0, 1), "E3")])
        results = check_stnu_dc(stnu)
        self.assertFalse(results[0])

    def test13_dc(self):
        stnu = STNU.from_links([("E1", (0, 2), "E2"),
                                ("E2", (0, 2), "E3")],
                               [("E1", (0, 4), "E3")])
        results = check_stnu_dc(stnu)
        self.assertTrue(results[0])

    def test14_non_dc(self):
        stnu = STNU.from_links([("E3", (0, 2), "E2")],
                               [("E1", (0, 10), "E2"),
                                ("E1", (0, 8), "E3")])
        results = check_stnu_dc(stnu)
        self.assertFalse(results[0])

    def test15_non_dc(self):
        stnu = STNU.from_links([("E2", (-1, 1000000), "E3"),
                                ("E4", (0, 1000000), "E5"),
                                ("E6", (0, 1000000), "E7"),
                                ("E2", (5, 18), "E7")],
                               [("E1", (3, 1000000), "E2"),
                                ("E3", (1, 5.5), "E4"),
                                ("E5", (10, 14.5), "E6")])
        results = check_stnu_dc(stnu)
        self.assertFalse(results[0])

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class TestPSTN(unittest.TestCase):

    def test01_relevant_probabilistic_activities(self):

        INF = 100000

        pstn = PSTN.from_requirements_and_activities([("A", (0, 15), "L"),
                                                      ("C", (0, 2), "K"),
                                                      ("H", (0, 2), "I"),
                                                      ("G", (0, 3), "I")],
                                                     [("C", (0, INF), "K"),
                                                      ("K", (0.5, 1), "L"),
                                                      ("E", (0.5, 1), "H"),
                                                      ("H", (0, INF), "I"),
                                                      ("I", (0.5, 1), "J"),
                                                      ("J", (0, INF), "K"),
                                                      ("G", (0, INF), "I"),
                                                      ("A", (0, INF), "B"),
                                                      ("A", (0, INF), "D"),
                                                      ("A", (0, INF), "F"),
                                                      ("B", ProbDistr.uniform((10, 14)), "C"),
                                                      ("D", ProbDistr.uniform((8, 11)), "E"),
                                                      ("F", ProbDistr.uniform((2, 4)), "G")])
  
        relevant_probabilistic_activities = determine_relevant_probabilistic_activities(pstn, ("H", "I"))
        self.assertSetEqual(set(relevant_probabilistic_activities), set([("D", "E"),
                                                                         ("F", "G")]))

        relevant_probabilistic_activities = determine_relevant_probabilistic_activities(pstn, ("G", "I"))
        self.assertSetEqual(set(relevant_probabilistic_activities), set([("D", "E"),
                                                                         ("F", "G")]))

        relevant_probabilistic_activities = determine_relevant_probabilistic_activities(pstn, ("C", "K"))
        self.assertSetEqual(set(relevant_probabilistic_activities), set([("B", "C"),
                                                                         ("D", "E"),
                                                                         ("F", "G")]))

        relevant_probabilistic_activities = determine_relevant_probabilistic_activities(pstn, ("A", "L"))
        self.assertSetEqual(set(relevant_probabilistic_activities), set([("B", "C"),
                                                                         ("D", "E"),
                                                                         ("F", "G")]))

    def test02_relevant_probabilistic_activities(self):

        INF = 100000

        pstn = PSTN.from_requirements_and_activities([("A", (0, 15), "L"),
                                                      ("C", (0, 2), "K"),
                                                      ("H", (0, 2), "I"),
                                                      ("G", (0, 3), "I")],
                                                     [("C", (0, INF), "K"),
                                                      ("K", (0.5, 1), "L"),
                                                      ("E", (0.5, 1), "H"),
                                                      ("H", (0, INF), "I"),
                                                      ("I", (0.5, 1), "J"),
                                                      ("J", (0, INF), "K"),
                                                      ("G", (0, INF), "I"),
                                                      ("A", (0, INF), "B"),
                                                      ("A", (0, INF), "D"),
                                                      ("A", (0, INF), "F"),
                                                      ("B", ProbDistr.uniform((10, 14)), "C"),
                                                      ("D", ProbDistr.uniform((8, 11)), "E"),
                                                      ("F", ProbDistr.uniform((2, 4)), "G"),
                                                      ("X", ProbDistr.uniform((2, 4)), "A")])
  
        relevant_probabilistic_activities = determine_relevant_probabilistic_activities(pstn, ("H", "I"))
        self.assertSetEqual(set(relevant_probabilistic_activities), set([("D", "E"),
                                                                         ("F", "G")]))

        relevant_probabilistic_activities = determine_relevant_probabilistic_activities(pstn, ("G", "I"))
        self.assertSetEqual(set(relevant_probabilistic_activities), set([("D", "E"),
                                                                         ("F", "G")]))

        relevant_probabilistic_activities = determine_relevant_probabilistic_activities(pstn, ("C", "K"))
        self.assertSetEqual(set(relevant_probabilistic_activities), set([("B", "C"),
                                                                         ("D", "E"),
                                                                         ("F", "G")]))

        relevant_probabilistic_activities = determine_relevant_probabilistic_activities(pstn, ("A", "L"))
        self.assertSetEqual(set(relevant_probabilistic_activities), set([("B", "C"),
                                                                         ("D", "E"),
                                                                         ("F", "G"),
                                                                         ("X", "A")]))

    def test03_relevant_probabilistic_activities(self):

        INF = 100000

        pstn = PSTN.from_requirements_and_activities([("A", (0, 15), "L"),
                                                      ("C", (0, 2), "K"),
                                                      ("H", (0, 2), "I"),
                                                      ("G", (0, 3), "I")],
                                                     [("C", (0, INF), "K"),
                                                      ("K", (0.5, 1), "L"),
                                                      ("E", (0.5, 1), "H"),
                                                      ("H", (0, INF), "I"),
                                                      ("I", (0.5, 1), "J"),
                                                      ("J", (0, INF), "K"),
                                                      ("G", (0, INF), "I"),
                                                      ("A", (0, INF), "B"),
                                                      ("Y", (0, INF), "D"),
                                                      ("Y", (0, INF), "F"),
                                                      ("B", ProbDistr.uniform((10, 14)), "C"),
                                                      ("D", ProbDistr.uniform((8, 11)), "E"),
                                                      ("F", ProbDistr.uniform((2, 4)), "G"),
                                                      ("X", ProbDistr.uniform((2, 4)), "A"),
                                                      ("A", ProbDistr.uniform((0.5, 1)), "Y")])
  
        relevant_probabilistic_activities = determine_relevant_probabilistic_activities(pstn, ("H", "I"))
        self.assertSetEqual(set(relevant_probabilistic_activities), set([("D", "E"),
                                                                         ("F", "G")]))

        relevant_probabilistic_activities = determine_relevant_probabilistic_activities(pstn, ("G", "I"))
        self.assertSetEqual(set(relevant_probabilistic_activities), set([("D", "E"),
                                                                         ("F", "G")]))

        relevant_probabilistic_activities = determine_relevant_probabilistic_activities(pstn, ("C", "K"))
        self.assertSetEqual(set(relevant_probabilistic_activities), set([("B", "C"),
                                                                         ("D", "E"),
                                                                         ("F", "G"),
                                                                         ("A", "Y")]))

        relevant_probabilistic_activities = determine_relevant_probabilistic_activities(pstn, ("A", "L"))
        self.assertSetEqual(set(relevant_probabilistic_activities), set([("B", "C"),
                                                                         ("D", "E"),
                                                                         ("F", "G"),
                                                                         ("X", "A"), 
                                                                         ("A", "Y")]))

    def test04_relevant_probabilistic_activities(self):

        INF = 100000

        pstn = PSTN.from_requirements_and_activities([("A", (0, 15), "L"),
                                                      ("C", (0, 2), "K"),
                                                      ("I", (0, 2), "H"),
                                                      ("G", (0, 3), "I")],
                                                     [("C", (0, INF), "K"),
                                                      ("K", (0.5, 1), "L"),
                                                      ("E", (0.5, 1), "H"),
                                                      ("I", (0, INF), "H"),
                                                      ("H", (0.5, 1), "J"),
                                                      ("J", (0, INF), "K"),
                                                      ("G", (0, INF), "I"),
                                                      ("A", (0, INF), "B"),
                                                      ("Y", (0, INF), "D"),
                                                      ("Y", (0, INF), "F"),
                                                      ("B", ProbDistr.uniform((10, 14)), "C"),
                                                      ("D", ProbDistr.uniform((8, 11)), "E"),
                                                      ("F", ProbDistr.uniform((2, 4)), "G"),
                                                      ("X", ProbDistr.uniform((2, 4)), "A"),
                                                      ("A", ProbDistr.uniform((0.5, 1)), "Y")])
  
        relevant_probabilistic_activities = determine_relevant_probabilistic_activities(pstn, ("I", "H"))
        self.assertSetEqual(set(relevant_probabilistic_activities), set([("D", "E"),
                                                                         ("F", "G")]))

        relevant_probabilistic_activities = determine_relevant_probabilistic_activities(pstn, ("G", "I"))
        self.assertSetEqual(set(relevant_probabilistic_activities), set([("F", "G"),
                                                                         ("A", "Y"),
                                                                         ("X", "A")]))

        relevant_probabilistic_activities = determine_relevant_probabilistic_activities(pstn, ("C", "K"))
        self.assertSetEqual(set(relevant_probabilistic_activities), set([("B", "C"),
                                                                         ("D", "E"),
                                                                         ("F", "G"),
                                                                         ("A", "Y")]))

        relevant_probabilistic_activities = determine_relevant_probabilistic_activities(pstn, ("A", "L"))
        self.assertSetEqual(set(relevant_probabilistic_activities), set([("B", "C"),
                                                                         ("D", "E"),
                                                                         ("F", "G"),
                                                                         ("X", "A"),
                                                                         ("A", "Y")]))

    def test05_reformulated_chance_constraints(self):

        INF = 100000

        pstn = PSTN.from_requirements_and_activities([("A", (0, 15), "L"),
                                                      ("C", (0, 2), "K"),
                                                      ("H", (0, 2), "I"),
                                                      ("G", (0, 3), "I")],
                                                     [("C", (0, INF), "K"),
                                                      ("K", (0.5, 1), "L"),
                                                      ("E", (0.5, 1), "H"),
                                                      ("H", (0, INF), "I"),
                                                      ("I", (0.5, 1), "J"),
                                                      ("J", (0, INF), "K"),
                                                      ("G", (0, INF), "I"),
                                                      ("A", (0, INF), "B"),
                                                      ("A", (0, INF), "D"),
                                                      ("A", (0, INF), "F"),
                                                      ("B", ProbDistr.uniform((10, 14)), "C"),
                                                      ("D", ProbDistr.uniform((8, 11)), "E"),
                                                      ("F", ProbDistr.uniform((2, 4)), "G")])
  
        reformulated_chance_constraints = reformulate_chance_constraints(pstn,
                                                                         [ChanceConstraint((("A", "L"),), 0.02)])

        self.assertTrue(len(reformulated_chance_constraints) == 1)
        self.assertSetEqual(set(reformulated_chance_constraints[0].ra_vars),
                            set([(('F', 'G'), False),
                                 (('F', 'G'), True),
                                 (('B', 'C'), False),
                                 (('B', 'C'), True),
                                 (('D', 'E'), False),
                                 (('D', 'E'), True)]))
        self.assertTrue(reformulated_chance_constraints[0].risk == 0.02)

        reformulated_chance_constraints = reformulate_chance_constraints(pstn,
                                                                         [ChanceConstraint((("A", "L"),), 0.05),
                                                                          ChanceConstraint((("H", "I"),), 0.02)])

        self.assertTrue(len(reformulated_chance_constraints) == 2)
        self.assertSetEqual(set(reformulated_chance_constraints[0].ra_vars),
                            set([(('F', 'G'), False),
                                 (('F', 'G'), True),
                                 (('B', 'C'), False),
                                 (('B', 'C'), True),
                                 (('D', 'E'), False),
                                 (('D', 'E'), True)]))
        self.assertTrue(reformulated_chance_constraints[0].risk == 0.05)
        self.assertSetEqual(set(reformulated_chance_constraints[1].ra_vars),
                            set([(('F', 'G'), False),
                                 (('F', 'G'), True),
                                 (('D', 'E'), False),
                                 (('D', 'E'), True)]))
        self.assertTrue(reformulated_chance_constraints[1].risk == 0.02)

        reformulated_chance_constraints = reformulate_chance_constraints(pstn,
                                                                         [ChanceConstraint((("G", "I"),), 0.02),
                                                                          ChanceConstraint((("H", "I"),), 0.05)])
        self.assertTrue(len(reformulated_chance_constraints) == 1)
        self.assertSetEqual(set(reformulated_chance_constraints[0].ra_vars),
                            set([(('F', 'G'), False),
                                 (('F', 'G'), True),
                                 (('D', 'E'), False),
                                 (('D', 'E'), True)]))
        self.assertTrue(reformulated_chance_constraints[0].risk == 0.02)

        reformulated_chance_constraints = reformulate_chance_constraints(pstn,
                                                                         [ChanceConstraint((("H", "I"),), 0.05),
                                                                          ChanceConstraint((("G", "I"),), 0.02)])
        self.assertTrue(len(reformulated_chance_constraints) == 1)
        self.assertSetEqual(set(reformulated_chance_constraints[0].ra_vars),
                            set([(('F', 'G'), False),
                                 (('F', 'G'), True),
                                 (('D', 'E'), False),
                                 (('D', 'E'), True)]))
        self.assertTrue(reformulated_chance_constraints[0].risk == 0.02)

    def test06_build_and_solve_ra_mip(self):

        INF = 100000

        pstn = PSTN.from_requirements_and_activities([("A", (0, 15), "L"),
                                                      ("C", (0, 2), "K"),
                                                      ("H", (0, 2), "I"),
                                                      ("G", (0, 3), "I")],
                                                     [("C", (0, INF), "K"),
                                                      ("K", (0.5, 1), "L"),
                                                      ("E", (0.5, 1), "H"),
                                                      ("H", (0, INF), "I"),
                                                      ("I", (0.5, 1), "J"),
                                                      ("J", (0, INF), "K"),
                                                      ("G", (0, INF), "I"),
                                                      ("A", (0, INF), "B"),
                                                      ("A", (0, INF), "D"),
                                                      ("A", (0, INF), "F"),
                                                      ("B", ProbDistr.uniform((10, 14)), "C"),
                                                      ("D", ProbDistr.uniform((8, 11)), "E"),
                                                      ("F", ProbDistr.uniform((2, 4)), "G")])
  
        reformulated_chance_constraints = reformulate_chance_constraints(pstn,
                                                                         [ChanceConstraint((("A", "L"),), 0.02)])

        (feasible,
         ra_vars_assignments,
         ra_conflict_resolution_constraint) = build_and_solve_ra_mip(pstn,
                                                                     reformulated_chance_constraints,
                                                                     [],
                                                                     [])
    def test07_check_cc_pstn_risk_bounded_dc(self):

        INF = 100000

        pstn = PSTN.from_requirements_and_activities([("A", (0, 15), "L"),
                                                      ("C", (0, 2), "K"),
                                                      ("H", (0, 2), "I"),
                                                      ("G", (0, 3), "I")],
                                                     [("C", (0, INF), "K"),
                                                      ("K", (0.5, 1), "L"),
                                                      ("E", (0.5, 1), "H"),
                                                      ("H", (0, INF), "I"),
                                                      ("I", (0.5, 1), "J"),
                                                      ("J", (0, INF), "K"),
                                                      ("G", (0, INF), "I"),
                                                      ("A", (0, INF), "B"),
                                                      ("A", (0, INF), "D"),
                                                      ("A", (0, INF), "F"),
                                                      ("B", ProbDistr.uniform((10, 14)), "C"),
                                                      ("D", ProbDistr.uniform((6, 11)), "E"),
                                                      ("F", ProbDistr.uniform((2, 4)), "G")])
  
        reformulated_chance_constraints = reformulate_chance_constraints(pstn,
                                                                         [ChanceConstraint((("A", "L"),), 0.1)])

        (ra_found,
         ra_vars_assignments,
         projection_stnu,
         constrs) = check_pstn_cc_dc(pstn,
                                                  reformulated_chance_constraints)

        self.assertFalse(ra_found)
        #print(ra_found)
        #print(ra_vars_assignments)
        #print(projection_stnu)
        #print(constrs)

        reformulated_chance_constraints = reformulate_chance_constraints(pstn,
                                                                         [ChanceConstraint((("A", "L"),), 0.6)])

        (ra_found,
         ra_vars_assignments,
         projection_stnu,
         constrs) = check_pstn_cc_dc(pstn,
                                                  reformulated_chance_constraints)

        self.assertTrue(ra_found)
        #print(ra_found)
        #print(ra_vars_assignments)
        #print(projection_stnu)
        #print(constrs)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
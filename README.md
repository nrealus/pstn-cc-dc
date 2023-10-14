# PSTN risk-bounded / chance-constrained dynamic controllability

An implementation of the algorithm for risk-bounded (aka chance-constrained) dynamic controllability checking of Probabilistic Simple Temporal Networks (PSTN) from [1].

In addition to risk-bounded dynamic controllability checking, risk-bounded dynamic execution policies for PSTNs, as defined in [1], could be easily extracted from the
outputs of the main function of the algorithm as well. However, this is descoped and not implemented in this repository.

On the other hand, we have implemented some improvements and extensions discussed in chapter 8 of [1]. The first is that we use a MIP solver
with piecewise linearization, instead of a NLP solver. The second is the support for multiple separate chance constraints, no
just a single global one. The fourth one is that we use a SAT solver to perform unit propagation when enumerating / branching over
linear constraints composing disjunctive linear constraints. Finally, when the selected (linear) constraints lead
to infeasibility of the MIP, [1] advises to use an irreducible infeasible set (IIS) of constraints, but descopes it,
since IIS computation for LP/MIP and especially NLP is often only available in the most high-end commercial solvers. We use a simple general
additive-deletion filter [2] to return an IIS of linear constraints selected for the MIP. Indeed, there is no need to include the reformulated chance constraint and linearization constraints in our IIS, because they must always hold at the top level.

Also, please note that our implementation differs a bit from the one presented in [1].
Indeed, instead of explicitly implementing the 3 layers described in [1], we "flattened" them into one function.

Final note: the code is still a bit unclean, and almost uncommented / undocumented. This should be improved. Also, it would be nice to have some more tests.

[1]: [Wang, A.J. Risk-Bounded Dynamic Scheduling of Temporal Plans, 2022 (PhD Thesis)](https://dspace.mit.edu/handle/1721.1/147542)
[2]: [Chinneck, J.W. Feasibility and Infeasibility in Optimization, 2007 (Tutorial for CP-AI-OR-07)](https://www.sce.carleton.ca/faculty/chinneck/docs/CPAIOR07InfeasibilityTutorial.pdf)

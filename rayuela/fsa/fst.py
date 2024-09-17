from frozendict import frozendict
from itertools import product


from rayuela.base.semiring import Boolean, Semiring
from rayuela.base.misc import epsilon_filter
from rayuela.base.symbol import Sym, ε, ε_1, ε_2

from rayuela.fsa.fsa import FSA
from rayuela.fsa.state import State, PairState


class FST(FSA):

    def __init__(self, R=Boolean):

        # DEFINITION
        # A weighted finite-state transducer is a 8-tuple <Σ, Δ, Q, F, I, δ, λ, ρ> where
        # • Σ is an alphabet of symbols;
        # • Δ is an alphabet of symbols;
        # • Q is a finite set of states;
        # • I ⊆ Q is a set of initial states;
        # • F ⊆ Q is a set of final states;
        # • δ is a finite relation Q × Σ × Δ × Q × R;
        # • λ is an initial weight function;
        # • ρ is a final weight function.

        # NOTATION CONVENTIONS
        # • single states (elements of Q) are denoted q
        # • multiple states not in sequence are denoted, p, q, r, ...
        # • multiple states in sequence are denoted i, j, k, ...
        # • symbols (elements of Σ and Δ) are denoted lowercase a, b, c, ...
        # • single weights (elements of R) are denoted w
        # • multiple weights (elements of R) are denoted u, v, w, ...

        super().__init__(R=R)

        # alphabet of output symbols
        self.Delta = set()

    def add_arc(self, i, a, b, j, w=None):
        if w is None: w = self.R.one

        if not isinstance(i, State): i = State(i)
        if not isinstance(j, State): j = State(j)
        if not isinstance(a, Sym): a = Sym(a)
        if not isinstance(b, Sym): b = Sym(b)
        if not isinstance(w, self.R): w = self.R(w)

        self.add_states([i, j])
        self.Sigma.add(a)
        self.Delta.add(b)
        self.δ[i][(a, b)][j] += w

    def set_arc(self, i, a, b, j, w=None):
        if w is None: w = self.R.one

        if not isinstance(i, State): i = State(i)
        if not isinstance(j, State): j = State(j)
        if not isinstance(a, Sym): a = Sym(a)
        if not isinstance(b, Sym): b = Sym(b)
        if not isinstance(w, self.R): w = self.R(w)

        self.add_states([i, j])
        self.Sigma.add(a)
        self.Delta.add(b)
        self.δ[i][(a, b)][j] = w

    def freeze(self):
        self.Sigma = frozenset(self.Sigma)
        self.Delta = frozenset(self.Delta)
        self.Q = frozenset(self.Q)
        self.δ = frozendict(self.δ)
        self.λ = frozendict(self.λ)
        self.ρ = frozendict(self.ρ)

    def arcs(self, i, no_eps=False):
        for ab, T in self.δ[i].items():
            if no_eps and ab == (ε, ε):
                continue
            for j, w in T.items():
                if w == self.R.zero:
                    continue
                yield ab, j, w
    
    def reverse(self):
        """creates a reversed machine"""

        # create the new machine
        Tr = self.spawn()

        # add the arcs in the reversed machine
        for i in self.Q:
            for (a, b), j, w in self.arcs(i):
                Tr.add_arc(j, a, b, i, w)

        # reverse the initial and final states
        for q, w in self.I:
            Tr.set_F(q, w)
        for q, w in self.F:
            Tr.set_I(q, w)

        return Tr

    def spawn(self, keep_init: bool = False, keep_final: bool = False):
        """returns a new FST in the same semiring"""
        F = FST(R=self.R)

        if keep_init:
            for q, w in self.I:
                F.set_I(q, w)
        if keep_final:
            for q, w in self.F:
                F.set_F(q, w)

        return F
    
    def trim(self):
        """trims the machine"""

        # compute accessible and co-accessible arcs
        A, C = self.accessible(), self.coaccessible()
        AC = A.intersection(C)

        # create a new F with only the pruned arcs
        T = self.spawn()
        for i in AC:
            for (a, b), j, w in self.arcs(i):
                if j in AC:
                    T.add_arc(i, a, b, j, w)

        # add initial state
        for q, w in self.I:
            if q in AC:
                T.set_I(q, w)

        # add final state
        for q, w in self.F:
            if q in AC:
                T.set_F(q, w)

        return T
    
    def _transform_arc(self, q: State, a: Sym, b: Sym, j: State, w, idx: int):
        if idx == 1:
            if b != ε:
                return a, b
            else:
                return a, ε_2
        else:
            if a != ε:
                return a, b
            else:
                return ε_1, b

    def augment_epsilon_transitions(self, idx: int):
        """Augments the FST by changing the appropriate epsilon transitions to
        epsilon_1 or epsilon_2 transitions to be able to perform the composition
        with epsilon transitions correctly.
        See also Fig. 7 in Mohri, Weighted Automata Algorithms, p. 17.

        Args:
            idx (int): 1 if the FST is the first one in the composition, 2 otherwise.

        Returns:
            FST: The augmented FST.
        """
        assert idx in [1, 2]

        T = self.spawn()

        for q in self.Q:
            if idx == 1:
                T.add_arc(q, ε, ε_1, q, w=self.R.one)
            else:
                T.add_arc(q, ε_2, ε, q, w=self.R.one)

            for (a, b), j, w in self.arcs(q):
                _a, _b = self._transform_arc(q, a, b, j, w, idx)
                T.add_arc(q, _a, _b, j, w=w)

        for q, w in self.I:
            T.set_I(q, w=w)

        for q, w in self.F:
            T.set_F(q, w=w)

        return T


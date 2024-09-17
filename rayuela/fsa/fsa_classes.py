"""Contains functions generating WFSAs of various types."""
from rayuela.base.semiring import Semiring
from rayuela.base.symbol import Sym, ε_1, ε_2
from rayuela.fsa.fsa import FSA
from rayuela.fsa.fst import FST
from rayuela.fsa.state import State

def get_epsilon_filter(R, Sigma):
        """Returns the epsilon filtered FST required for the correct composition of WFSTs
        with epsilon transitions.
        
        Returns:
        FST: The 3-state epsilon filtered WFST.
        """
        
        F = FST(R)
        
        # 0 ->
        for a in Sigma:
            F.add_arc(State(0), a, a, State(0), R.one)
        F.add_arc(State(0), ε_2, ε_1, State(0), R.one)
        F.add_arc(State(0), ε_1, ε_1, State(1), R.one)
        F.add_arc(State(0), ε_2, ε_2, State(2), R.one)
        
        # 1 ->
        for a in Sigma:
            F.add_arc(State(1), a, a, State(0), R.one)
        F.add_arc(State(1), ε_1, ε_1, State(1), R.one)
        
        # 2 ->
        for a in Sigma:
            F.add_arc(State(2), a, a, State(0), R.one)
        F.add_arc(State(2), ε_2, ε_2, State(2), R.one)
        
        F.set_I(State(0), R.one)
        F.set_F(State(0), R.one)
        F.set_F(State(1), R.one)
        F.set_F(State(2), R.one)
        
        return F
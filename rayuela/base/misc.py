from math import sqrt
import random
import string
import numpy as np
from fractions import Fraction


def _random_weight(semiring, **kwargs):
    from rayuela.base.semiring import (
        Real,
        Tropical,
        String
    )

    if semiring is String:
        str_len = int(random.random() * 8 + 1)
        return semiring(
            "".join(random.choice(string.ascii_lowercase) for _ in range(str_len))
        )

    elif semiring is Real:
        tol = 1e-3
        s = kwargs.get("divide_by", 6)
        random_weight = round(random.random() / s, 3)
        while random_weight < sqrt(tol):
            random_weight = round(random.random() / s, 3)
        return semiring(random_weight)

    elif semiring is Tropical:
        return semiring(random.randint(0, 50))


def random_weight_negative(semiring):
    from rayuela.base.semiring import Tropical

    if semiring is Tropical:
        return semiring(random.randint(-50, 50))
    else:
        raise AssertionError("Unsupported Semiring")


def epsilon_filter(a1, a2, q3):
    """
    Filter for composition with epsilon transitions
    """
    from rayuela.fsa.state import State
    from rayuela.base.symbol import ε_1, ε_2

    if a1 == a2 and a1 not in [ε_1, ε_2]:
        return State("0")
    elif a1 == ε_2 and a2 == ε_1 and q3 == State("0"):
        return State("0")
    elif a1 == ε_1 and a2 == ε_1 and q3 != State("2"):
        return State("1")
    elif a1 == ε_2 and a2 == ε_2 and q3 != State("1"):
        return State("2")
    else:
        return State("⊥")


def is_pathsum_positive(fsa):
    from rayuela.fsa.fsa import FSA
    from rayuela.fsa.fst import FST

    assert isinstance(fsa, FSA) or isinstance(fsa, FST)

    return fsa.pathsum() > fsa.R.zero


def filter_negative_pathsums(list_of_fsas):
    return [fsa for fsa in list_of_fsas if is_pathsum_positive(fsa)]


def compare_fsas(original_fsa, student_fsa) -> bool:
    from rayuela.fsa.fsa import FSA
    from rayuela.fsa.fst import FST

    assert isinstance(original_fsa, FSA) or isinstance(original_fsa, FST)
    assert isinstance(student_fsa, FSA) or isinstance(student_fsa, FST)

    if is_pathsum_positive(original_fsa):
        # TODO: Change check for: there is no an arbitrary number of initial states
        same_number_initial_states = len(list(original_fsa.I)) == len(
            list(student_fsa.I)
        )
        return np.allclose(
            float(original_fsa.pathsum()), float(student_fsa.pathsum()), atol=1e-3
        )  # and same_number_initial_states --> This would break some correct implementations
    # Skip non-convergent pathsums
    return True


def ansi(color=None, light=None, bg=3):
    return "\x1b[%s;%s%sm" % (light, bg, color)


_reset = "\x1b[0m"


def colorstring(s, c):
    return c + s + _reset


class colors:
    black, red, green, yellow, blue, magenta, cyan, white = [
        colorstring("%s", ansi(c, 0)) for c in range(8)
    ]

    class light:
        black, red, green, yellow, blue, magenta, cyan, white = [
            colorstring("%s", ansi(c, 1)) for c in range(8)
        ]

    class dark:
        black, red, green, yellow, blue, magenta, cyan, white = [
            colorstring("%s", ansi(c, 2)) for c in range(8)
        ]

    class bg:
        black, red, green, yellow, blue, magenta, cyan, white = [
            colorstring("%s", ansi(c, 0, bg=4)) for c in range(8)
        ]

    normal = "\x1b[0m%s\x1b[0m"
    bold = "\x1b[1m%s\x1b[0m"
    italic = "\x1b[3m%s\x1b[0m"
    underline = "\x1b[4m%s\x1b[0m"
    strike = "\x1b[9m%s\x1b[0m"
    # overline = lambda x: (u''.join(unicode(c) + u'\u0305' for c in unicode(x))).encode('utf-8')

    leftarrow = "←"
    rightarrow = "→"
    reset = _reset

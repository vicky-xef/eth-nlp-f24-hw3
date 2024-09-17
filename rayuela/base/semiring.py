from typing import Optional, Type

import numpy as np

from fractions import Fraction
from collections import defaultdict as dd
from frozendict import frozendict


# base code from https://github.com/timvieira/hypergraphs/blob/master/hypergraphs/semirings/boolean.py
class Semiring:

    zero: "Semiring"
    one: "Semiring"
    idempotent = False

    def __init__(self, score):
        self.score = score

    @classmethod
    def zeros(cls, mat_shape):
        import numpy as np

        return np.full(mat_shape, cls.zero)

    @classmethod
    def chart(cls, default=None):
        if default is None:
            default = cls.zero
        return dd(lambda: default)

    @classmethod
    def diag(cls, N):
        W = cls.zeros((N, N))
        for n in range(N):
            W[n, n] = cls.one

        return W

    def __add__(self, other):
        raise NotImplementedError

    def __mul__(self, other):
        raise NotImplementedError

    def __eq__(self, other):
        return self.score == other.score

    def __hash__(self):
        return hash(self.score)


class Boolean(Semiring):
    def __init__(self, score):
        super().__init__(score)

    def star(self):
        return Boolean.one

    def __add__(self, other):
        return Boolean(self.score or other.score)

    def __mul__(self, other):
        if other.score is self.one:
            return self.score
        if self.score is self.one:
            return other.score
        if other.score is self.zero:
            return self.zero
        if self.score is self.zero:
            return self.zero
        return Boolean(other.score and self.score)

    # TODO: is this correct?
    def __invert__(self):
        return Boolean.one

    def __truediv__(self, other):
        return Boolean.one

    def __eq__(self, other):
        return self.score == other.score

    def __lt__(self, other):
        return self.score < other.score

    def __repr__(self):
        return f"{self.score}"

    def __str__(self):
        return str(self.score)

    def __hash__(self):
        return hash(self.score)


Boolean.zero = Boolean(False)
Boolean.one = Boolean(True)
Boolean.idempotent = True
# TODO: check
Boolean.cancellative = True


class String(Semiring):
    def __init__(self, score):
        super().__init__(score)

    def star(self):
        return String.one

    def __add__(self, other):
        from rayuela.base.misc import lcp

        if other is self.zero:
            return self
        if self is self.zero:
            return other
        return String(lcp(self.score, other.score))

    def __mul__(self, other):
        if other is self.one:
            return self
        if self is self.one:
            return other
        if other is self.zero:
            return self.zero
        if self is self.zero:
            return self.zero
        return String(self.score + other.score)

    def __truediv__(self, other):
        from rayuela.base.misc import lcp

        prefix = lcp(self.score, other.score)
        return String(self.score[len(prefix) :])

    def __eq__(self, other):
        return self.score == other.score

    def __repr__(self):
        return f"{self.score}"

    def __hash__(self):
        return hash(self.score)

# unique "infinity" string
String.zero = String("∞")
# empty string
String.one = String("")
String.idempotent = False
String.cancellative = False


class Tropical(Semiring):
    def __init__(self, score):
        self.score = score

    def star(self):
        return self.one

    def __float__(self):
        return float(self.score)

    def __int__(self):
        return int(self.score)

    def __add__(self, other):
        return Tropical(min(self.score, other.score))

    def __mul__(self, other):
        if other is self.one:
            return self
        if self is self.one:
            return other
        if other is self.zero:
            return self.zero
        if self is self.zero:
            return self.zero
        return Tropical(self.score + other.score)

    def __invert__(self):
        return Tropical(-self.score)

    def __truediv__(self, other):
        return Tropical(self.score - other.score)

    def __lt__(self, other):
        return self.score < other.score

    def __repr__(self):
        return f"Tropical({self.score})"

    def __str__(self):
        return str(self.score)


Tropical.zero = Tropical(float("inf"))
Tropical.one = Tropical(0.0)
Tropical.idempotent = True
Tropical.superior = True
Tropical.cancellative = True


class Real(Semiring):
    def __init__(self, score):
        # TODO: this is hack to deal with the fact
        # that we have to hash weights
        self.score = score

    def star(self):
        return Real(1.0 / (1.0 - self.score))

    def __float__(self):
        return float(self.score)

    def __add__(self, other):
        return Real(self.score + other.score)

    def __sub__(self, other):
        return Real(self.score - other.score)

    def __mul__(self, other):
        if other is self.one:
            return self
        if self is self.one:
            return other
        if other is self.zero:
            return self.zero
        if self is self.zero:
            return self.zero
        return Real(self.score * other.score)

    def __invert__(self):
        return Real(1.0 / self.score)

    def __truediv__(self, other):
        return Real(self.score / other.score)

    def __lt__(self, other):
        return self.score < other.score

    def __repr__(self):
        # return f'Real({self.score})'
        return f"{round(self.score, 15)}"

    def __eq__(self, other):
        # return float(self.score) == float(other.score)
        return np.allclose(float(self.score), float(other.score), atol=1e-3)

    # TODO: find out why this wasn't inherited
    def __hash__(self):
        return hash(self.score)


Real.zero = Real(0.0)
Real.one = Real(1.0)
Real.idempotent = False
Real.cancellative = True

class ProductSemiring(Semiring):
    def __init__(self, x, y):
        super().__init__((x, y))

    def star(self):
        raise NotImplemented

    def __add__(self, other):
        w1, w2 = self.score[0], other.score[0]
        v1, v2 = self.score[1], other.score[1]
        return ProductSemiring(w1 + w2, v1 + v2)

    def __mul__(self, other):
        w1, w2 = self.score[0], other.score[0]
        v1, v2 = self.score[1], other.score[1]
        return ProductSemiring(w1 * w2, v1 * v2)

    def __truediv__(self, other):
        w1, w2 = self.score[0], other.score[0]
        v1, v2 = self.score[1], other.score[1]
        return ProductSemiring(w1 / w2, v1 / v2)

    def __invert__(self):
        return ProductSemiring(~self.score[0], ~self.score[1])

    def __eq__(self, other):
        return self.score == other.score

    def __repr__(self):
        if isinstance(self.score[0], String):
            # the imporant special case of encoding transducers
            if len(self.score[0].score) > 0:
                return f"{self.score[0]} / {self.score[1]}"
            else:
                return f"{self.score[1]}"
        return f"〈{self.score[0]}, {self.score[1]}〉"

    def __hash__(self):
        return hash(self.score)

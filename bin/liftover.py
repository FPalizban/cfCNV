"""Minimal UCSC chain-file liftOver for segments (no external dependency).
Segment ends are mapped through the highest-scoring chain block containing them (or the nearest
block inside the segment); segments whose ends land on different chromosomes or whose length changes
by more than 2-fold are dropped, as with liftOver -minMatch-style filtering."""
import bisect, gzip
from collections import defaultdict

class Lift:
    def __init__(self, chain):
        self.b = defaultdict(list)
        with gzip.open(chain, "rt") as f:
            hdr = None
            for l in f:
                x = l.split()
                if not x: continue
                if x[0] == "chain":
                    sc, tN, tS, tSt, tB = float(x[1]), x[2], int(x[3]), x[4], int(x[5])
                    qN, qS, qSt, qB = x[7], int(x[8]), x[9], int(x[10]); t, q = tB, qB; hdr = (sc, tN, qN, qS, qSt); continue
                sc, tN, qN, qS, qSt = hdr; n = int(x[0])
                self.b[tN].append((t, t + n, q, qN, qS, qSt, sc))
                if len(x) == 3: t += n + int(x[1]); q += n + int(x[2])
        for k in list(self.b):
            self.b[k].sort(); self.b[k + "_s"] = [y[0] for y in self.b[k]]

    def _blocks(self, c, a, z):  # blocks overlapping [a, z)
        B, S = self.b.get(c, []), self.b.get(c + "_s", [])
        i = max(0, bisect.bisect_right(S, a) - 50); out = []
        while i < len(B) and B[i][0] < z:
            if B[i][1] > a: out.append(B[i])
            i += 1
        return out

    def _pos(self, blk, p):
        t0, t1, q0, qN, qS, qSt, _ = blk; q = q0 + (p - t0)
        return qN, (q if qSt == "+" else qS - q - 1)

    def seg(self, c, s, e):
        bl = self._blocks(c, s, e)
        if not bl: return None
        top = max(b[6] for b in bl); bl = [b for b in bl if b[6] == top]   # best chain only
        first = min(bl, key=lambda b: b[0]); last = max(bl, key=lambda b: b[1])
        n1, a = self._pos(first, max(s, first[0])); n2, z = self._pos(last, min(e, last[1]) - 1)
        if n1 != n2: return None
        a, z = min(a, z), max(a, z) + 1
        if not (0.5 <= (z - a) / max(1, e - s) <= 2): return None
        return n1, a, z

def _pieces(self, c, s, e):
    """aligned pieces of [s, e) in the target build: one per chain block; where chains overlap, the
    higher-scoring chain wins (covered regions tracked in a sorted list -> O(n log n))"""
    out, cs, ce = [], [], []
    for b in sorted(self._blocks(c, s, e), key=lambda b: (-b[6], b[0])):
        a, z = max(s, b[0]), min(e, b[1])
        if z <= a: continue
        i = bisect.bisect_left(ce, a + 1)
        if i < len(cs) and cs[i] < z: continue
        cs.insert(i, a); ce.insert(i, z)
        n1, qa = self._pos(b, a); _, qz = self._pos(b, z - 1)
        out.append((n1, min(qa, qz), max(qa, qz) + 1))
    return out
Lift.pieces = _pieces

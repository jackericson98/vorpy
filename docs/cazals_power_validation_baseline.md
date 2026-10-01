# Cazals / Power validation baseline

This is the frozen validation baseline for VorPy's weighted Power/Laguerre
interface analyzer. It records tested observations and the unresolved whole-
interface discrepancy; it does not claim exact reproduction of every result
in Cazals et al. (2006).

## Condition β interpretation

For two probe-expanded atom balls (B_i=(p_i,R_i)) and
(B_j=(p_j,R_j)), equality of their power distances at a candidate center
(x) gives

\[
q(x)=\|x-p_i\|^2-R_i^2=\|x-p_j\|^2-R_j^2=m^2.
\]

If (q>0), a real orthogonal ball exists with radius (m=\sqrt q). If
(q=0), the only solution has radius zero. If (q<0), no real ball can
satisfy the orthogonality equations at that center. For a bounded convex
Power facet, (q) is convex on its plane, so its maximum over the polygon is
attained at a polygon vertex. The 2KAI audit evaluates the maximum at all
incident power vertices.

The paper states the rejection condition (m/r>M), where (r) is the smaller
expanded atom-ball radius and (M=5), but does not spell out the case where
all candidate (q\le0). The mathematically literal interpretation is that no
real positive-radius orthogonal ball exists there, so the facet cannot meet a
rejection condition requiring such a ball to exceed the threshold. The
production implementation represents this as `m=0` for the comparison; this
is a decision-equivalent reporting placeholder, not a claim that a real
radius-zero orthogonal ball exists when (q<0). We retain this behavior and
do not change the selected interface.

For the 362 final 2KAI facets: 88 have a positive maximum (q), none have
exactly zero maximum, and 274 have a negative maximum with no real candidate.
Among the 88 positive candidates, the minimum (m/r) is 0.0825127; the maximum
over all facets is 3.413599. Rejected facets at (M=5): 0.

## Synthetic checks

- A known 30-degree singleton angle and 5 Å edge gives exact
  \(\beta=\pi/6\) and \(l\beta=5\pi/6\).
- Swapping groups reverses signed curvature and preserves unsigned curvature.
- Generator ordering does not change the result.
- Condition-β tests cover positive, zero, and negative (m^2), values on both
  sides of (m/r=5), translation, rotation, generator swapping, and uniform
  scaling.

## 1UDI local-angle baseline

- Direct singleton-angle beta audit: 945/945.
- Covalent same-group pairs: N=294; mean 17.386°; median 17.414°;
  284/294 lie in 12–24°.
- Noncovalent same-group pairs: N=651; mean 47.216°; median 45.993°;
  638/651 lie in 20–80°.

## 2KAI whole-interface baseline

Groups are protease/kallikrein chains A+B versus inhibitor/BPTI chain I.
The Cazals/Chothia radii and 1.40 Å probe are used, with alpha=0 and M=5.

- Protein atoms retained: 2,236.
- Interface atoms: protease 131; inhibitor 65; total 196.
- \(r_{AB}=2.01538\).
- Final facets: 362; VIA=873.762685 Å².
- Connected components=1; significant components=1.
- Interior curvature edges=864; boundary interface edges=111.
- Polygonal/combinatorial interior-edge agreement: 864/864.
- Beta singleton-angle agreement: 864/864.
- \(L=953.524161\) Å.
- \(C_{signed}=-127.050015\) Å·rad; \(C_{unsigned}=607.364973\) Å·rad.
- With A=protease and B=inhibitor, signed \(s_H=-7.634237°\).
- Reversing A/B gives \(+7.634237°\), changing sign only.
- Published Cazals observation: approximately \(+17°\).

VorPy has **not reproduced the published +17-degree whole-interface result**.
The local beta, selected-facet incidence, condition-β outcome, and geometric
edge lengths have been audited independently. The remaining magnitude
discrepancy is unresolved. The Power analyzer is validated for its measured
local geometry and current interface construction, not as an exact reproduction
of all Cazals results.

The raw Cazals turning-angle sum is \(C_{raw}=\sum_e l_e\beta_e\). The
conventional edge contribution for \(H=(k_1+k_2)/2\) is
\(C_H=\frac12\sum_e l_e\beta_e\). These are both retained and labeled
separately.

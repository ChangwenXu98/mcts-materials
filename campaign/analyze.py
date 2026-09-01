"""
Summarise a running (or finished) campaign into JSON for the slide deck.

Reads the QE evaluator's live CSV cache plus the search's progress snapshot,
ranks candidates by the Tc fit, and reports the structural characteristics that
distinguish the top of the list from the rest.

    python analyze.py <results_dir> [top_n]
"""

import json
import sys
from pathlib import Path

import pandas as pd

from mcts_framework.superhydride import belli2025_tc
from mcts_framework.superhydride.elements import classify_host
from mcts_framework.superhydride.rewards import PHI_STAR_OPTIMUM

#: The compression the campaign is defined at, and the window the pressure
#: match accepts. Candidates outside it were measured elsewhere on the
#: isotherm and are flagged rather than silently ranked alongside the rest.
TARGET_PRESSURE_GPA = 150.0
PRESSURE_TOLERANCE_GPA = 10.0

#: Pauling electronegativities for the host elements this campaign can reach.
#: Used only to describe the winners, never to score them.
ELECTRONEGATIVITY = {
    "Li": 0.98, "Na": 0.93, "K": 0.82, "Rb": 0.82, "Cs": 0.79,
    "Be": 1.57, "Mg": 1.31, "Ca": 1.00, "Sr": 0.95, "Ba": 0.89,
    "Sc": 1.36, "Y": 1.22, "La": 1.10, "Lu": 1.27,
    "Ti": 1.54, "Zr": 1.33, "Hf": 1.30, "V": 1.63, "Nb": 1.60, "Ta": 1.50,
    "Al": 1.61, "B": 2.04, "C": 2.55, "N": 3.04,
    "Si": 1.90, "P": 2.19, "S": 2.58,
    "Ga": 1.81, "Ge": 2.01, "As": 2.18, "Se": 2.55,
    "In": 1.78, "Sn": 1.96, "Sb": 2.05, "Te": 2.10,
    "Tl": 1.62, "Pb": 2.33, "Bi": 2.02,
}


def split_hosts(formula: str):
    """Return the non-hydrogen element symbols in a formula like 'BeLaH8'."""
    import re

    return [e for e, _ in re.findall(r"([A-Z][a-z]?)(\d*)", formula) if e and e != "H"]


def load(results_dir: Path) -> pd.DataFrame:
    cache = results_dir / "qe_descriptors.csv"
    if not cache.exists():
        return pd.DataFrame()
    frame = pd.read_csv(cache)
    ok = frame[frame["status"] == "ok"].copy()
    if ok.empty:
        return ok
    ok["tc_k"] = [
        belli2025_tc(r.phi, r.phi_star, r.h_f, r.h_dos) for r in ok.itertuples()
    ]
    ok["hosts"] = ok["formula"].map(split_hosts)
    ok["classes"] = ok["hosts"].map(
        lambda hs: sorted({classify_host(_z(h)) or "other" for h in hs})
    )
    ok["chi_min"] = ok["hosts"].map(
        lambda hs: min((ELECTRONEGATIVITY.get(h, float("nan")) for h in hs), default=float("nan"))
    )
    ok["chi_max"] = ok["hosts"].map(
        lambda hs: max((ELECTRONEGATIVITY.get(h, float("nan")) for h in hs), default=float("nan"))
    )
    ok["phi_gap"] = ok["phi_star"] - ok["phi"]
    # The pressure match has a finite SCF budget; a candidate that ran out of
    # it was measured at a different compression from the rest and must not be
    # compared to them without saying so.
    ok["pressure_off_gpa"] = (ok["pressure_gpa"] - TARGET_PRESSURE_GPA).abs()
    ok["pressure_ok"] = ok["pressure_off_gpa"] <= PRESSURE_TOLERANCE_GPA
    ok["phi_star_offset"] = (ok["phi_star"] - PHI_STAR_OPTIMUM).abs()
    return ok.sort_values("tc_k", ascending=False).reset_index(drop=True)


def _z(symbol: str) -> int:
    from ase.data import atomic_numbers

    return atomic_numbers.get(symbol, 0)


#: Below this many completed candidates, a top-vs-bottom comparison is noise.
MIN_FOR_CHARACTERISATION = 9


def characterise(ranked: pd.DataFrame, top_n: int) -> dict:
    """
    Compare the top of the ranking against the bottom, descriptor by descriptor.

    The two groups are the top and bottom third, so the contrast is between
    comparable sample sizes. Ranking eight candidates against one - which a
    naive head/tail split gives early in a campaign - says nothing.
    """
    if len(ranked) < MIN_FOR_CHARACTERISATION:
        return {}
    group = max(3, min(top_n, len(ranked) // 3))
    top = ranked.head(group)
    rest = ranked.tail(group)

    def gap(column):
        return {
            "top_mean": round(float(top[column].mean()), 4),
            "rest_mean": round(float(rest[column].mean()), 4),
        }

    classes = {}
    for row in top.itertuples():
        classes["+".join(row.classes)] = classes.get("+".join(row.classes), 0) + 1

    return {
        "phi": gap("phi"),
        "phi_star": gap("phi_star"),
        "phi_gap": gap("phi_gap"),
        "phi_star_offset": gap("phi_star_offset"),
        "h_dos": gap("h_dos"),
        "chi_min": gap("chi_min"),
        "class_counts_top": classes,
        "n_top": int(len(top)),
        "n_rest": int(len(rest)),
    }


def main() -> int:
    results_dir = Path(sys.argv[1])
    top_n = int(sys.argv[2]) if len(sys.argv) > 2 else 10

    ranked = load(results_dir)
    progress = {}
    progress_file = results_dir / "progress.json"
    if progress_file.exists():
        progress = json.loads(progress_file.read_text())

    cache = results_dir / "qe_descriptors.csv"
    attempted = failed = 0
    if cache.exists():
        allrows = pd.read_csv(cache)
        attempted = int(len(allrows))
        failed = int((allrows["status"] != "ok").sum())

    payload = {
        "progress": progress,
        "attempted": attempted,
        "succeeded": int(len(ranked)),
        "failed": failed,
        "off_pressure": int((~ranked["pressure_ok"]).sum()) if len(ranked) else 0,
        "top": [
            {
                "formula": r.formula,
                "phi": round(r.phi, 4),
                "phi_star": round(r.phi_star, 4),
                "h_f": round(r.h_f, 3),
                "h_dos": round(r.h_dos, 4),
                "pressure_gpa": round(r.pressure_gpa, 1),
                "tc_k": round(r.tc_k, 1),
                "classes": r.classes,
                "pressure_ok": bool(r.pressure_ok),
            }
            for r in ranked.head(top_n).itertuples()
        ],
        "characteristics": characterise(ranked, top_n),
    }
    print(json.dumps(payload, indent=2))
    (results_dir / "analysis.json").write_text(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())

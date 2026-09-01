"""
Render the live campaign results into the slide deck.

Rewrites the marked regions of the deck HTML from analysis.json, so the deck
can be republished at any point while the job is still running.

    python update_deck.py <analysis.json> <deck.html>
"""

import html
import json
import sys
from pathlib import Path

from mcts_framework.superhydride.rewards import PHI_STAR_OPTIMUM, TC_MAX_K


def sub_formula(formula: str) -> str:
    """LaBeH8 -> LaBeH<sub>8</sub>."""
    out, digits = [], ""
    for ch in formula:
        if ch.isdigit():
            digits += ch
        else:
            if digits:
                out.append(f"<sub>{digits}</sub>")
                digits = ""
            out.append(ch)
    if digits:
        out.append(f"<sub>{digits}</sub>")
    return "".join(out)


def tree_payload(results_dir: Path) -> dict:
    """
    Flatten tree.json into what the deck's viewer needs.

    Reward is carried through as-is (0.0 marks a candidate whose funnel failed)
    and the descriptors come along so a node can be inspected without a second
    lookup.
    """
    path = results_dir / "tree.json"
    if not path.exists():
        return {"nodes": [], "root": None}
    raw = json.loads(path.read_text())
    nodes = []
    for n in raw.get("nodes", []):
        props = n.get("properties") or {}
        nodes.append({
            "id": n["id"],
            "parent": n["parent"],
            "name": str(n.get("identifier", "?")).split("|")[0],
            "reward": n.get("own_reward"),
            "visits": n.get("visits", 0),
            "terminated": bool(n.get("terminated")),
            "phi": props.get("phi"),
            "phi_star": props.get("phi_star"),
            "h_dos": props.get("h_dos"),
            "p": props.get("pressure_gpa"),
        })
    return {"nodes": nodes, "root": raw.get("root_id")}


ALKALI = {"Li", "Na", "K", "Rb", "Cs"}


def points_payload(results_dir: Path) -> dict:
    """Per-candidate scatter data: the descriptors, the fit, and one category."""
    import csv
    import math
    import re

    from mcts_framework.superhydride.rewards import belli2025_tc

    path = results_dir / "qe_descriptors.csv"
    if not path.exists():
        return {"points": []}
    points = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            if not str(row.get("status", "")).startswith("ok"):
                continue
            try:
                phi, phi_star = float(row["phi"]), float(row["phi_star"])
                h_f, h_dos = float(row["h_f"]), float(row["h_dos"])
                pressure = float(row["pressure_gpa"])
            except (TypeError, ValueError):
                continue
            hosts = {e for e, _ in re.findall(r"([A-Z][a-z]?)(\d*)", row["formula"])
                     if e and e != "H"}
            points.append({
                "f": row["formula"],
                "phi": round(phi, 4),
                "ps": round(phi_star, 4),
                "hdos": round(h_dos, 4),
                "p": round(pressure, 1),
                "tc": round(belli2025_tc(phi, phi_star, h_f, h_dos), 1),
                "alkali": bool(hosts & ALKALI),
                "pok": abs(pressure - 150.0) <= 10.0,
            })
    points.sort(key=lambda d: -d["tc"])
    return {"points": points}


def replace_region(text: str, name: str, body: str) -> str:
    start, end = f"<!--{name}_START-->", f"<!--{name}_END-->"
    head, _, tail = text.partition(start)
    _, _, tail = tail.partition(end)
    return f"{head}{start}\n{body}\n      {end}{tail}"


def results_region(data: dict) -> str:
    progress = data.get("progress", {})
    top = data.get("top", [])
    finished = progress.get("finished", False)
    succeeded, failed = data.get("succeeded", 0), data.get("failed", 0)
    off = data.get("off_pressure", 0)

    heading = (
        "Campaign complete" if finished
        else ("Campaign in flight" if not top else "Ranking so far")
    )
    best = top[0] if top else None
    best_value = f"{best['tc_k']:.0f} K" if best else "&mdash;"
    best_sub = (
        f"{sub_formula(best['formula'])} &middot; &plusmn;41 K fit RMSE"
        if best else "&plusmn;41 K fit RMSE"
    )

    if top:
        rows = []
        for i, r in enumerate(top):
            cls = ' class="best"' if i == 0 else (' class="flag"' if not r["pressure_ok"] else "")
            flag = "" if r["pressure_ok"] else ' <span class="pill warn">off P</span>'
            rows.append(
                f'<tr{cls}><td class="name">{sub_formula(r["formula"])}{flag}</td>'
                f'<td class="num">{r["phi"]:.3f}</td>'
                f'<td class="num">{r["phi_star"]:.3f}</td>'
                f'<td class="num">{r["h_f"]:.2f}</td>'
                f'<td class="num">{r["h_dos"]:.3f}</td>'
                f'<td class="num">{r["pressure_gpa"]:.0f}</td>'
                f'<td class="num"><strong>{r["tc_k"]:.0f}</strong></td></tr>'
            )
        body_rows = "\n            ".join(rows)
    else:
        body_rows = ('<tr><td colspan="7" style="text-align:center; color:var(--muted); '
                     'padding:26px">Awaiting first completions&hellip;</td></tr>')

    note = (
        f"{succeeded} candidate{'s' if succeeded != 1 else ''} through the funnel, "
        f"{failed} failed (mostly missing pseudopotentials on the lanthanide chain). "
    )
    if off:
        note += (
            f"<strong>{off}</strong> landed outside &plusmn;10&nbsp;GPa of the 150&nbsp;GPa "
            f"target and are flagged &mdash; those were measured at a different compression. "
        )
    note += ("The QE evaluator appends every completed candidate to its CSV cache, so this "
             "is a live read of the campaign rather than a post-hoc summary.")

    return f"""      <h2 id="results-heading">{heading}</h2>
      <div class="cols c3" style="margin:6px 0 18px">
        <div class="card"><span class="k">Candidates evaluated</span><div class="v cool">{succeeded}</div><div class="sub">{failed} failed &middot; {progress.get('iterations', 0)} search iterations</div></div>
        <div class="card"><span class="k">Best estimate</span><div class="v warm">{best_value}</div><div class="sub">{best_sub}</div></div>
        <div class="card"><span class="k">Tree size</span><div class="v">{progress.get('unique_materials', 0)}</div><div class="sub">unique compositions attached</div></div>
      </div>
      <div class="tablewrap">
        <table>
          <thead><tr><th>compound</th><th class="num">&phi;</th><th class="num">&phi;*</th><th class="num">H<sub>f</sub></th><th class="num">H<sub>DOS</sub></th><th class="num">P<sub>DFT</sub> (GPa)</th><th class="num">T<sub>c</sub> (K)</th></tr></thead>
          <tbody>
            {body_rows}
          </tbody>
        </table>
      </div>
      <p class="note" style="margin-top:16px">{note}</p>"""


def struct_region(data: dict) -> str:
    ch = data.get("characteristics") or {}
    top = data.get("top", [])
    if not ch or len(top) < 4:
        return """      <h2 id="struct-heading">Structural read &mdash; pending completions</h2>
      <p>This slide fills in once enough candidates have cleared the funnel to say something that is not noise. The three things worth reading off the ranking:</p>
      <ul>
        <li><strong>&phi;* against the 2/3 optimum.</strong> Candidates whose molecularity index drifts toward 1 carry intact H&sub2; and are penalised by construction, not by evidence.</li>
        <li><strong>The &phi;*&nbsp;&minus;&nbsp;&phi; gap.</strong> The two diverge exactly when the hydrogen interactions are inhomogeneous &mdash; a uniform network opens at one threshold.</li>
        <li><strong>Host electronegativity.</strong> Whether the winners cluster in the electropositive or the covalent class, and what that implies about the pressure they would need.</li>
      </ul>"""

    def card(label, key, note, fmt="{:.3f}"):
        top_v, rest_v = ch[key]["top_mean"], ch[key]["rest_mean"]
        arrow = "higher" if top_v > rest_v else "lower"
        colour = "cool" if top_v > rest_v else "warm"
        return f"""        <div class="card"><span class="k">{label}</span>
          <div class="v {colour}">{fmt.format(top_v)}</div>
          <div class="sub">vs {fmt.format(rest_v)} in the rest &mdash; {arrow} in the top {ch['n_top']}</div>
          <p style="margin-top:10px">{note}</p></div>"""

    classes = ch.get("class_counts_top", {})
    class_line = ", ".join(
        f"<strong>{html.escape(k)}</strong> &times;{v}" for k, v in sorted(classes.items(), key=lambda kv: -kv[1])
    ) or "&mdash;"

    return f"""      <h2 id="struct-heading">What separates the top {ch['n_top']} from the other {ch['n_rest']}</h2>
      <div class="cols c2" style="margin-top:8px">
{card("&phi;* &mdash; molecularity", "phi_star", "The fit peaks at &phi;* = 2/3 and vanishes at 1. A top candidate sitting near the optimum has hydrogen stretched but still bonded; one near 1 carries intact H&sub2;.")}
{card("|&phi;* &minus; 2/3| &mdash; distance from optimum", "phi_star_offset", "The single most direct read of whether the winners win on molecularity or in spite of it.")}
{card("&phi; &mdash; networking", "phi", "How high an ELF threshold still spans the crystal. Enters the fit only as a cube root, so it discriminates weakly compared with &phi;*.")}
{card("H<sub>DOS</sub> &mdash; H share at E<sub>F</sub>", "h_dos", "Also a cube-root term, but the one that separates a metallic hydrogen network from a host-dominated Fermi surface.")}
      </div>
      <p class="note" style="margin-top:18px">Host classes among the top {ch['n_top']}: {class_line}. Mean lowest host electronegativity {ch['chi_min']['top_mean']:.2f} against {ch['chi_min']['rest_mean']:.2f} in the rest &mdash; the campaign's own read on whether electropositive or covalent hosts are winning. &phi;*&nbsp;&minus;&nbsp;&phi; gap: {ch['phi_gap']['top_mean']:.3f} vs {ch['phi_gap']['rest_mean']:.3f}.</p>"""


def status_region(data: dict) -> str:
    progress = data.get("progress", {})
    finished = progress.get("finished", False)
    succeeded = data.get("succeeded", 0)
    best = data.get("top", [{}])[0].get("tc_k") if data.get("top") else None
    label = (
        f"complete &middot; {succeeded} evaluated" if finished
        else f"running &middot; {succeeded} evaluated"
    )
    if best:
        label += f" &middot; best {best:.0f} K"
    dot = "dot" if finished else "dot run"
    return (f'      <div class="status"><span class="{dot}"></span>'
            f'<span id="status-text">{label}</span></div>')


def main() -> int:
    data = json.loads(Path(sys.argv[1]).read_text())
    deck_path = Path(sys.argv[2])
    text = deck_path.read_text()
    text = replace_region(text, "RESULTS", results_region(data))
    text = replace_region(text, "STRUCT", struct_region(data))
    text = replace_region(text, "STATUS", status_region(data))
    results_dir = Path(sys.argv[1]).parent
    text = replace_region(
        text, "POINTS",
        '      <script id="points-data" type="application/json">'
        + json.dumps(points_payload(results_dir), separators=(",", ":")) + "</script>",
    )
    tree = tree_payload(results_dir)
    text = replace_region(
        text, "TREE",
        '      <script id="tree-data" type="application/json">'
        + json.dumps(tree, separators=(",", ":")) + "</script>",
    )
    deck_path.write_text(text)
    print(f"deck updated: {len(tree['nodes'])} tree nodes, "
          f"{data.get('succeeded', 0)} evaluated, "
          f"{len(data.get('top', []))} in the table "
          f"(fit max {TC_MAX_K:.1f} K, phi* optimum {PHI_STAR_OPTIMUM:.3f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())

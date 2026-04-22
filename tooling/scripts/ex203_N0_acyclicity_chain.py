#!/usr/bin/env python3
"""
ex203_N0_acyclicity_chain.py
============================
Verifies that the logical chain of N0 derivations (sek10_N0_wyprowadzenie)
has NO circular dependencies.

Builds a directed acyclic graph (DAG) of all derivation nodes, checks for
cycles via topological sort, validates parentage, and detects floating nodes.

Author : TGP automated verification
Date   : 2026-04-12
"""

import sys
import io
from collections import defaultdict, deque

# --- UTF-8 stdout -----------------------------------------------------------
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

# =============================================================================
# 1.  Define the derivation graph
# =============================================================================

# Each entry: (source, target)
# "source -> target" means "source is used to derive target"

EDGES = [
    # --- Core N0 chain (sek10) -----------------------------------------------
    ("A1_space_from_matter",          "AXIOM_FOUNDATION"),
    ("A2_substrate_Gamma_Z2",         "AXIOM_FOUNDATION"),

    ("A2_substrate_Gamma_Z2",         "H_Gamma"),
    ("H_Gamma",                       "Lem_K_phi2"),
    ("Lem_K_phi2",                    "K0_eq_0"),

    ("A2_substrate_Gamma_Z2",         "GL_u4_neg"),
    ("GL_u4_neg",                     "Thm_psi6_necessity"),
    ("Thm_psi6_necessity",            "Cor_psi3_tgp"),

    ("Lem_K_phi2",                    "Thm_beta_gamma"),

    ("A2_substrate_Gamma_Z2",         "Thm_N0_unique"),

    ("K0_eq_0",                       "Thm_N0_instability"),
    ("Cor_psi3_tgp",                  "Thm_N0_instability"),

    # --- Phi0 determinations ------------------------------------------------
    ("Lambda_obs",                    "Phi0_from_Lambda_NUM"),
    ("Phi0_Nf2_algebraic",           "Phi0_at_Nc3_AN"),
    ("ERG_self_consistent",           "Phi0_from_ERG_AN_HYP"),

    # --- CG bridge chain ----------------------------------------------------
    ("A2_substrate_Gamma_Z2",         "CG1_block_averaging"),
    ("CG1_block_averaging",           "LemA1_compactness_RK"),
    ("CG1_block_averaging",           "LemA2_local_functional"),
    ("W4_correlation_decay",          "LemA2_local_functional"),
    ("Lem_K_phi2",                    "LemA3_alpha_eff_2"),
    ("A2_substrate_Gamma_Z2",         "LemA4_beta_gamma_id"),
    ("LemA3_alpha_eff_2",             "LemA4_beta_gamma_id"),
    ("LemA1_compactness_RK",          "LemA5_weak_continuum"),
    ("LemA2_local_functional",        "LemA5_weak_continuum"),
    ("LemA3_alpha_eff_2",             "LemA5_weak_continuum"),
    ("LemA4_beta_gamma_id",           "LemA5_weak_continuum"),

    # --- Constants chain ----------------------------------------------------
    ("Phi0_from_Lambda_NUM",          "kappa_from_Phi0"),
    ("G0_from_soliton",              "kappa_from_Phi0"),

    ("Phi0_from_Lambda_NUM",          "G0_from_soliton"),
    ("soliton_profile",              "G0_from_soliton"),

    ("Phi0_from_Lambda_NUM",          "Lambda_eff_dark_energy"),
    ("Thm_beta_gamma",               "Lambda_eff_dark_energy"),

    ("Phi0_from_Lambda_NUM",          "lepton_masses"),
    ("soliton_spectrum",              "lepton_masses"),
]

# Axiom / external-input nodes (have no parents by design)
AXIOMS = {
    "A1_space_from_matter",
    "A2_substrate_Gamma_Z2",
    "AXIOM_FOUNDATION",       # sink collecting axiom tags
    "Lambda_obs",             # observational input [NUM]
    "Phi0_Nf2_algebraic",    # algebraic identity [AN]
    "ERG_self_consistent",    # ERG closure [AN+HYP]
    "W4_correlation_decay",   # assumption on correlation decay
    "soliton_profile",        # soliton solution input
    "soliton_spectrum",       # soliton spectrum input
}

# =============================================================================
# 2.  Build adjacency structures
# =============================================================================

def build_graph(edges):
    """Return children dict, parents dict, and full node set."""
    children = defaultdict(set)
    parents  = defaultdict(set)
    nodes    = set()
    for src, tgt in edges:
        children[src].add(tgt)
        parents[tgt].add(src)
        nodes.add(src)
        nodes.add(tgt)
    return children, parents, nodes

# =============================================================================
# 3.  Topological sort  (Kahn's algorithm) — cycle detection
# =============================================================================

def topological_sort(nodes, parents, children):
    """
    Returns (sorted_list, has_cycle).
    If has_cycle is True, sorted_list contains only the nodes that could be
    sorted; remaining nodes are involved in a cycle.
    """
    in_degree = {n: len(parents.get(n, set())) for n in nodes}
    queue = deque(n for n in nodes if in_degree[n] == 0)
    order = []
    while queue:
        node = queue.popleft()
        order.append(node)
        for child in children.get(node, set()):
            in_degree[child] -= 1
            if in_degree[child] == 0:
                queue.append(child)
    has_cycle = len(order) != len(nodes)
    return order, has_cycle

# =============================================================================
# 4.  Ancestor check — no derived result feeds back into its own ancestor
# =============================================================================

def ancestors_of(node, parents):
    """BFS upward to collect all ancestors of a node."""
    visited = set()
    queue = deque(parents.get(node, set()))
    while queue:
        n = queue.popleft()
        if n in visited:
            continue
        visited.add(n)
        for p in parents.get(n, set()):
            queue.append(p)
    return visited

def check_no_ancestor_circularity(nodes, parents, children):
    """For every edge src->tgt, verify tgt is NOT an ancestor of src."""
    violations = []
    for src in children:
        src_ancestors = ancestors_of(src, parents)
        for tgt in children[src]:
            if tgt in src_ancestors:
                violations.append((src, tgt))
    return violations

# =============================================================================
# 5.  Floating-node detection
# =============================================================================

def find_reachable_from(start_nodes, children):
    """BFS forward from start_nodes, return all reachable."""
    visited = set()
    queue = deque(start_nodes)
    while queue:
        n = queue.popleft()
        if n in visited:
            continue
        visited.add(n)
        for c in children.get(n, set()):
            queue.append(c)
    return visited

# =============================================================================
# 6.  Run all checks
# =============================================================================

def main():
    children, parents, nodes = build_graph(EDGES)

    # --- basic stats ---------------------------------------------------------
    n_nodes = len(nodes)
    n_edges = len(EDGES)
    axiom_nodes = AXIOMS & nodes
    n_axioms = len(axiom_nodes)

    print("=" * 72)
    print("  ex203  N\u2080 ACYCLICITY CHAIN VERIFICATION")
    print("=" * 72)
    print(f"\n  Total nodes  : {n_nodes}")
    print(f"  Total edges  : {n_edges}")
    print(f"  Axiom nodes  : {n_axioms}")
    print()

    results = {}

    # --- CHECK 1: cycle detection via topological sort -----------------------
    order, has_cycle = topological_sort(nodes, parents, children)
    if has_cycle:
        cycle_nodes = nodes - set(order)
        print("[FAIL] CHECK 1 — Cycle detected!")
        print(f"       Nodes involved in cycle: {sorted(cycle_nodes)}")
        results["cycle_free"] = False
    else:
        print("[PASS] CHECK 1 — No cycles detected (topological sort succeeded).")
        results["cycle_free"] = True

    # --- CHECK 2: every non-axiom node has >= 1 parent -----------------------
    orphans = []
    for n in nodes:
        if n not in AXIOMS and len(parents.get(n, set())) == 0:
            orphans.append(n)
    if orphans:
        print(f"[FAIL] CHECK 2 — Orphan (parentless non-axiom) nodes: {sorted(orphans)}")
        results["all_have_parents"] = False
    else:
        print("[PASS] CHECK 2 — Every non-axiom node has at least one parent.")
        results["all_have_parents"] = True

    # --- CHECK 3: no ancestor circularity ------------------------------------
    violations = check_no_ancestor_circularity(nodes, parents, children)
    if violations:
        print(f"[FAIL] CHECK 3 — Ancestor circularity found:")
        for s, t in violations:
            print(f"       {s} -> {t}  (but {t} is ancestor of {s})")
        results["no_ancestor_circ"] = False
    else:
        print("[PASS] CHECK 3 — No derived result feeds back into its own ancestor.")
        results["no_ancestor_circ"] = True

    # --- CHECK 4: floating nodes (disconnected from axiom base) --------------
    reachable = find_reachable_from(axiom_nodes, children)
    floating = nodes - reachable
    if floating:
        print(f"[WARN] CHECK 4 — Floating nodes (not reachable from axioms):")
        for f in sorted(floating):
            print(f"       {f}")
        results["no_floating"] = False
    else:
        print("[PASS] CHECK 4 — All nodes reachable from axiom base.")
        results["no_floating"] = True

    # --- Topological order ---------------------------------------------------
    print("\n" + "-" * 72)
    print("  TOPOLOGICAL ORDER (derivation sequence):")
    print("-" * 72)
    for i, node in enumerate(order, 1):
        tag = " [AXIOM]" if node in AXIOMS else ""
        print(f"  {i:3d}. {node}{tag}")

    # --- Summary -------------------------------------------------------------
    all_pass = all(results.values())
    print("\n" + "=" * 72)
    print("  SUMMARY")
    print("=" * 72)
    for name, ok in results.items():
        status = "PASS" if ok else "FAIL"
        print(f"  {name:25s} : {status}")
    print()
    if all_pass:
        print("  >>> ALL CHECKS PASSED — derivation chain is acyclic. <<<")
    else:
        print("  >>> SOME CHECKS FAILED — review issues above. <<<")
    print("=" * 72)

    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()

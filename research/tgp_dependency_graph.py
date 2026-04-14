#!/usr/bin/env python3
"""
tgp_dependency_graph.py — Generuje graf zależności TGP_v1
=========================================================

Tworzy:
1. Graf zależności .tex (main.tex → sekcje → dodatki)
2. Graf zależności Python (moduły nbody + scripts)
3. Graf konceptualny (twierdzenia → predykcje → eksperymenty)
4. Eksport: DOT, PNG (matplotlib), GEXF (Gephi-compatible)

Uruchomienie: python research/tgp_dependency_graph.py
"""

import os
import re
import sys
import json
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


# ======================================================================
# 1. GRAF .TEX — include tree z main.tex
# ======================================================================

def build_tex_graph():
    """Parse main.tex for \\input{} statements and build include tree."""
    G = nx.DiGraph()
    main_tex = os.path.join(REPO, "main.tex")

    if not os.path.exists(main_tex):
        print("  main.tex not found, skipping tex graph")
        return G

    with open(main_tex, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()

    # Find all \input{...} statements (not commented out)
    for line in content.split("\n"):
        line_stripped = line.strip()
        if line_stripped.startswith("%"):
            continue
        matches = re.findall(r'\\input\{([^}]+)\}', line)
        for m in matches:
            name = m.replace(".tex", "")
            G.add_edge("main", name)

            # Classify node
            if name.startswith("sek"):
                G.nodes[name]["type"] = "section"
                G.nodes[name]["color"] = "#4CAF50"
            elif name.startswith("dodatek"):
                G.nodes[name]["type"] = "appendix"
                G.nodes[name]["color"] = "#2196F3"
            elif name.startswith("nbody/"):
                G.nodes[name]["type"] = "nbody_tex"
                G.nodes[name]["color"] = "#FF9800"
            else:
                G.nodes[name]["type"] = "other"
                G.nodes[name]["color"] = "#9E9E9E"

    G.nodes["main"]["type"] = "root"
    G.nodes["main"]["color"] = "#F44336"

    # Now scan each section .tex for \input{} of sub-files
    for node in list(G.nodes()):
        if node == "main":
            continue
        tex_path = os.path.join(REPO, node + ".tex")
        if not os.path.exists(tex_path):
            continue
        with open(tex_path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.strip().startswith("%"):
                    continue
                for m in re.findall(r'\\input\{([^}]+)\}', line):
                    child = m.replace(".tex", "")
                    if child != node:
                        G.add_edge(node, child)
                        if "type" not in G.nodes.get(child, {}):
                            G.nodes[child]["type"] = "sub-include"
                            G.nodes[child]["color"] = "#FF9800"

    return G


# ======================================================================
# 2. GRAF PYTHON — import dependencies
# ======================================================================

def build_python_graph():
    """Parse nbody/*.py for inter-module imports."""
    G = nx.DiGraph()
    nbody_dir = os.path.join(REPO, "nbody")

    if not os.path.isdir(nbody_dir):
        print("  nbody/ not found, skipping Python graph")
        return G

    # Core modules
    py_files = [f for f in os.listdir(nbody_dir)
                if f.endswith(".py") and not f.startswith("__")]

    for py_file in py_files:
        module_name = py_file.replace(".py", "")
        G.add_node(module_name, type="module", color="#4CAF50")

        filepath = os.path.join(nbody_dir, py_file)
        with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                line = line.strip()
                if line.startswith("#"):
                    continue
                # from nbody.X import ...
                m = re.match(r'from\s+(?:nbody\.)?(\w+)\s+import', line)
                if m:
                    dep = m.group(1)
                    if dep != module_name and dep + ".py" in [x for x in py_files]:
                        G.add_edge(module_name, dep)
                # from .X import ...
                m = re.match(r'from\s+\.(\w+)\s+import', line)
                if m:
                    dep = m.group(1)
                    if dep != module_name:
                        G.add_edge(module_name, dep)

    # Classify by layer
    layer0 = {"tgp_field", "configurations", "dynamics_backends"}
    layer1 = {"pairwise", "bridge_nbody", "yukawa_from_defect",
              "three_body_terms", "three_body_force_exact",
              "multipole_triple_overlap", "soliton_interaction"}
    layer2 = {"nbody_energy", "equilibria", "stability", "lyapunov",
              "dynamics_v2", "eom_tgp"}

    for n in G.nodes():
        if n in layer0:
            G.nodes[n]["layer"] = 0
            G.nodes[n]["color"] = "#4CAF50"  # green
        elif n in layer1:
            G.nodes[n]["layer"] = 1
            G.nodes[n]["color"] = "#FF9800"  # orange
        elif n in layer2:
            G.nodes[n]["layer"] = 2
            G.nodes[n]["color"] = "#2196F3"  # blue
        else:
            G.nodes[n]["layer"] = 3
            G.nodes[n]["color"] = "#9C27B0"  # purple

    return G


# ======================================================================
# 3. GRAF KONCEPTUALNY — twierdzenia → predykcje → eksperymenty
# ======================================================================

def build_concept_graph():
    """Build conceptual graph: axioms → theorems → predictions → experiments."""
    G = nx.DiGraph()

    # AXIOMS (inputs)
    axioms = {
        "AK1: Z₂ substrat": {"type": "axiom", "color": "#F44336"},
        "AK2: H_Γ = -J·Σ(φᵢφⱼ)²": {"type": "axiom", "color": "#F44336"},
        "AK3: h(Φ)=Φ (p=1)": {"type": "axiom", "color": "#F44336", "status": "POSTULAT"},
        "PARAM: g₀ᵉ = 0.86941": {"type": "parameter", "color": "#E91E63"},
        "PARAM: Ω_Λ = 0.6847": {"type": "parameter", "color": "#E91E63"},
    }
    for name, attrs in axioms.items():
        G.add_node(name, **attrs)

    # THEOREMS (derived)
    theorems = {
        "TW: α=2 (kinetic)": {"type": "theorem", "color": "#4CAF50", "status": "SŁABE TW."},
        "TW: K(ρ)=ρ²": {"type": "theorem", "color": "#4CAF50", "status": "NUM"},
        "TW: β=γ (vacuum)": {"type": "theorem", "color": "#4CAF50", "status": "TWIERDZENIE"},
        "TW: d=3 (wymiar)": {"type": "theorem", "color": "#4CAF50", "status": "TWIERDZENIE"},
        "TW: m∝A⁴": {"type": "theorem", "color": "#CDDC39", "status": "POSTULAT+HEUR"},
        "TW: K(ℓ)=2/3 (Koide)": {"type": "theorem", "color": "#CDDC39", "status": "NUM"},
        "TW: B=√2 (Brannen)": {"type": "theorem", "color": "#CDDC39", "status": "NUM"},
        "TW: N=3 generacji": {"type": "theorem", "color": "#FF9800", "status": "HEURYSTYKA"},
        "TW: |GL(3,F₂)|=168": {"type": "theorem", "color": "#4CAF50", "status": "TWIERDZENIE"},
    }
    for name, attrs in theorems.items():
        G.add_node(name, **attrs)

    # PREDICTIONS
    predictions = {
        "F1: α_s = 0.1190": {"type": "prediction", "color": "#2196F3", "sigma": "1.1σ"},
        "F2: λ_C = 0.2282": {"type": "prediction", "color": "#FF5722", "sigma": "4.8σ"},
        "F3: K(ℓ)=2/3": {"type": "prediction", "color": "#2196F3", "sigma": "0σ"},
        "F6: Ω_DM = 0.262": {"type": "prediction", "color": "#2196F3", "sigma": "0.3σ"},
        "F7: Σm_ν = 59.6 meV": {"type": "prediction", "color": "#2196F3", "sigma": "testable"},
        "F9: n_s = 0.967": {"type": "prediction", "color": "#2196F3", "sigma": "0.4σ"},
        "F10: m_W = 80.354": {"type": "prediction", "color": "#2196F3", "sigma": "0.01σ"},
        "F11: m_H = 125.31": {"type": "prediction", "color": "#2196F3", "sigma": "0.3σ"},
        "P: r = 0.0033": {"type": "prediction", "color": "#2196F3", "sigma": "testable"},
        "P: w_DE ≥ -1": {"type": "prediction", "color": "#2196F3", "sigma": "testable"},
        "P: Normal ordering": {"type": "prediction", "color": "#2196F3", "sigma": "testable"},
        "P: GW breathing mode": {"type": "prediction", "color": "#2196F3", "sigma": "testable"},
    }
    for name, attrs in predictions.items():
        G.add_node(name, **attrs)

    # EXPERIMENTS (future tests)
    experiments = {
        "JUNO (~2028-30)": {"type": "experiment", "color": "#9C27B0"},
        "DESI DR3 (~2027)": {"type": "experiment", "color": "#9C27B0"},
        "LiteBIRD (~2028)": {"type": "experiment", "color": "#9C27B0"},
        "Hyper-K (~2026-28)": {"type": "experiment", "color": "#9C27B0"},
        "LISA/ET (~2035)": {"type": "experiment", "color": "#9C27B0"},
        "CMB-S4 (~2028)": {"type": "experiment", "color": "#9C27B0"},
    }
    for name, attrs in experiments.items():
        G.add_node(name, **attrs)

    # EDGES: axiom → theorem
    G.add_edge("AK1: Z₂ substrat", "TW: α=2 (kinetic)")
    G.add_edge("AK2: H_Γ = -J·Σ(φᵢφⱼ)²", "TW: α=2 (kinetic)")
    G.add_edge("AK2: H_Γ = -J·Σ(φᵢφⱼ)²", "TW: K(ρ)=ρ²")
    G.add_edge("AK1: Z₂ substrat", "TW: β=γ (vacuum)")
    G.add_edge("AK1: Z₂ substrat", "TW: d=3 (wymiar)")
    G.add_edge("AK3: h(Φ)=Φ (p=1)", "TW: m∝A⁴")

    # theorem → theorem
    G.add_edge("TW: α=2 (kinetic)", "TW: m∝A⁴")
    G.add_edge("TW: d=3 (wymiar)", "TW: m∝A⁴")
    G.add_edge("TW: m∝A⁴", "TW: K(ℓ)=2/3 (Koide)")
    G.add_edge("TW: m∝A⁴", "TW: B=√2 (Brannen)")
    G.add_edge("TW: B=√2 (Brannen)", "TW: K(ℓ)=2/3 (Koide)")
    G.add_edge("TW: d=3 (wymiar)", "TW: N=3 generacji")
    G.add_edge("TW: N=3 generacji", "TW: |GL(3,F₂)|=168")

    # theorem + param → prediction
    G.add_edge("PARAM: g₀ᵉ = 0.86941", "F1: α_s = 0.1190")
    G.add_edge("PARAM: Ω_Λ = 0.6847", "F1: α_s = 0.1190")
    G.add_edge("TW: |GL(3,F₂)|=168", "F1: α_s = 0.1190")

    G.add_edge("PARAM: Ω_Λ = 0.6847", "F2: λ_C = 0.2282")
    G.add_edge("TW: N=3 generacji", "F2: λ_C = 0.2282")

    G.add_edge("TW: K(ℓ)=2/3 (Koide)", "F3: K(ℓ)=2/3")

    G.add_edge("PARAM: Ω_Λ = 0.6847", "F6: Ω_DM = 0.262")
    G.add_edge("TW: N=3 generacji", "F6: Ω_DM = 0.262")

    G.add_edge("TW: K(ℓ)=2/3 (Koide)", "F7: Σm_ν = 59.6 meV")
    G.add_edge("TW: B=√2 (Brannen)", "F7: Σm_ν = 59.6 meV")

    G.add_edge("TW: α=2 (kinetic)", "F9: n_s = 0.967")
    G.add_edge("TW: α=2 (kinetic)", "P: r = 0.0033")

    G.add_edge("PARAM: g₀ᵉ = 0.86941", "F10: m_W = 80.354")
    G.add_edge("PARAM: g₀ᵉ = 0.86941", "F11: m_H = 125.31")

    G.add_edge("TW: α=2 (kinetic)", "P: w_DE ≥ -1")
    G.add_edge("TW: β=γ (vacuum)", "P: w_DE ≥ -1")

    G.add_edge("TW: K(ℓ)=2/3 (Koide)", "P: Normal ordering")
    G.add_edge("TW: α=2 (kinetic)", "P: GW breathing mode")

    # prediction → experiment
    G.add_edge("F7: Σm_ν = 59.6 meV", "JUNO (~2028-30)")
    G.add_edge("P: Normal ordering", "JUNO (~2028-30)")
    G.add_edge("P: w_DE ≥ -1", "DESI DR3 (~2027)")
    G.add_edge("P: r = 0.0033", "LiteBIRD (~2028)")
    G.add_edge("F9: n_s = 0.967", "CMB-S4 (~2028)")
    G.add_edge("P: GW breathing mode", "LISA/ET (~2035)")

    return G


# ======================================================================
# 4. RESEARCH PROBLEMS GRAPH
# ======================================================================

def build_research_graph():
    """Map research problems R1-R7 to their dependencies."""
    G = nx.DiGraph()

    problems = {
        "R1: Cabibbo\n(Ω_Λ/N)²": {"color": "#F44336", "priority": "NOW"},
        "R2: CG-1/3/4\ncontinuum": {"color": "#F44336", "priority": "LONG"},
        "R3: Why N=3": {"color": "#FF9800", "priority": "LONG"},
        "R4: h(Φ)=Φ\nmetric": {"color": "#F44336", "priority": "NOW"},
        "R5: m∝A⁴\nscaling": {"color": "#FF9800", "priority": "MED"},
        "R6: B=√2\nBrannen": {"color": "#FF9800", "priority": "MED"},
        "R7: UV\ncompletion": {"color": "#9E9E9E", "priority": "LOW"},
    }
    for name, attrs in problems.items():
        G.add_node(name, type="problem", **attrs)

    tools = {
        "galois\nGL(n,GF(2))": {"color": "#4CAF50", "type": "tool"},
        "Cadabra2\ntensors": {"color": "#4CAF50", "type": "tool"},
        "Lean 4\nproofs": {"color": "#4CAF50", "type": "tool"},
        "SymPy\nCAS": {"color": "#4CAF50", "type": "tool"},
        "SciPy\nODE": {"color": "#4CAF50", "type": "tool"},
    }
    for name, attrs in tools.items():
        G.add_node(name, **attrs)

    # problem → tool
    G.add_edge("R1: Cabibbo\n(Ω_Λ/N)²", "galois\nGL(n,GF(2))")
    G.add_edge("R3: Why N=3", "galois\nGL(n,GF(2))")
    G.add_edge("R4: h(Φ)=Φ\nmetric", "Cadabra2\ntensors")
    G.add_edge("R2: CG-1/3/4\ncontinuum", "Lean 4\nproofs")
    G.add_edge("R6: B=√2\nBrannen", "Lean 4\nproofs")
    G.add_edge("R5: m∝A⁴\nscaling", "SymPy\nCAS")
    G.add_edge("R5: m∝A⁴\nscaling", "SciPy\nODE")
    G.add_edge("R6: B=√2\nBrannen", "SciPy\nODE")
    G.add_edge("R7: UV\ncompletion", "SymPy\nCAS")

    # inter-problem dependencies
    G.add_edge("R2: CG-1/3/4\ncontinuum", "R5: m∝A⁴\nscaling",
               label="CG proves α=2 → k=4")
    G.add_edge("R5: m∝A⁴\nscaling", "R6: B=√2\nBrannen",
               label="k=4 → mass formula → B")
    G.add_edge("R6: B=√2\nBrannen", "R1: Cabibbo\n(Ω_Λ/N)²",
               label="Koide → CKM structure")

    return G


# ======================================================================
# VISUALIZATION
# ======================================================================

def draw_graph(G, title, filename, layout="spring", figsize=(16, 12)):
    """Draw a networkx graph and save as PNG."""
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    if layout == "spring":
        pos = nx.spring_layout(G, k=2.5, iterations=100, seed=42)
    elif layout == "kamada":
        pos = nx.kamada_kawai_layout(G)
    elif layout == "shell":
        pos = nx.shell_layout(G)
    elif layout == "hierarchical":
        pos = nx.multipartite_layout(G, subset_key="layer")
    else:
        pos = nx.spring_layout(G, seed=42)

    colors = [G.nodes[n].get("color", "#9E9E9E") for n in G.nodes()]
    sizes = []
    for n in G.nodes():
        t = G.nodes[n].get("type", "")
        if t in ("root", "axiom", "parameter"):
            sizes.append(1200)
        elif t in ("theorem", "problem"):
            sizes.append(900)
        elif t in ("prediction",):
            sizes.append(700)
        else:
            sizes.append(500)

    nx.draw_networkx_nodes(G, pos, node_color=colors, node_size=sizes,
                           alpha=0.9, ax=ax)
    nx.draw_networkx_edges(G, pos, edge_color="#BDBDBD",
                           arrows=True, arrowsize=15,
                           connectionstyle="arc3,rad=0.1", ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=6, font_weight="bold", ax=ax)

    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.axis("off")
    plt.tight_layout()
    plt.savefig(os.path.join(REPO, "research", filename), dpi=150,
                bbox_inches="tight")
    plt.close()
    print(f"  Saved: research/{filename}")


def export_dot(G, filename):
    """Export graph in DOT format."""
    path = os.path.join(REPO, "research", filename)
    nx.drawing.nx_pydot.write_dot(G, path)
    print(f"  Saved: research/{filename}")


def export_gexf(G, filename):
    """Export graph in GEXF format (Gephi-compatible)."""
    path = os.path.join(REPO, "research", filename)
    nx.write_gexf(G, path)
    print(f"  Saved: research/{filename}")


def print_stats(G, name):
    """Print basic graph statistics."""
    print(f"\n  {name}:")
    print(f"    Nodes: {G.number_of_nodes()}")
    print(f"    Edges: {G.number_of_edges()}")
    if nx.is_directed(G):
        if nx.is_directed_acyclic_graph(G):
            print(f"    DAG: YES")
            longest = nx.dag_longest_path(G)
            print(f"    Longest path: {len(longest)} nodes")
        else:
            print(f"    DAG: NO (contains cycles)")
            cycles = list(nx.simple_cycles(G))
            print(f"    Cycles: {len(cycles)}")


# ======================================================================
# MAIN
# ======================================================================

def main():
    print("=" * 60)
    print("TGP Dependency Graph Generator")
    print("=" * 60)

    # 1. TeX graph
    print("\n1. Building TeX include graph...")
    tex_g = build_tex_graph()
    print_stats(tex_g, "TeX Include Tree")
    draw_graph(tex_g, "TGP TeX Include Tree", "graph_tex_includes.png",
               layout="kamada", figsize=(20, 14))
    export_gexf(tex_g, "graph_tex_includes.gexf")

    # 2. Python graph
    print("\n2. Building Python module graph...")
    py_g = build_python_graph()
    print_stats(py_g, "Python Module Dependencies")
    draw_graph(py_g, "TGP Python Module Dependencies (nbody/)",
               "graph_python_modules.png", layout="kamada", figsize=(14, 10))
    export_gexf(py_g, "graph_python_modules.gexf")

    # 3. Concept graph
    print("\n3. Building conceptual graph...")
    concept_g = build_concept_graph()
    print_stats(concept_g, "Conceptual: Axioms -> Predictions -> Experiments")
    draw_graph(concept_g, "TGP: Axioms → Theorems → Predictions → Experiments",
               "graph_concept_flow.png", layout="kamada", figsize=(20, 14))
    export_gexf(concept_g, "graph_concept_flow.gexf")

    # 4. Research graph
    print("\n4. Building research problems graph...")
    research_g = build_research_graph()
    print_stats(research_g, "Research Problems R1-R7")
    draw_graph(research_g, "TGP Research Problems → Tools",
               "graph_research_problems.png", layout="spring", figsize=(14, 10))
    export_gexf(research_g, "graph_research_problems.gexf")

    # Summary
    print("\n" + "=" * 60)
    print("DONE — 4 graphs generated in research/")
    print("=" * 60)
    print("\nFiles:")
    print("  PNG: graph_tex_includes.png, graph_python_modules.png,")
    print("       graph_concept_flow.png, graph_research_problems.png")
    print("  GEXF (Gephi): *.gexf")


if __name__ == "__main__":
    main()

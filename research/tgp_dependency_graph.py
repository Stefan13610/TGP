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
    """Map ALL research folders (R1-R7, Q0-Q7, Galaxy, Cosmo) with dependencies."""
    G = nx.DiGraph()

    # --- Status colors ---
    # green (#4CAF50) = COMPLETE/SOLID, lime (#CDDC39) = BREAKTHROUGH,
    # orange (#FF9800) = PARTIAL/OPEN, red (#F44336) = RESOLVED/CLOSED,
    # grey (#9E9E9E) = EMPTY placeholder, blue (#2196F3) = OUTSIDE SCOPE

    # ===================== R-series: Core theory =====================
    r_nodes = {
        "R1: cabibbo\nRESOLVED 0.75σ": {
            "color": "#F44336", "status": "RESOLVED", "folder": "cabibbo_correction"},
        "R2: continuum\nOPEN 6-12mo": {
            "color": "#FF9800", "status": "OPEN", "folder": "continuum_limit"},
        "R3: why_n3\nSOLID 13/13": {
            "color": "#4CAF50", "status": "SOLID", "folder": "why_n3"},
        "R4: metric\nCOMPLETE 11/11": {
            "color": "#4CAF50", "status": "COMPLETE", "folder": "metric_ansatz"},
        "R5: mass_scaling\nBREAKTHROUGH 7/7": {
            "color": "#CDDC39", "status": "BREAKTHROUGH", "folder": "mass_scaling_k4"},
        "R6: brannen_sqrt2\nPARTIAL": {
            "color": "#FF9800", "status": "PARTIAL", "folder": "brannen_sqrt2"},
        "R7: uv_completion\nEMPTY": {
            "color": "#9E9E9E", "status": "EMPTY", "folder": "uv_completion"},
    }
    for name, attrs in r_nodes.items():
        G.add_node(name, type="problem_R", **attrs)

    # ===================== Q-series: Quantum mechanics =====================
    q_nodes = {
        "Q0: foundations\nCOMPLETE 12/12": {
            "color": "#4CAF50", "status": "COMPLETE", "folder": "qm_foundations"},
        "Q1: measurement\nCOMPLETE 22/22": {
            "color": "#4CAF50", "status": "COMPLETE", "folder": "qm_measurement"},
        "Q2: born_rule\nEMPTY placeholder": {
            "color": "#9E9E9E", "status": "EMPTY", "folder": "qm_born_rule"},
        "Q3: superposition\nCOMPLETE 7/7": {
            "color": "#4CAF50", "status": "COMPLETE", "folder": "qm_superposition"},
        "Q4: entanglement\nCOMPLETE 10/10": {
            "color": "#4CAF50", "status": "COMPLETE", "folder": "qm_entanglement"},
        "Q5: spin\nCOMPLETE 7/7": {
            "color": "#4CAF50", "status": "COMPLETE", "folder": "qm_spin"},
        "Q6: statistics\nCOMPLETE 8/8": {
            "color": "#4CAF50", "status": "COMPLETE", "folder": "qm_statistics"},
        "Q7: decoherence\nCOMPLETE 8/8": {
            "color": "#4CAF50", "status": "COMPLETE", "folder": "qm_decoherence"},
    }
    for name, attrs in q_nodes.items():
        G.add_node(name, type="problem_Q", **attrs)

    # ===================== Galaxy scaling =====================
    gal_node = "Galaxy: scaling\n46+ scripts PASS"
    G.add_node(gal_node, type="problem_G", color="#4CAF50",
               status="PASS", folder="galaxy_scaling")

    # ===================== Cosmological =====================
    cosmo_nodes = {
        "Cosmo: tensions\nCLOSED no mechanism": {
            "color": "#F44336", "status": "CLOSED", "folder": "cosmo_tensions"},
        "Cosmo: H0\nOUTSIDE SCOPE": {
            "color": "#2196F3", "status": "OUTSIDE_SCOPE", "folder": "hubble_tension"},
        "Cosmo: S8\nNEGLIGIBLE": {
            "color": "#2196F3", "status": "NEGLIGIBLE", "folder": "s8_tension"},
        "Cosmo: DESI w(z)\nINCOMPATIBLE": {
            "color": "#F44336", "status": "INCOMPATIBLE", "folder": "desi_dark_energy"},
    }
    for name, attrs in cosmo_nodes.items():
        G.add_node(name, type="problem_C", **attrs)

    # ===================== Tools =====================
    tools = {
        "galois\nGL(n,GF(2))": {"color": "#81D4FA", "type": "tool"},
        "Cadabra2\ntensors": {"color": "#81D4FA", "type": "tool"},
        "Lean 4\nproofs": {"color": "#81D4FA", "type": "tool"},
        "SymPy\nCAS": {"color": "#81D4FA", "type": "tool"},
        "SciPy\nODE": {"color": "#81D4FA", "type": "tool"},
    }
    for name, attrs in tools.items():
        G.add_node(name, **attrs)

    # ===================== R-series dependencies =====================
    # R2 (continuum) → ALL (foundational)
    for rn in r_nodes:
        if "R2:" not in rn:
            G.add_edge("R2: continuum\nOPEN 6-12mo", rn, label="foundational")
    # R5 (mass_scaling) → R3 (why_n3) via K^2 mechanism
    G.add_edge("R5: mass_scaling\nBREAKTHROUGH 7/7", "R3: why_n3\nSOLID 13/13",
               label="K^2 mechanism")
    # R6 (brannen) → R3 (why_n3) via Koide K=2/3
    G.add_edge("R6: brannen_sqrt2\nPARTIAL", "R3: why_n3\nSOLID 13/13",
               label="Koide K=2/3")
    # R4 (metric) → galaxy_scaling via metric ansatz
    G.add_edge("R4: metric\nCOMPLETE 11/11", gal_node,
               label="metric ansatz")

    # ===================== Q-series dependencies =====================
    # Q0 (foundations) → all Q modules (hbar derivation)
    for qn in q_nodes:
        if "Q0:" not in qn:
            G.add_edge("Q0: foundations\nCOMPLETE 12/12", qn, label="hbar")
    # Q5 (spin) → Q6 (statistics) via winding number B
    G.add_edge("Q5: spin\nCOMPLETE 7/7", "Q6: statistics\nCOMPLETE 8/8",
               label="winding number B")
    # Q1 (measurement) → Q4 (entanglement) via Born rule
    G.add_edge("Q1: measurement\nCOMPLETE 22/22", "Q4: entanglement\nCOMPLETE 10/10",
               label="Born rule")
    # Q3 (superposition) → Q7 (decoherence) via NL mixing
    G.add_edge("Q3: superposition\nCOMPLETE 7/7", "Q7: decoherence\nCOMPLETE 8/8",
               label="NL mixing")

    # ===================== Galaxy → Cosmo scope boundary =====================
    G.add_edge(gal_node, "Cosmo: tensions\nCLOSED no mechanism",
               label="scope boundary")
    G.add_edge(gal_node, "Cosmo: H0\nOUTSIDE SCOPE",
               label="predictions")
    G.add_edge(gal_node, "Cosmo: S8\nNEGLIGIBLE",
               label="predictions")
    G.add_edge(gal_node, "Cosmo: DESI w(z)\nINCOMPATIBLE",
               label="predictions")

    # ===================== Problem → tool edges =====================
    G.add_edge("R1: cabibbo\nRESOLVED 0.75σ", "galois\nGL(n,GF(2))")
    G.add_edge("R3: why_n3\nSOLID 13/13", "galois\nGL(n,GF(2))")
    G.add_edge("R4: metric\nCOMPLETE 11/11", "Cadabra2\ntensors")
    G.add_edge("R2: continuum\nOPEN 6-12mo", "Lean 4\nproofs")
    G.add_edge("R6: brannen_sqrt2\nPARTIAL", "Lean 4\nproofs")
    G.add_edge("R5: mass_scaling\nBREAKTHROUGH 7/7", "SymPy\nCAS")
    G.add_edge("R5: mass_scaling\nBREAKTHROUGH 7/7", "SciPy\nODE")
    G.add_edge("R6: brannen_sqrt2\nPARTIAL", "SciPy\nODE")
    G.add_edge("R7: uv_completion\nEMPTY", "SymPy\nCAS")

    return G


# ======================================================================
# 5. CONCEPT FLOW GRAPH — QM emergence chain + cosmological scope
# ======================================================================

def build_concept_flow():
    """Build concept flow: QM emergence chain and cosmological scope boundary."""
    G = nx.DiGraph()

    # --- QM emergence chain ---
    qm_chain = [
        ("Q0: hbar\nfoundations", "#4CAF50"),
        ("Q1: measurement\nuncertainty", "#4CAF50"),
        ("Q2: Born rule\n(in Q0/Q1)", "#9E9E9E"),
        ("Q3: superposition\nNL mixing", "#4CAF50"),
        ("Q4: entanglement\nBell/CHSH", "#4CAF50"),
        ("Q5: spin-1/2\ntopology", "#4CAF50"),
        ("Q6: spin-statistics\nconnection", "#4CAF50"),
        ("Q7: decoherence\nemergent", "#4CAF50"),
    ]
    for name, color in qm_chain:
        G.add_node(name, type="qm", color=color)

    # QM emergence edges
    G.add_edge(qm_chain[0][0], qm_chain[1][0], label="hbar → uncertainty")
    G.add_edge(qm_chain[0][0], qm_chain[2][0], label="contains Born rule")
    G.add_edge(qm_chain[1][0], qm_chain[3][0], label="measurement basis")
    G.add_edge(qm_chain[1][0], qm_chain[4][0], label="Born rule → Bell")
    G.add_edge(qm_chain[3][0], qm_chain[6][0], label="NL → decoherence")
    G.add_edge(qm_chain[4][0], qm_chain[5][0], label="winding number")
    G.add_edge(qm_chain[5][0], qm_chain[6][0])

    # --- Core theory chain ---
    core_nodes = [
        ("R2: continuum\nlimit", "#FF9800"),
        ("R4: metric\nh(Phi)=Phi", "#4CAF50"),
        ("R5: mass k=4\nscaling", "#CDDC39"),
        ("R3: N=3\ngenerations", "#4CAF50"),
        ("R6: Brannen\nsqrt2", "#FF9800"),
        ("R1: Cabibbo\nresolved", "#F44336"),
    ]
    for name, color in core_nodes:
        G.add_node(name, type="core", color=color)

    G.add_edge(core_nodes[0][0], core_nodes[1][0], label="CG theorems")
    G.add_edge(core_nodes[0][0], core_nodes[2][0], label="alpha=2")
    G.add_edge(core_nodes[2][0], core_nodes[3][0], label="K^2")
    G.add_edge(core_nodes[4][0], core_nodes[3][0], label="Koide K=2/3")
    G.add_edge(core_nodes[3][0], core_nodes[5][0], label="CKM")

    # --- Galaxy → Cosmo scope boundary ---
    scope_nodes = [
        ("Galaxy scaling\n46+ scripts", "#4CAF50"),
        ("SCOPE\nBOUNDARY", "#FFD600"),
        ("Cosmo tensions\nCLOSED", "#F44336"),
        ("H0 tension\nOUTSIDE", "#2196F3"),
        ("S8 tension\nNEGLIGIBLE", "#2196F3"),
        ("DESI w(z)\nINCOMPATIBLE", "#F44336"),
    ]
    for name, color in scope_nodes:
        G.add_node(name, type="scope", color=color)

    G.add_edge(core_nodes[1][0], scope_nodes[0][0], label="metric ansatz")
    G.add_edge(scope_nodes[0][0], scope_nodes[1][0])
    G.add_edge(scope_nodes[1][0], scope_nodes[2][0])
    G.add_edge(scope_nodes[1][0], scope_nodes[3][0])
    G.add_edge(scope_nodes[1][0], scope_nodes[4][0])
    G.add_edge(scope_nodes[1][0], scope_nodes[5][0])

    # Connect QM to core via foundations
    G.add_edge(qm_chain[0][0], core_nodes[0][0], label="Z2 substrate")

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

    # 4. Research graph (all 20 folders: R1-R7, Q0-Q7, Galaxy, Cosmo)
    print("\n4. Building research problems graph (R1-R7, Q0-Q7, Galaxy, Cosmo)...")
    research_g = build_research_graph()
    print_stats(research_g, "Research Problems (20 folders)")
    draw_graph(research_g,
               "TGP Research: R1-R7 + Q0-Q7 + Galaxy + Cosmo → Tools",
               "graph_research_problems.png", layout="kamada", figsize=(24, 18))
    export_gexf(research_g, "graph_research_problems.gexf")

    # 5. Concept flow (QM emergence + cosmological scope boundary)
    print("\n5. Building concept flow graph (QM chain + scope boundary)...")
    flow_g = build_concept_flow()
    print_stats(flow_g, "Concept Flow (QM + Core + Cosmo scope)")
    draw_graph(flow_g,
               "TGP Concept Flow: QM Emergence → Core Theory → Galaxy → Cosmo Scope",
               "graph_concept_flow_v2.png", layout="kamada", figsize=(22, 16))
    export_gexf(flow_g, "graph_concept_flow_v2.gexf")

    # Summary
    print("\n" + "=" * 60)
    print("DONE — 5 graphs generated in research/")
    print("=" * 60)
    print("\nFiles:")
    print("  PNG: graph_tex_includes.png, graph_python_modules.png,")
    print("       graph_concept_flow.png, graph_research_problems.png,")
    print("       graph_concept_flow_v2.png")
    print("  GEXF (Gephi): *.gexf")


if __name__ == "__main__":
    main()

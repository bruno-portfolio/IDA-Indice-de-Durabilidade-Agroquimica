import logging
from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd

from .config import LOG_FORMAT, LOG_DATE_FORMAT

logging.basicConfig(format=LOG_FORMAT, datefmt=LOG_DATE_FORMAT, level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class FitnessCostResult:
    position: int
    residue_name: str

    conservation_cost: float
    functional_cost: float
    frac_evidence_cost: float

    fitness_cost_proxy: float

    fixation_probability: str
    interpretation: str

@dataclass
class FRACEvidence:
    position: int
    mutation_id: str
    original_aa: str
    mutant_aa: str
    fixed: bool
    appeared_not_fixed: bool
    frequency: str

FUNCTIONAL_ROLE_COSTS = {
    "catalytic": 1.0,
    "direct_contact": 0.8,
    "structural": 0.5,
    "peripheral": 0.2,
    "unknown": 0.5,
}

def estimate_fitness_cost_proxy(
    conservation_score: float,
    functional_role: str,
    frac_evidence: Optional[FRACEvidence] = None,
    weights: Optional[dict] = None,
) -> FitnessCostResult:
    if weights is None:
        weights = {
            "conservation": 0.4,
            "functional": 0.3,
            "frac": 0.3,
        }

    conservation_cost = conservation_score ** 2

    functional_cost = FUNCTIONAL_ROLE_COSTS.get(functional_role, 0.5)

    if frac_evidence is not None:
        if frac_evidence.fixed:
            frac_cost = 0.2
        elif frac_evidence.appeared_not_fixed:
            frac_cost = 0.9
        else:
            frac_cost = 0.5
    else:
        frac_cost = 0.5

    fitness_cost_proxy = (
        weights["conservation"] * conservation_cost +
        weights["functional"] * functional_cost +
        weights["frac"] * frac_cost
    )

    fixation_prob, interpretation = _interpret_fitness_cost(
        fitness_cost_proxy,
        conservation_score,
        functional_role,
        frac_evidence,
    )

    return FitnessCostResult(
        position=frac_evidence.position if frac_evidence else 0,
        residue_name="",
        conservation_cost=conservation_cost,
        functional_cost=functional_cost,
        frac_evidence_cost=frac_cost,
        fitness_cost_proxy=fitness_cost_proxy,
        fixation_probability=fixation_prob,
        interpretation=interpretation,
    )

def _interpret_fitness_cost(
    fitness_cost: float,
    conservation: float,
    functional_role: str,
    frac_evidence: Optional[FRACEvidence],
) -> tuple[str, str]:
    if fitness_cost >= 0.7:
        fixation_prob = "low"
    elif fitness_cost >= 0.4:
        fixation_prob = "medium"
    else:
        fixation_prob = "high"

    parts = []

    if conservation >= 0.95:
        parts.append("posição extremamente conservada")
    elif conservation >= 0.8:
        parts.append("posição altamente conservada")

    if functional_role == "catalytic":
        parts.append("resíduo catalítico")
    elif functional_role == "direct_contact":
        parts.append("contato direto com ligante")

    if frac_evidence is not None:
        if frac_evidence.fixed:
            parts.append("mutação fixou (custo tolerável)")
        elif frac_evidence.appeared_not_fixed:
            parts.append("mutação apareceu mas não fixou")

    if not parts:
        interpretation = f"Custo de fitness moderado (score={fitness_cost:.2f})"
    else:
        interpretation = "; ".join(parts)

    return fixation_prob, interpretation

def calculate_fitness_cost_batch(
    residue_positions: list[int],
    conservation_scores: dict[int, float],
    functional_roles: dict[int, str],
    frac_data: Optional[dict[int, FRACEvidence]] = None,
    residue_names: Optional[dict[int, str]] = None,
) -> dict[int, FitnessCostResult]:
    if frac_data is None:
        frac_data = {}

    if residue_names is None:
        residue_names = {}

    results = {}

    for pos in residue_positions:
        conservation = conservation_scores.get(pos, 0.5)
        role = functional_roles.get(pos, "unknown")
        frac = frac_data.get(pos)

        result = estimate_fitness_cost_proxy(
            conservation_score=conservation,
            functional_role=role,
            frac_evidence=frac,
        )

        result.position = pos
        result.residue_name = residue_names.get(pos, "UNK")

        results[pos] = result

    logger.info(f"Fitness cost calculado para {len(results)} posições")

    return results

def get_fitness_cost_scores(results: dict[int, FitnessCostResult]) -> dict[int, float]:
    return {pos: r.fitness_cost_proxy for pos, r in results.items()}

def validate_fitness_cost_vs_frac(
    fitness_results: dict[int, FitnessCostResult],
    frac_mutations: list[FRACEvidence],
) -> pd.DataFrame:
    validation_data = []

    for mut in frac_mutations:
        pos = mut.position
        fitness = fitness_results.get(pos)

        if fitness is None:
            continue

        if mut.fixed and mut.frequency == "high":
            consistent = fitness.fitness_cost_proxy < 0.8
            expected = "Fitness cost < 0.8 (mutação fixou)"
        elif mut.appeared_not_fixed:
            consistent = fitness.fitness_cost_proxy > 0.5
            expected = "Fitness cost > 0.5 (não fixou)"
        else:
            consistent = True
            expected = "N/A"

        validation_data.append({
            "position": pos,
            "mutation": f"{mut.original_aa}{pos}{mut.mutant_aa}",
            "fixed": mut.fixed,
            "frequency": mut.frequency,
            "fitness_cost_proxy": fitness.fitness_cost_proxy,
            "expected": expected,
            "consistent": consistent,
        })

    df = pd.DataFrame(validation_data)

    if len(df) > 0:
        n_consistent = df["consistent"].sum()
        logger.info(
            f"Validação fitness cost vs FRAC: {n_consistent}/{len(df)} consistentes"
        )

    return df

def summarize_fitness_costs(
    results: dict[int, FitnessCostResult]
) -> dict:
    if not results:
        return {}

    scores = [r.fitness_cost_proxy for r in results.values()]
    conservation_costs = [r.conservation_cost for r in results.values()]
    functional_costs = [r.functional_cost for r in results.values()]

    fixation_counts = {"low": 0, "medium": 0, "high": 0}
    for r in results.values():
        fixation_counts[r.fixation_probability] += 1

    return {
        "n_positions": len(results),
        "mean_fitness_cost": np.mean(scores),
        "std_fitness_cost": np.std(scores),
        "min_fitness_cost": min(scores),
        "max_fitness_cost": max(scores),
        "mean_conservation_component": np.mean(conservation_costs),
        "mean_functional_component": np.mean(functional_costs),
        "fixation_probability_distribution": fixation_counts,
        "high_cost_positions": [
            pos for pos, r in results.items()
            if r.fitness_cost_proxy >= 0.7
        ],
        "low_cost_positions": [
            pos for pos, r in results.items()
            if r.fitness_cost_proxy < 0.3
        ],
    }

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Estimador de Fitness Cost")
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Executa demonstração com dados sintéticos",
    )

    args = parser.parse_args()

    if args.demo:
        positions = [145, 230, 272, 273]
        conservation_scores = {145: 0.92, 230: 0.95, 272: 0.72, 273: 0.80}
        functional_roles = {
            145: "direct_contact",
            230: "catalytic",
            272: "direct_contact",
            273: "structural"
        }
        residue_names = {145: "TYR", 230: "HIS", 272: "PRO", 273: "ILE"}

        frac_data = {
            272: FRACEvidence(
                position=272,
                mutation_id="FRAC_001",
                original_aa="P",
                mutant_aa="L",
                fixed=True,
                appeared_not_fixed=False,
                frequency="high"
            )
        }

        results = calculate_fitness_cost_batch(
            residue_positions=positions,
            conservation_scores=conservation_scores,
            functional_roles=functional_roles,
            frac_data=frac_data,
            residue_names=residue_names,
        )

        print("\n" + "=" * 60)
        print("ESTIMATIVA DE FITNESS COST")
        print("=" * 60)

        for pos, result in results.items():
            print(f"\nPosição {pos} ({result.residue_name}):")
            print(f"  Conservação cost: {result.conservation_cost:.3f}")
            print(f"  Funcional cost:   {result.functional_cost:.3f}")
            print(f"  FRAC cost:        {result.frac_evidence_cost:.3f}")
            print(f"  FITNESS COST:     {result.fitness_cost_proxy:.3f}")
            print(f"  Prob. fixação:    {result.fixation_probability}")
            print(f"  Interpretação:    {result.interpretation}")

        summary = summarize_fitness_costs(results)
        print("\n" + "-" * 60)
        print("SUMÁRIO:")
        print(f"  Média: {summary['mean_fitness_cost']:.3f}")
        print(f"  Desvio: {summary['std_fitness_cost']:.3f}")
        print(f"  Range: [{summary['min_fitness_cost']:.3f}, {summary['max_fitness_cost']:.3f}]")
        print(f"  Distribuição fixação: {summary['fixation_probability_distribution']}")
        print("=" * 60)
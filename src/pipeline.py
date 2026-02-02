import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from .config import (
    STRUCTURES_DIR, SEQUENCES_DIR, ALIGNMENTS_DIR, RESULTS_DIR,
    LOG_FORMAT, LOG_DATE_FORMAT,
)
from .pdb_parser import PDBParser, Structure, find_binding_site_residues
from .sasa_calculator import calculate_structure_sasa, calculate_burial_score
from .pocket_detector import detect_pockets, calculate_pocket_srs
from .confidence_extractor import load_confidence, analyze_binding_site_confidence
from .msa_runner import run_msa, validate_alignment, map_alignment_to_structure
from .conservation_scorer import calculate_conservation, calculate_ecs_score

from .fitness_cost_estimator import (
    calculate_fitness_cost_batch, get_fitness_cost_scores, FRACEvidence,
)
from .robustness_scanner import (
    calculate_robustness_heuristic, get_robustness_scores,
)
from .ida_molecular import (
    calculate_ida_molecular, IDAMolecularResult, IDAMolecularWeights,
    export_ida_molecular_json, print_ida_molecular_summary,
)
from .ida_regional import (
    calculate_ida_regional, IDARegionalResult,
    calculate_mutation_penalty, calculate_selection_pressure, calculate_target_saturation,
    RegionalMutation, export_ida_regional_results, print_ida_regional_summary,
)
from .ida_programa import (
    calculate_ida_programa, IDAProgramaResult, ManagementProgram, Application,
    export_ida_programa_results, print_ida_programa_summary,
)
from .confidence_reporter import (
    assess_confidence, ConfidenceReport, create_ida_with_confidence,
    print_confidence_summary,
)

logging.basicConfig(format=LOG_FORMAT, datefmt=LOG_DATE_FORMAT, level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class PipelineConfigV2:
    pdb_path: Path
    confidence_json_path: Optional[Path] = None
    sequences_fasta_path: Optional[Path] = None

    binding_site_residues: Optional[list[int]] = None
    binding_site_chain: str = "A"

    ligand_residues: Optional[list[int]] = None
    ligand_chain: Optional[str] = None
    binding_distance: float = 5.0

    query_id: Optional[str] = None
    run_msa: bool = True

    region_id: str = ""
    region_name: str = ""
    season: str = ""
    mutations_csv_path: Optional[Path] = None
    region_data: Optional[dict] = None
    moa_usage_data: Optional[dict] = None

    program: Optional[ManagementProgram] = None

    functional_roles: Optional[dict[int, str]] = None
    frac_positions: Optional[set[int]] = None

    output_dir: Path = RESULTS_DIR

    weights: Optional[IDAMolecularWeights] = None

@dataclass
class PipelineResultV2:
    ida_molecular_result: Optional[IDAMolecularResult] = None
    ida_regional_result: Optional[IDARegionalResult] = None
    ida_programa_result: Optional[IDAProgramaResult] = None

    confidence_report: Optional[ConfidenceReport] = None

    success: bool = False
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    structure: Optional[Structure] = None
    binding_site: Optional[list[int]] = None

    def get_final_ida(self) -> float:
        if self.ida_programa_result:
            return self.ida_programa_result.ida_programa
        elif self.ida_regional_result:
            return self.ida_regional_result.ida_regional
        elif self.ida_molecular_result:
            return self.ida_molecular_result.ida_molecular_aggregate
        return 0.0

class IDAPipelineV2:

    def __init__(self, config: PipelineConfigV2) -> None:
        self.config = config
        self.errors: list[str] = []
        self.warnings: list[str] = []

        self.structure: Optional[Structure] = None
        self.binding_site: list[int] = []
        self.residue_names: dict[int, str] = {}
        self.chain_ids: dict[int, str] = {}

        self.srs_scores: dict[int, float] = {}
        self.sce_scores: dict[int, float] = {}
        self.robustness_scores: dict[int, float] = {}
        self.fitness_cost_scores: dict[int, float] = {}
        self.sce_weights: dict[int, float] = {}
        self.plddt_scores: dict[int, float] = {}
        self.burial_scores: dict[int, float] = {}

        self.structural_metrics: dict = {}
        self.evolutionary_metrics: dict = {}
        self.regional_metrics: dict = {}

    def _log_step(self, step: str) -> None:
        logger.info(f"{'='*20} {step} {'='*20}")

    def _step_load_structure(self) -> bool:
        self._log_step("CARREGAR ESTRUTURA")

        if not self.config.pdb_path.exists():
            self.errors.append(f"Arquivo PDB não encontrado: {self.config.pdb_path}")
            return False

        try:
            parser = PDBParser()
            self.structure = parser.parse(self.config.pdb_path)

            if self.structure.n_residues == 0:
                self.errors.append("Estrutura não contém resíduos")
                return False

            logger.info(
                f"Estrutura carregada: {self.structure.n_chains} cadeias, "
                f"{self.structure.n_residues} resíduos"
            )
            return True

        except Exception as e:
            self.errors.append(f"Erro ao carregar PDB: {e}")
            return False

    def _step_define_binding_site(self) -> bool:
        self._log_step("DEFINIR SÍTIO DE LIGAÇÃO")

        if self.config.binding_site_residues:
            self.binding_site = self.config.binding_site_residues
            logger.info(f"Sítio de ligação fornecido: {len(self.binding_site)} resíduos")

        elif self.config.ligand_residues and self.config.ligand_chain:
            self.binding_site = [
                r.res_seq for r in find_binding_site_residues(
                    self.structure,
                    self.config.ligand_residues,
                    self.config.ligand_chain,
                    self.config.binding_distance,
                )
            ]
            logger.info(f"Sítio calculado a partir do ligante: {len(self.binding_site)} resíduos")

        else:
            logger.info("Detectando cavidades automaticamente...")
            pocket_result = detect_pockets(self.config.pdb_path)

            if pocket_result.success and pocket_result.pockets:
                best_pocket = pocket_result.get_best_pocket()
                if best_pocket:
                    self.binding_site = [r[1] for r in best_pocket.residues]
                    logger.info(f"Melhor cavidade: {len(self.binding_site)} resíduos")

        if not self.binding_site:
            self.errors.append("Não foi possível definir o sítio de ligação")
            return False

        for residue in self.structure.iter_residues():
            if residue.res_seq in self.binding_site:
                self.residue_names[residue.res_seq] = residue.res_name
                self.chain_ids[residue.res_seq] = residue.chain_id

        return True

    def _step_structural_analysis(self) -> bool:
        self._log_step("ANÁLISE ESTRUTURAL")

        try:
            sasa_results = calculate_structure_sasa(self.structure, self.binding_site)
            self.burial_scores = calculate_burial_score(sasa_results, self.binding_site)

            for res_seq in self.binding_site:
                self.srs_scores[res_seq] = self.burial_scores.get(res_seq, 0.5)

            logger.info(f"SASA calculada para {len(sasa_results)} resíduos")

        except Exception as e:
            self.warnings.append(f"Erro no cálculo de SASA: {e}")
            for res_seq in self.binding_site:
                self.srs_scores[res_seq] = 0.5
                self.burial_scores[res_seq] = 0.5

        try:
            pocket_result = detect_pockets(self.config.pdb_path)
            if pocket_result.success:
                best_pocket = pocket_result.get_best_pocket()
                if best_pocket:
                    pocket_srs = calculate_pocket_srs(best_pocket)
                    for res_seq in self.binding_site:
                        base_srs = self.srs_scores.get(res_seq, 0.5)
                        self.srs_scores[res_seq] = (base_srs + pocket_srs) / 2

        except Exception as e:
            self.warnings.append(f"Erro na detecção de cavidades: {e}")

        n_low_conf = 0
        try:
            confidence = load_confidence(
                pdb_path=self.config.pdb_path,
                json_path=self.config.confidence_json_path,
            )

            if confidence:
                for res_seq in self.binding_site:
                    res_conf = confidence.get_residue(res_seq)
                    if res_conf:
                        self.plddt_scores[res_seq] = res_conf.plddt
                        self.sce_weights[res_seq] = res_conf.sce_weight
                        if res_conf.plddt < 70:
                            n_low_conf += 1
                    else:
                        self.plddt_scores[res_seq] = 70.0
                        self.sce_weights[res_seq] = 0.8

                logger.info(f"pLDDT carregado, média={confidence.mean_plddt:.1f}")

                self.structural_metrics = {
                    "mean_plddt": confidence.mean_plddt,
                    "min_plddt": min(self.plddt_scores.values()) if self.plddt_scores else 70,
                    "n_low_confidence_residues": n_low_conf,
                }
            else:
                self.warnings.append("Dados de confiança não disponíveis, usando defaults")
                for res_seq in self.binding_site:
                    self.plddt_scores[res_seq] = 70.0
                    self.sce_weights[res_seq] = 0.8
                self.structural_metrics = {
                    "mean_plddt": 70.0,
                    "min_plddt": 70.0,
                    "n_low_confidence_residues": 0,
                }

        except Exception as e:
            self.warnings.append(f"Erro ao carregar confiança: {e}")
            for res_seq in self.binding_site:
                self.plddt_scores[res_seq] = 70.0
                self.sce_weights[res_seq] = 0.8
            self.structural_metrics = {
                "mean_plddt": 70.0,
                "min_plddt": 70.0,
                "n_low_confidence_residues": 0,
            }

        return True

    def _step_evolutionary_analysis(self) -> bool:
        self._log_step("ANÁLISE EVOLUTIVA")

        sequence = self.structure.get_sequence()
        if not sequence:
            self.warnings.append("Sequência não extraída da estrutura")
            sequence = "X" * max(self.binding_site)

        n_sequences = 0
        alignment_quality = 0.5
        avg_gap_fraction = 0.0

        if self.config.run_msa and self.config.sequences_fasta_path:
            try:
                msa_result = run_msa(self.config.sequences_fasta_path)

                if msa_result.success:
                    validation = validate_alignment(msa_result)
                    if not validation["valid"]:
                        self.warnings.extend(validation.get("issues", []))

                    conservation = calculate_conservation(
                        msa_result,
                        query_id=self.config.query_id,
                    )

                    self.sce_scores = calculate_ecs_score(
                        conservation,
                        self.binding_site,
                    )

                    n_sequences = msa_result.n_sequences
                    alignment_quality = conservation.mean_conservation
                    avg_gap_fraction = validation.get("avg_gap_fraction", 0)

                    logger.info(f"Conservação calculada, média={conservation.mean_conservation:.3f}")
                else:
                    self.warnings.append(f"MSA falhou: {msa_result.error_message}")

            except Exception as e:
                self.warnings.append(f"Erro na análise de MSA: {e}")

        if not self.sce_scores:
            self.warnings.append("Usando valores default para conservação")
            for res_seq in self.binding_site:
                self.sce_scores[res_seq] = 0.7

        self.evolutionary_metrics = {
            "n_sequences": n_sequences,
            "alignment_quality": alignment_quality,
            "avg_gap_fraction": avg_gap_fraction,
        }

        functional_roles = self.config.functional_roles or {
            r: "unknown" for r in self.binding_site
        }

        try:
            fitness_results = calculate_fitness_cost_batch(
                residue_positions=self.binding_site,
                conservation_scores=self.sce_scores,
                functional_roles=functional_roles,
                frac_data=None,
                residue_names=self.residue_names,
            )
            self.fitness_cost_scores = get_fitness_cost_scores(fitness_results)
            logger.info(f"Fitness cost calculado para {len(self.fitness_cost_scores)} resíduos")

        except Exception as e:
            self.warnings.append(f"Erro no cálculo de fitness cost: {e}")
            for res_seq in self.binding_site:
                self.fitness_cost_scores[res_seq] = 0.5

        try:
            robustness_result = calculate_robustness_heuristic(
                binding_site_residues=self.binding_site,
                residue_names=self.residue_names,
                conservation_scores=self.sce_scores,
                burial_scores=self.burial_scores,
            )
            self.robustness_scores = get_robustness_scores(robustness_result)
            logger.info(f"Robustez calculada, média={robustness_result.mean_robustness:.3f}")

        except Exception as e:
            self.warnings.append(f"Erro no cálculo de robustez: {e}")
            for res_seq in self.binding_site:
                self.robustness_scores[res_seq] = 0.5

        return True

    def _step_calculate_ida_molecular(self) -> Optional[IDAMolecularResult]:
        self._log_step("CALCULAR IDA MOLECULAR")

        try:
            result = calculate_ida_molecular(
                binding_site_residues=self.binding_site,
                residue_names=self.residue_names,
                chain_ids=self.chain_ids,
                srs_scores=self.srs_scores,
                sce_scores=self.sce_scores,
                robustness_scores=self.robustness_scores,
                fitness_cost_scores=self.fitness_cost_scores,
                sce_weights=self.sce_weights,
                plddt_scores=self.plddt_scores,
                functional_roles=self.config.functional_roles,
                frac_positions=self.config.frac_positions,
                weights=self.config.weights,
            )

            return result

        except Exception as e:
            self.errors.append(f"Erro no cálculo do IDA Molecular: {e}")
            return None

    def _step_calculate_ida_regional(
        self,
        ida_molecular_result: IDAMolecularResult,
    ) -> Optional[IDARegionalResult]:
        self._log_step("CALCULAR IDA REGIONAL")

        if not self.config.region_id:
            logger.info("Dados regionais não fornecidos, pulando IDA Regional")
            return None

        try:
            ida_per_residue = ida_molecular_result.get_per_residue_scores()

            mutations = []
            if self.config.mutations_csv_path and self.config.mutations_csv_path.exists():
                pass

            mutation_penalty = calculate_mutation_penalty(mutations, ida_per_residue)

            region_data = self.config.region_data or {
                "region_id": self.config.region_id,
                "region_name": self.config.region_name,
                "season": self.config.season,
            }
            selection_pressure = calculate_selection_pressure(region_data)

            target_saturation = calculate_target_saturation(
                self.config.moa_usage_data or {}
            )

            self.regional_metrics = {
                "n_mutations_documented": mutation_penalty.n_mutations,
                "years_of_monitoring": region_data.get("sdhi_years_use", 0),
                "n_data_sources": 1 if self.config.mutations_csv_path else 0,
            }

            result = calculate_ida_regional(
                ida_molecular=ida_molecular_result.ida_molecular_aggregate,
                mutation_penalty=mutation_penalty,
                selection_pressure=selection_pressure,
                target_saturation=target_saturation,
            )

            return result

        except Exception as e:
            self.warnings.append(f"Erro no cálculo do IDA Regional: {e}")
            return None

    def _step_calculate_ida_programa(
        self,
        ida_regional_result: IDARegionalResult,
    ) -> Optional[IDAProgramaResult]:
        self._log_step("CALCULAR IDA PROGRAMA")

        if not self.config.program:
            logger.info("Programa de manejo não fornecido, pulando IDA Programa")
            return None

        try:
            result = calculate_ida_programa(
                ida_regional_result=ida_regional_result,
                program=self.config.program,
            )

            return result

        except Exception as e:
            self.warnings.append(f"Erro no cálculo do IDA Programa: {e}")
            return None

    def run(self) -> PipelineResultV2:
        logger.info("Iniciando pipeline IDA v2.0...")

        if not self._step_load_structure():
            return PipelineResultV2(
                success=False,
                errors=self.errors,
                warnings=self.warnings,
            )

        if not self._step_define_binding_site():
            return PipelineResultV2(
                success=False,
                errors=self.errors,
                warnings=self.warnings,
                structure=self.structure,
            )

        self._step_structural_analysis()

        self._step_evolutionary_analysis()

        ida_molecular = self._step_calculate_ida_molecular()

        if not ida_molecular:
            return PipelineResultV2(
                success=False,
                errors=self.errors,
                warnings=self.warnings,
                structure=self.structure,
                binding_site=self.binding_site,
            )

        ida_regional = None
        if self.config.region_id:
            ida_regional = self._step_calculate_ida_regional(ida_molecular)

        ida_programa = None
        if ida_regional and self.config.program:
            ida_programa = self._step_calculate_ida_programa(ida_regional)

        confidence = assess_confidence(
            self.structural_metrics,
            self.evolutionary_metrics,
            self.regional_metrics,
        )

        output_dir = self.config.output_dir
        output_dir.mkdir(parents=True, exist_ok=True)

        export_ida_molecular_json(ida_molecular, output_dir / "ida_molecular.json")

        if ida_regional:
            export_ida_regional_results([ida_regional], output_dir / "ida_regional.json")

        if ida_programa:
            export_ida_programa_results([ida_programa], output_dir / "ida_programa.json")

        final_result = {
            "ida_molecular": ida_molecular.ida_molecular_aggregate,
            "ida_regional": ida_regional.ida_regional if ida_regional else None,
            "ida_programa": ida_programa.ida_programa if ida_programa else None,
            "confidence": confidence.to_dict(),
        }

        with open(output_dir / "ida_complete.json", "w", encoding="utf-8") as f:
            json.dump(final_result, f, indent=2, ensure_ascii=False)

        logger.info("Pipeline v2.0 concluído com sucesso!")

        return PipelineResultV2(
            ida_molecular_result=ida_molecular,
            ida_regional_result=ida_regional,
            ida_programa_result=ida_programa,
            confidence_report=confidence,
            success=True,
            errors=self.errors,
            warnings=self.warnings,
            structure=self.structure,
            binding_site=self.binding_site,
        )

def run_pipeline_v2(
    pdb_path: Path,
    binding_site_residues: Optional[list[int]] = None,
    sequences_fasta: Optional[Path] = None,
    confidence_json: Optional[Path] = None,
    region_id: str = "",
    region_name: str = "",
    season: str = "",
    region_data: Optional[dict] = None,
    program: Optional[ManagementProgram] = None,
    output_dir: Optional[Path] = None,
) -> PipelineResultV2:
    config = PipelineConfigV2(
        pdb_path=Path(pdb_path),
        binding_site_residues=binding_site_residues,
        sequences_fasta_path=Path(sequences_fasta) if sequences_fasta else None,
        confidence_json_path=Path(confidence_json) if confidence_json else None,
        region_id=region_id,
        region_name=region_name,
        season=season,
        region_data=region_data,
        program=program,
        output_dir=Path(output_dir) if output_dir else RESULTS_DIR,
    )

    pipeline = IDAPipelineV2(config)
    return pipeline.run()

def run_demo() -> PipelineResultV2:
    from .config import EXAMPLE_DATA_DIR

    example_pdb = EXAMPLE_DATA_DIR / "estrutura_exemplo.pdb"
    example_fasta = EXAMPLE_DATA_DIR / "homologos_exemplo.fasta"
    example_confidence = EXAMPLE_DATA_DIR / "confidence_exemplo.json"

    if not example_pdb.exists():
        logger.error("Arquivo de exemplo não encontrado. Execute com dados reais.")
        return PipelineResultV2(
            success=False,
            errors=["Arquivo de exemplo não encontrado"],
        )

    demo_program = ManagementProgram(
        program_id="DEMO_001",
        name="Programa Demo - Rotação Balanceada",
        applications=[
            Application("V8", "multisite", "Mancozeb", None),
            Application("R1", "SDHI", "Fluxapiroxade", "Piraclostrobina"),
            Application("R3", "DMI", "Protioconazol", "Trifloxistrobina"),
            Application("R5.1", "multisite", "Clorotalonil", None),
        ]
    )

    return run_pipeline_v2(
        pdb_path=example_pdb,
        binding_site_residues=[145, 146, 147, 230, 231, 272, 273],
        sequences_fasta=example_fasta if example_fasta.exists() else None,
        confidence_json=example_confidence if example_confidence.exists() else None,
        region_id="BR-MT",
        region_name="Mato Grosso",
        season="2024/25",
        region_data={
            "soy_area_ha": 12_000_000,
            "sdhi_years_use": 10,
            "disease_severity": "high",
        },
        program=demo_program,
    )

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Pipeline IDA v2.0")
    parser.add_argument("pdb_file", type=Path, nargs="?", help="Arquivo PDB")
    parser.add_argument(
        "-r", "--residues",
        type=int,
        nargs="+",
        help="Resíduos do sítio de ligação",
    )
    parser.add_argument(
        "-s", "--sequences",
        type=Path,
        help="Arquivo FASTA com sequências homólogas",
    )
    parser.add_argument(
        "-c", "--confidence",
        type=Path,
        help="Arquivo JSON de confiança AlphaFold",
    )
    parser.add_argument(
        "--region",
        type=str,
        default="",
        help="ID da região (ex: BR-MT)",
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=RESULTS_DIR,
        help="Diretório de saída",
    )
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Executa demonstração com dados sintéticos",
    )

    args = parser.parse_args()

    if args.demo:
        result = run_demo()
    elif args.pdb_file:
        result = run_pipeline_v2(
            pdb_path=args.pdb_file,
            binding_site_residues=args.residues,
            sequences_fasta=args.sequences,
            confidence_json=args.confidence,
            region_id=args.region,
            output_dir=args.output,
        )
    else:
        parser.print_help()
        exit(1)

    if result.success:
        print("\n" + "=" * 70)
        print("RESULTADO DO PIPELINE IDA v2.0")
        print("=" * 70)

        if result.ida_molecular_result:
            print_ida_molecular_summary(result.ida_molecular_result)

        if result.ida_regional_result:
            print_ida_regional_summary(result.ida_regional_result)

        if result.ida_programa_result:
            print_ida_programa_summary(result.ida_programa_result)

        if result.confidence_report:
            print_confidence_summary(result.confidence_report)
    else:
        print("\nPipeline falhou:")
        for error in result.errors:
            print(f"  ERRO: {error}")

    if result.warnings:
        print("\nAvisos:")
        for warning in result.warnings:
            print(f"  AVISO: {warning}")
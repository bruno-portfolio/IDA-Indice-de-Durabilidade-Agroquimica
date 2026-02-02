__version__ = "2.0.0"
__author__ = "Bruno"

from .config import (
    PROJECT_ROOT,
    DATA_DIR,
    STRUCTURES_DIR,
    SEQUENCES_DIR,
    RESULTS_DIR,
    ALPHAFOLD_CONFIG,
    UNIPROT_CONFIG,
    SIDRA_CONFIG,
    SDH_TARGET,
    ANALYSIS_PARAMS,
)

from .alphafold_client import AlphaFoldClient, AlphaFoldStructure
from .uniprot_client import UniProtClient, FastaSequence, SearchResult
from .sidra_client import SidraClient, ProductivityData, TerritorialLevel

from .pdb_parser import PDBParser, Structure, Residue, Atom
from .sasa_calculator import calculate_structure_sasa, SASAResult
from .pocket_detector import detect_pockets, Pocket, PocketDetectionResult
from .confidence_extractor import load_confidence, ConfidenceData

from .msa_runner import run_msa, MSAResult
from .conservation_scorer import calculate_conservation, ConservationResult

from .fitness_cost_estimator import (
    estimate_fitness_cost_proxy,
    calculate_fitness_cost_batch,
    FitnessCostResult,
    FRACEvidence,
)
from .robustness_scanner import (
    calculate_robustness_heuristic,
    RobustnessResult,
    ResidueRobustness,
)

from .ida_molecular import (
    calculate_ida_molecular,
    IDAMolecularResult,
    IDAMolecularWeights,
    export_ida_molecular_json,
    print_ida_molecular_summary,
)
from .ida_regional import (
    calculate_ida_regional,
    IDARegionalResult,
    RegionalMutation,
    calculate_mutation_penalty,
    calculate_selection_pressure,
    calculate_target_saturation,
    export_ida_regional_results,
    print_ida_regional_summary,
)
from .ida_programa import (
    calculate_ida_programa,
    IDAProgramaResult,
    ManagementProgram,
    Application,
    export_ida_programa_results,
    print_ida_programa_summary,
)

from .confidence_reporter import (
    assess_confidence,
    ConfidenceReport,
    create_ida_with_confidence,
    print_confidence_summary,
)

from .pipeline import (
    IDAPipelineV2,
    PipelineConfigV2,
    PipelineResultV2,
    run_pipeline_v2,
)

from .stability_scorer import calculate_stability_scores, StabilityResult
from .ida_calculator import (
    calculate_ida,
    IDAResult,
    IDAWeights,
    export_ida_json,
    export_ida_csv,
    print_ida_summary,
)

__all__ = [
    "__version__",
    "__author__",
    "PROJECT_ROOT",
    "DATA_DIR",
    "STRUCTURES_DIR",
    "SEQUENCES_DIR",
    "RESULTS_DIR",
    "SDH_TARGET",
    "ANALYSIS_PARAMS",
    "AlphaFoldClient",
    "UniProtClient",
    "SidraClient",
    "PDBParser",
    "Structure",
    "calculate_structure_sasa",
    "detect_pockets",
    "load_confidence",
    "run_msa",
    "calculate_conservation",
    "estimate_fitness_cost_proxy",
    "calculate_fitness_cost_batch",
    "FitnessCostResult",
    "calculate_robustness_heuristic",
    "RobustnessResult",
    "calculate_ida_molecular",
    "IDAMolecularResult",
    "IDAMolecularWeights",
    "calculate_ida_regional",
    "IDARegionalResult",
    "RegionalMutation",
    "calculate_ida_programa",
    "IDAProgramaResult",
    "ManagementProgram",
    "Application",
    "assess_confidence",
    "ConfidenceReport",
    "create_ida_with_confidence",
    "IDAPipelineV2",
    "PipelineConfigV2",
    "PipelineResultV2",
    "run_pipeline_v2",
    "calculate_stability_scores",
    "calculate_ida",
    "IDAResult",
    "IDAWeights",
    "export_ida_json",
    "print_ida_summary",
]
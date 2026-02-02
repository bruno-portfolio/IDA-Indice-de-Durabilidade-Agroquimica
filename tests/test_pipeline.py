"""Testes para o pipeline principal."""

import pytest
from pathlib import Path
import json


class TestPipelineImports:
    """Testa se os módulos podem ser importados."""

    def test_import_pipeline(self) -> None:
        """Pipeline deve ser importável."""
        from src.pipeline import run_pipeline_v2, PipelineConfigV2
        assert run_pipeline_v2 is not None
        assert PipelineConfigV2 is not None

    def test_import_pdb_parser(self) -> None:
        """PDB parser deve ser importável."""
        from src.pdb_parser import PDBParser, Structure
        assert PDBParser is not None
        assert Structure is not None

    def test_import_ida_molecular(self) -> None:
        """IDA Molecular deve ser importável."""
        from src.ida_molecular import calculate_ida_molecular
        assert calculate_ida_molecular is not None

    def test_import_sasa_calculator(self) -> None:
        """SASA calculator deve ser importável."""
        from src.sasa_calculator import calculate_structure_sasa
        assert calculate_structure_sasa is not None

    def test_import_confidence_extractor(self) -> None:
        """Confidence extractor deve ser importável."""
        from src.confidence_extractor import load_confidence_json
        assert load_confidence_json is not None


class TestPDBParser:
    """Testes para o parser de PDB."""

    def test_parse_pdb_sdh(self, example_pdb: Path) -> None:
        """Deve parsear PDB do SDH corretamente."""
        if not example_pdb.exists():
            pytest.skip("PDB de exemplo não encontrado")

        from src.pdb_parser import PDBParser

        parser = PDBParser()
        structure = parser.parse(example_pdb)

        assert structure is not None
        assert len(structure.chains) >= 1
        assert structure.n_residues > 0
        assert structure.n_atoms > 0

    def test_parse_pdb_returns_structure(self, example_pdb: Path) -> None:
        """parse_pdb deve retornar objeto Structure."""
        if not example_pdb.exists():
            pytest.skip("PDB de exemplo não encontrado")

        from src.pdb_parser import PDBParser, Structure

        parser = PDBParser()
        structure = parser.parse(example_pdb)
        assert isinstance(structure, Structure)


class TestSASACalculator:
    """Testes para o calculador de SASA."""

    def test_calculate_sasa(self, example_pdb: Path, sdh_binding_site: list) -> None:
        """Deve calcular SASA para resíduos do binding site."""
        if not example_pdb.exists():
            pytest.skip("PDB de exemplo não encontrado")

        from src.pdb_parser import PDBParser
        from src.sasa_calculator import calculate_structure_sasa

        parser = PDBParser()
        structure = parser.parse(example_pdb)
        sasa_results = calculate_structure_sasa(structure, sdh_binding_site)

        assert sasa_results is not None
        assert len(sasa_results) > 0


class TestConfidenceExtractor:
    """Testes para extração de confiança."""

    def test_load_confidence_json(self, example_confidence: Path) -> None:
        """Deve carregar JSON de confiança."""
        if not example_confidence.exists():
            pytest.skip("JSON de confiança não encontrado")

        from src.confidence_extractor import load_confidence_json

        confidence = load_confidence_json(example_confidence)

        assert confidence is not None
        assert len(confidence.residue_scores) > 0

    def test_confidence_plddt_range(self, example_confidence: Path) -> None:
        """pLDDT deve estar no range 0-100."""
        if not example_confidence.exists():
            pytest.skip("JSON de confiança não encontrado")

        from src.confidence_extractor import load_confidence_json

        confidence = load_confidence_json(example_confidence)

        for score in confidence.residue_scores:
            assert 0 <= score.plddt <= 100


class TestPipelineIntegration:
    """Testes de integração do pipeline."""

    def test_pipeline_sdh(
        self,
        example_pdb: Path,
        example_fasta: Path,
        example_confidence: Path,
        sdh_binding_site: list,
    ) -> None:
        """Pipeline deve rodar para SDH."""
        if not all(p.exists() for p in [example_pdb, example_fasta, example_confidence]):
            pytest.skip("Arquivos de exemplo não encontrados")

        from src.pipeline import run_pipeline_v2

        results = run_pipeline_v2(
            pdb_path=example_pdb,
            binding_site_residues=sdh_binding_site,
            sequences_fasta=example_fasta,
            confidence_json=example_confidence,
        )

        assert results is not None
        assert results.success or len(results.errors) == 0
        assert results.ida_molecular_result is not None

    def test_pipeline_ida_range(
        self,
        example_pdb: Path,
        example_fasta: Path,
        example_confidence: Path,
        sdh_binding_site: list,
    ) -> None:
        """IDA deve estar no range 0-1."""
        if not all(p.exists() for p in [example_pdb, example_fasta, example_confidence]):
            pytest.skip("Arquivos de exemplo não encontrados")

        from src.pipeline import run_pipeline_v2

        results = run_pipeline_v2(
            pdb_path=example_pdb,
            binding_site_residues=sdh_binding_site,
            sequences_fasta=example_fasta,
            confidence_json=example_confidence,
        )

        if results.ida_molecular_result:
            ida = results.ida_molecular_result.ida_molecular_aggregate
            assert 0 <= ida <= 1


class TestResultsFiles:
    """Testes para arquivos de resultados."""

    def test_sdh_results_exist(self, project_root: Path) -> None:
        """Resultados do SDH devem existir."""
        results_file = project_root / "resultados" / "ida_molecular_sdh.json"
        if not results_file.exists():
            pytest.skip("Resultados do SDH não encontrados")

        data = json.loads(results_file.read_text(encoding="utf-8"))

        assert "metadata" in data
        assert "aggregate_score" in data
        assert data["metadata"]["target"] == "Succinate dehydrogenase (SDH-B)"

    def test_cyp51_results_exist(self, project_root: Path) -> None:
        """Resultados do CYP51 devem existir."""
        results_file = project_root / "resultados" / "ida_molecular_cyp51.json"
        if not results_file.exists():
            pytest.skip("Resultados do CYP51 não encontrados")

        data = json.loads(results_file.read_text(encoding="utf-8"))

        assert "metadata" in data
        assert "aggregate_score" in data
        assert "CYP51" in data["metadata"]["target"]

    def test_cytb_results_exist(self, project_root: Path) -> None:
        """Resultados do Cyt b devem existir."""
        results_file = project_root / "resultados" / "ida_molecular_cytb.json"
        if not results_file.exists():
            pytest.skip("Resultados do Cyt b não encontrados")

        data = json.loads(results_file.read_text(encoding="utf-8"))

        assert "metadata" in data
        assert "aggregate_score" in data
        assert "Cytochrome" in data["metadata"]["target"]

    def test_results_ida_valid(self, project_root: Path) -> None:
        """IDA nos resultados deve ser válido."""
        for target in ["sdh", "cyp51", "cytb"]:
            results_file = project_root / "resultados" / f"ida_molecular_{target}.json"
            if not results_file.exists():
                continue

            data = json.loads(results_file.read_text(encoding="utf-8"))
            ida = data["aggregate_score"]["ida_molecular"]

            assert 0 <= ida <= 1, f"IDA inválido para {target}: {ida}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

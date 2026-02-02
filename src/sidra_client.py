import logging
from dataclasses import dataclass
from enum import IntEnum
from pathlib import Path
from typing import Optional

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from .config import SIDRA_CONFIG, FIELD_DATA_DIR, LOG_FORMAT, LOG_DATE_FORMAT

logging.basicConfig(format=LOG_FORMAT, datefmt=LOG_DATE_FORMAT, level=logging.INFO)
logger = logging.getLogger(__name__)

class TerritorialLevel(IntEnum):
    BRASIL = 1
    REGIAO = 2
    UF = 3
    MESORREGIAO = 4
    MICRORREGIAO = 6
    MUNICIPIO = 7

@dataclass
class ProductivityData:
    year: int
    territorial_code: str
    territorial_name: str
    yield_kg_ha: Optional[float]
    area_ha: Optional[float]
    production_ton: Optional[float]

    @property
    def yield_ton_ha(self) -> Optional[float]:
        return self.yield_kg_ha / 1000 if self.yield_kg_ha else None

    @property
    def safra(self) -> str:
        return f"{self.year - 1}/{str(self.year)[-2:]}"

@dataclass
class SidraResult:
    data: list[ProductivityData]
    success: bool
    query_url: str
    error_message: Optional[str] = None

class SidraClient:

    def __init__(
        self,
        output_dir: Path = FIELD_DATA_DIR,
        timeout: int = SIDRA_CONFIG.timeout_seconds,
    ) -> None:
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.timeout = timeout
        self.session = self._create_session()

    def _create_session(self) -> requests.Session:
        session = requests.Session()
        retry_strategy = Retry(
            total=3,
            backoff_factor=2,
            status_forcelist=[429, 500, 502, 503, 504],
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        session.mount("https://", adapter)
        return session

    def _build_url(
        self,
        table_id: int,
        territorial_level: TerritorialLevel,
        periods: str,
        variables: list[int],
        product_code: int,
    ) -> str:
        base = SIDRA_CONFIG.base_url
        vars_str = ",".join(str(v) for v in variables)

        url = f"{base}/t/{table_id}/n{territorial_level}/all/v/{vars_str}/p/{periods}/c81/{product_code}"

        return url

    def _parse_response(self, data: list[dict]) -> list[ProductivityData]:
        results = []

        if len(data) <= 1:
            return results

        for row in data[1:]:
            try:
                value_str = row.get("V", "")
                value = float(value_str) if value_str and value_str != "-" else None

                var_code = row.get("D2C")
                yield_val = value if var_code == "112" else None
                area_val = value if var_code == "109" else None
                prod_val = value if var_code == "214" else None

                year_str = row.get("D3N", "")
                year = int(year_str) if year_str.isdigit() else 0

                if year > 0:
                    results.append(ProductivityData(
                        year=year,
                        territorial_code=row.get("D1C", ""),
                        territorial_name=row.get("D1N", ""),
                        yield_kg_ha=yield_val,
                        area_ha=area_val,
                        production_ton=prod_val,
                    ))

            except (ValueError, KeyError) as e:
                logger.debug(f"Erro ao processar linha: {e}")
                continue

        return results

    def get_soy_productivity(
        self,
        territorial_level: TerritorialLevel = TerritorialLevel.UF,
        periods: str = "last 15",
        include_area: bool = False,
        include_production: bool = False,
    ) -> SidraResult:
        variables = [SIDRA_CONFIG.yield_variable]
        if include_area:
            variables.append(109)
        if include_production:
            variables.append(214)

        url = self._build_url(
            table_id=SIDRA_CONFIG.table_pam,
            territorial_level=territorial_level,
            periods=periods,
            variables=variables,
            product_code=SIDRA_CONFIG.soy_product_code,
        )

        logger.info(f"Consultando SIDRA: {url}")

        try:
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()

            data = response.json()
            productivity_data = self._parse_response(data)

            logger.info(f"Obtidos {len(productivity_data)} registros")

            return SidraResult(
                data=productivity_data,
                success=True,
                query_url=url,
            )

        except requests.exceptions.HTTPError as e:
            logger.error(f"Erro HTTP: {e}")
            return SidraResult(
                data=[],
                success=False,
                query_url=url,
                error_message=f"HTTP {e.response.status_code}",
            )
        except requests.exceptions.Timeout:
            logger.error("Timeout na consulta SIDRA")
            return SidraResult(
                data=[],
                success=False,
                query_url=url,
                error_message="Timeout",
            )
        except requests.exceptions.JSONDecodeError:
            logger.error("Erro ao decodificar resposta JSON")
            return SidraResult(
                data=[],
                success=False,
                query_url=url,
                error_message="JSON inválido",
            )

    def save_csv(
        self,
        data: list[ProductivityData],
        filename: str,
        output_dir: Optional[Path] = None,
    ) -> Path:
        output_dir = output_dir or self.output_dir
        output_path = output_dir / filename

        if not filename.endswith(".csv"):
            output_path = output_dir / f"{filename}.csv"

        lines = ["ano,safra,cod_territorial,nome_territorial,rendimento_kg_ha,rendimento_t_ha"]

        for d in data:
            yield_t = f"{d.yield_ton_ha:.3f}" if d.yield_ton_ha else ""
            yield_kg = f"{d.yield_kg_ha:.1f}" if d.yield_kg_ha else ""

            line = f"{d.year},{d.safra},{d.territorial_code},{d.territorial_name},{yield_kg},{yield_t}"
            lines.append(line)

        output_path.write_text("\n".join(lines), encoding="utf-8")
        logger.info(f"Salvo: {output_path}")

        return output_path

def filter_productivity_data(
    data: list[ProductivityData],
    min_yield: float = 500.0,
    max_yield: float = 6000.0,
) -> list[ProductivityData]:
    filtered = []

    for d in data:
        if d.yield_kg_ha is None:
            continue
        if d.yield_kg_ha < min_yield or d.yield_kg_ha > max_yield:
            continue
        filtered.append(d)

    removed = len(data) - len(filtered)
    if removed > 0:
        logger.info(f"Filtrados {removed} registros (outliers ou nulos)")

    return filtered

def calculate_productivity_stats(data: list[ProductivityData]) -> dict:
    if not data:
        return {"count": 0}

    yields = [d.yield_kg_ha for d in data if d.yield_kg_ha is not None]

    if not yields:
        return {"count": len(data), "valid_yields": 0}

    years = sorted(set(d.year for d in data))
    territories = set(d.territorial_name for d in data)

    return {
        "count": len(data),
        "valid_yields": len(yields),
        "min_yield_kg_ha": min(yields),
        "max_yield_kg_ha": max(yields),
        "avg_yield_kg_ha": sum(yields) / len(yields),
        "avg_yield_t_ha": sum(yields) / len(yields) / 1000,
        "year_range": f"{min(years)}-{max(years)}",
        "n_territories": len(territories),
    }

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Coleta dados de produtividade de soja do SIDRA/IBGE"
    )
    parser.add_argument(
        "-l", "--level",
        choices=["uf", "micro", "meso"],
        default="uf",
        help="Nível territorial",
    )
    parser.add_argument(
        "-p", "--periods",
        default="last 15",
        help="Períodos (ex: 'last 15', '2010-2023')",
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=FIELD_DATA_DIR,
        help="Diretório de saída",
    )

    args = parser.parse_args()

    level_map = {
        "uf": TerritorialLevel.UF,
        "micro": TerritorialLevel.MICRORREGIAO,
        "meso": TerritorialLevel.MESORREGIAO,
    }

    client = SidraClient(output_dir=args.output)
    result = client.get_soy_productivity(
        territorial_level=level_map[args.level],
        periods=args.periods,
    )

    if result.success:
        filtered = filter_productivity_data(result.data)
        output_path = client.save_csv(filtered, f"soja_produtividade_{args.level}.csv")

        stats = calculate_productivity_stats(filtered)
        print(f"\nEstatísticas:")
        for key, value in stats.items():
            if isinstance(value, float):
                print(f"  {key}: {value:.2f}")
            else:
                print(f"  {key}: {value}")
    else:
        print(f"Erro: {result.error_message}")
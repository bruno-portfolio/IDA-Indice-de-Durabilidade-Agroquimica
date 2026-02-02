import json
from pathlib import Path
from datetime import datetime

import streamlit as st
import pandas as pd
import plotly.graph_objects as go

PROJECT_ROOT = Path(__file__).parent
RESULTS_DIR = PROJECT_ROOT / "resultados"

VERSION = "2.4"

st.set_page_config(
    page_title="IDA - Índice de Durabilidade Agroquímica",
    page_icon="🌿",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown("""
<style>
    .risk-indicator {
        display: inline-block;
        width: 12px;
        height: 12px;
        border-radius: 50%;
        margin-right: 8px;
    }
    .risk-high { background: #722f37; }
    .risk-elevated { background: #fd7e14; }
    .risk-moderate { background: #ffc107; }
    .risk-low { background: #28a745; }
    .context-bar {
        background: linear-gradient(90deg, #1a1a2e 0%, #16213e 100%);
        padding: 12px 20px;
        border-radius: 6px;
        margin-bottom: 20px;
        border-left: 3px solid #0f3460;
        font-size: 13px;
        color: #e8e8e8;
    }
    .metric-card {
        background: #f8f9fa;
        padding: 20px;
        border-radius: 8px;
        border-left: 4px solid;
        margin: 8px 0;
    }
    .section-header {
        color: #1a1a2e;
        border-bottom: 2px solid #0f3460;
        padding-bottom: 8px;
        margin-bottom: 16px;
    }
</style>
""", unsafe_allow_html=True)

RISK_CONFIG = {
    "high": {"color": "#722f37", "bg": "#f5e6e8", "label": "Alto", "css": "risk-high"},
    "moderate_high": {"color": "#fd7e14", "bg": "#ffe0b3", "label": "Elevado", "css": "risk-elevated"},
    "moderate": {"color": "#ffc107", "bg": "#fff3cd", "label": "Moderado", "css": "risk-moderate"},
    "low_moderate": {"color": "#90EE90", "bg": "#e8f5e9", "label": "Moderado-Baixo", "css": "risk-low"},
    "low": {"color": "#28a745", "bg": "#d4edda", "label": "Baixo", "css": "risk-low"},
}

REGIOES = {
    "Brasil - Centro-Oeste": {"pressure": 0.75, "mutations": 2, "saturation": 0.45},
    "Brasil - Sul": {"pressure": 0.82, "mutations": 1, "saturation": 0.38},
    "Brasil - Sudeste": {"pressure": 0.78, "mutations": 1, "saturation": 0.40},
    "Brasil - Norte": {"pressure": 0.90, "mutations": 0, "saturation": 0.25},
    "Brasil - Nordeste": {"pressure": 0.88, "mutations": 0, "saturation": 0.30},
    "Brasil - MATOPIBA": {"pressure": 0.70, "mutations": 2, "saturation": 0.50},
}

WEIGHTS = {"srs": 0.30, "sce": 0.35, "robustness": 0.35}

FUNCTIONAL_ROLES_PT = {
    "catalytic": "Catalítico",
    "direct_contact": "Contato direto",
    "structural": "Estrutural",
    "peripheral": "Periférico",
    "unknown": "—",
}

def get_risk_from_score(score: float) -> str:
    if score >= 0.7:
        return "low"
    elif score >= 0.5:
        return "low_moderate"
    elif score >= 0.3:
        return "moderate"
    elif score >= 0.2:
        return "moderate_high"
    return "high"

def get_barrier_label(score: float) -> str:
    if score >= 0.7:
        return "Barreira alta"
    elif score >= 0.5:
        return "Barreira moderada-alta"
    elif score >= 0.3:
        return "Barreira moderada"
    elif score >= 0.2:
        return "Barreira baixa-moderada"
    return "Barreira baixa"

def create_gauge_chart(value: float, title: str) -> go.Figure:
    risk = get_risk_from_score(value)
    color = RISK_CONFIG[risk]["color"]

    fig = go.Figure(go.Indicator(
        mode="gauge+number",
        value=value,
        domain={"x": [0, 1], "y": [0, 1]},
        title={"text": title, "font": {"size": 18, "color": "#fff", "family": "Inter, sans-serif"}},
        number={"font": {"size": 42, "color": "#fff", "family": "Inter, sans-serif"}, "valueformat": ".2f"},
        gauge={
            "axis": {"range": [0, 1], "tickwidth": 1, "tickcolor": "#555", "tickvals": [0, 0.3, 0.5, 0.7, 1.0], "tickfont": {"color": "#888"}},
            "bar": {"color": color, "thickness": 0.75},
            "bgcolor": "#2a2a3e",
            "borderwidth": 1,
            "bordercolor": "#3a3a4e",
            "steps": [
                {"range": [0, 0.3], "color": "#3d2a2e"},
                {"range": [0.3, 0.5], "color": "#3d3528"},
                {"range": [0.5, 0.7], "color": "#3d3d28"},
                {"range": [0.7, 1.0], "color": "#2a3d2e"},
            ],
        }
    ))
    fig.update_layout(
        height=240,
        margin=dict(l=20, r=20, t=70, b=10),
        paper_bgcolor="#1a1a2e",
        plot_bgcolor="#1a1a2e",
        font={"family": "Inter, sans-serif"}
    )
    return fig

def load_json_safe(path: Path) -> dict | None:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None

def get_demo_molecular_data() -> dict:
    return {
        "metadata": {
            "target": "Succinato Desidrogenase (SDH)",
            "organism": "Phakopsora pachyrhizi",
            "defensivo_class": "Carboxamidas (SDHI)",
            "frac_code": "7",
            "date_generated": datetime.now().isoformat(),
            "version": VERSION,
        },
        "aggregate_score": {
            "ida_molecular": 0.42,
            "interpretation": {
                "barrier_level": "moderate_barrier",
                "risk_level": "moderate",
            }
        },
        "component_summary": {
            "avg_srs": 0.47,
            "avg_sce": 0.54,
            "avg_robustness": 0.50,
            "avg_fitness_cost": 0.66
        },
        "confidence": {
            "avg_plddt": 88.6,
            "min_plddt": 82.7,
            "n_residues": 6,
            "confidence_level": "alta"
        },
        "residue_details": [
            {"res_seq": 216, "res_name": "HIS", "srs": 0.00, "sce": 0.85, "robustness": 1.00, "ida_molecular": 0.65, "functional_role": "catalytic", "frac_documented": True},
            {"res_seq": 218, "res_name": "ILE", "srs": 0.40, "sce": 0.65, "robustness": 0.53, "ida_molecular": 0.52, "functional_role": "direct_contact", "frac_documented": True},
            {"res_seq": 169, "res_name": "PRO", "srs": 0.72, "sce": 0.70, "robustness": 0.09, "ida_molecular": 0.48, "functional_role": "direct_contact", "frac_documented": True},
            {"res_seq": 170, "res_name": "SER", "srs": 1.00, "sce": 0.45, "robustness": 0.00, "ida_molecular": 0.16, "functional_role": "peripheral", "frac_documented": False},
            {"res_seq": 172, "res_name": "TRP", "srs": 0.62, "sce": 0.80, "robustness": 0.18, "ida_molecular": 0.51, "functional_role": "structural", "frac_documented": False},
            {"res_seq": 173, "res_name": "TRP", "srs": 0.03, "sce": 0.82, "robustness": 0.96, "ida_molecular": 0.62, "functional_role": "structural", "frac_documented": False},
        ]
    }

def get_regional_data(ida_molecular: float, region: str) -> dict:
    reg_config = REGIOES.get(region, REGIOES["Brasil - Centro-Oeste"])
    pressure = reg_config["pressure"]
    mutations = reg_config["mutations"]
    saturation = reg_config["saturation"]

    mutation_factor = 1.0 - (mutations * 0.05)
    saturation_factor = 1.0 - (saturation * 0.15)
    ida_regional = ida_molecular * pressure * mutation_factor * saturation_factor

    if pressure < 0.8 or mutations >= 2:
        status = "Atenção"
    elif pressure < 0.85:
        status = "Monitorar"
    else:
        status = "Favorável"

    reasons = []
    if pressure < 0.8:
        reasons.append("pressão de seleção alta")
    if mutations >= 2:
        reasons.append(f"{mutations} mutações detectadas")
    if saturation > 0.4:
        reasons.append("saturação de MoA elevada")

    return {
        "ida_molecular_base": ida_molecular,
        "ida_regional": ida_regional,
        "pressure_factor": pressure,
        "region": region,
        "factors": {
            "selection_pressure": {
                "value": "alta" if pressure < 0.8 else "moderada" if pressure < 0.9 else "baixa",
                "score": pressure,
                "reason": f"Fator de pressão: {pressure:.0%}"
            },
            "mutations_detected": {
                "value": mutations,
                "score": mutation_factor,
                "reason": f"{mutations} mutação(ões) documentada(s) no alvo"
            },
            "moa_saturation": {
                "value": f"{saturation:.0%}",
                "score": saturation_factor,
                "reason": f"{saturation:.0%} das aplicações usam este MoA"
            },
        },
        "status": status,
        "status_reason": " + ".join(reasons) if reasons else "Indicadores dentro do esperado",
    }

def get_program_data(ida_regional: float) -> dict:
    return {
        "ida_regional_base": ida_regional,
        "programs": [
            {
                "name": "Rotação Completa",
                "ida_programa": min(ida_regional * 1.25, 1.0),
                "bonus": +0.25,
                "effect": "protetor",
                "sequence": ["SDHI", "DMI", "Multi", "QoI", "Multi"],
                "details": {"Rotação": +0.12, "Multissítios": +0.08, "Sem repetição": +0.05},
                "issues": [],
            },
            {
                "name": "Uso Moderado",
                "ida_programa": ida_regional,
                "bonus": 0.0,
                "effect": "neutro",
                "sequence": ["SDHI", "DMI", "SDHI", "QoI"],
                "details": {"Rotação parcial": +0.05, "Repetição SDHI": -0.05},
                "issues": ["SDHI usado 2x"],
            },
            {
                "name": "Risco Elevado",
                "ida_programa": ida_regional * 0.70,
                "bonus": -0.30,
                "effect": "acelerador",
                "sequence": ["SDHI", "SDHI", "SDHI", "DMI"],
                "details": {"Repetição sequencial": -0.20, "Sem multissítio": -0.10},
                "issues": ["SDHI 3x sequencial", "Sem multissítio"],
            },
        ]
    }

def show_context_header(mol_data: dict, use_demo: bool, selected_region: str):
    meta = mol_data.get("metadata", {})
    source = "Demo" if use_demo else "Pipeline"
    date_str = meta.get("date_generated", datetime.now().isoformat())[:10]

    st.markdown(
        f"""
        <div class="context-bar">
            <b>🎯 Alvo:</b> {meta.get('target', 'SDH')} &nbsp;&nbsp;│&nbsp;&nbsp;
            <b>🦠 Patógeno:</b> <i>{meta.get('organism', 'P. pachyrhizi')}</i> &nbsp;&nbsp;│&nbsp;&nbsp;
            <b>💊 MoA:</b> {meta.get('defensivo_class', 'SDHI')} (FRAC {meta.get('frac_code', '7')}) &nbsp;&nbsp;│&nbsp;&nbsp;
            <b>📍 Região:</b> {selected_region} &nbsp;&nbsp;│&nbsp;&nbsp;
            <b>Fonte:</b> {source} &nbsp;&nbsp;│&nbsp;&nbsp;
            {date_str} &nbsp;│&nbsp; v{VERSION}
        </div>
        """,
        unsafe_allow_html=True
    )

def show_risk_indicator(score: float) -> str:
    risk = get_risk_from_score(score)
    css = RISK_CONFIG[risk]["css"]
    label = RISK_CONFIG[risk]["label"]
    return f'<span class="risk-indicator {css}"></span>{label}'

def show_score_card(label: str, score: float, subtitle: str = ""):
    risk = get_risk_from_score(score)
    config = RISK_CONFIG[risk]
    barrier = get_barrier_label(score)

    st.markdown(
        f"""
        <div style="background: {config['bg']}; padding: 16px; border-radius: 6px;
                    border-left: 4px solid {config['color']}; margin: 4px 0;">
            <div style="font-size: 11px; color: #666; text-transform: uppercase; letter-spacing: 0.5px;">{label}</div>
            <div style="font-size: 32px; font-weight: 600; color: #1a1a2e; margin: 4px 0;">{score:.2f}</div>
            <div style="font-size: 12px; color: #555;">
                <span class="risk-indicator {config['css']}"></span>
                {barrier} · Risco {config['label']}
            </div>
            {f'<div style="font-size: 11px; color: #888; margin-top: 4px;">{subtitle}</div>' if subtitle else ''}
        </div>
        """,
        unsafe_allow_html=True
    )

def get_action_recommendation(ida_score: float) -> dict:
    if ida_score >= 0.7:
        return {
            "nivel": "✓ Baixo Risco",
            "cor": "#28a745",
            "acao": "Monitoramento padrão",
            "detalhes": [
                "Manter programa atual de aplicação",
                "Monitorar eficácia a cada safra",
                "Sem necessidade de mudanças imediatas"
            ]
        }
    elif ida_score >= 0.5:
        return {
            "nivel": "◐ Risco Moderado",
            "cor": "#d4a017",
            "acao": "Atenção ao manejo",
            "detalhes": [
                "Implementar rotação de mecanismos de ação",
                "Incluir multissítios no programa",
                "Aumentar frequência de monitoramento"
            ]
        }
    elif ida_score >= 0.3:
        return {
            "nivel": "◑ Risco Elevado",
            "cor": "#e65100",
            "acao": "Rotação obrigatória",
            "detalhes": [
                "Rotação de MoA é obrigatória",
                "Máximo 2 aplicações do mesmo grupo/safra",
                "Usar misturas com multissítios",
                "Monitorar resistência ativamente"
            ]
        }
    else:
        return {
            "nivel": "✕ Alto Risco",
            "cor": "#722f37",
            "acao": "Revisão urgente",
            "detalhes": [
                "Revisar programa imediatamente",
                "Limitar uso deste MoA a 1x/safra",
                "Priorizar multissítios e biológicos",
                "Consultar assistência técnica"
            ]
        }

def identify_limiting_factors(comp: dict) -> list[tuple[str, float, str]]:
    factors = [
        ("Restrição Estrutural", comp.get("avg_srs", 0.5)),
        ("Conservação Evolutiva", comp.get("avg_sce", 0.5)),
        ("Robustez do Pocket", comp.get("avg_robustness", 0.5)),
    ]
    limiting = [(label, val, "baixo" if val < 0.5 else "moderado")
                for label, val in sorted(factors, key=lambda x: x[1])
                if val < 0.7]
    return limiting[:2]

st.sidebar.markdown("## 🌿 IDA")
st.sidebar.caption(f"Índice de Durabilidade Agroquímica · v{VERSION}")
st.sidebar.divider()

page = st.sidebar.radio(
    "Navegação",
    ["🏠 Visão Geral", "🔬 Análise Molecular", "🗺️ Análise Regional", "📋 Comparar Programas", "📖 Documentação"],
    label_visibility="collapsed"
)

st.sidebar.divider()

st.sidebar.markdown("##### Configuração")

data_source = st.sidebar.radio(
    "Fonte de dados",
    ["Demonstração", "Pipeline"],
    index=0,
    help="Demonstração: dados simulados. Pipeline: dados reais processados."
)
use_demo = data_source == "Demonstração"

selected_region = st.sidebar.selectbox(
    "Região",
    list(REGIOES.keys()),
    index=0,
    help="Selecione a região para ajuste do IDA Regional"
)

FITOPATOGENOS = {
    "Ferrugem da Soja (Phakopsora) - SDHI": "phakopsora_sdh",
    "Brusone (Pyricularia) - QoI": "pyricularia_cytb",
    "Septoriose do Trigo (Zymoseptoria) - DMI": "zymoseptoria_cyp51",
    "Antracnose (Colletotrichum) - SDHI": "colletotrichum_sdh",
}

selected_pathogen = st.sidebar.selectbox(
    "Patógeno/Alvo",
    list(FITOPATOGENOS.keys()),
    index=0,
)

if use_demo:
    mol_data = get_demo_molecular_data()
else:
    pathogen_file = FITOPATOGENOS.get(selected_pathogen, "phakopsora_sdh")
    mol_data = load_json_safe(RESULTS_DIR / f"ida_{pathogen_file}.json")
    if mol_data is None:
        mol_data = get_demo_molecular_data()

ida_mol = mol_data.get("aggregate_score", {}).get("ida_molecular", 0.5)
reg_data = get_regional_data(ida_mol, selected_region)
prog_data = get_program_data(reg_data["ida_regional"])

if page == "🏠 Visão Geral":
    st.markdown("## IDA · Índice de Durabilidade Agroquímica")
    st.caption("Ferramenta de apoio à decisão para manejo de resistência a fungicidas")

    show_context_header(mol_data, use_demo, selected_region)

    ida_reg = reg_data.get("ida_regional", ida_mol * 0.85)
    programs = prog_data.get("programs", [])
    best_prog = max(programs, key=lambda x: x.get("ida_programa", 0)) if programs else {"ida_programa": ida_reg}
    ida_prog = best_prog.get("ida_programa", ida_reg)

    col1, col2, col3 = st.columns(3)

    with col1:
        st.plotly_chart(create_gauge_chart(ida_mol, "IDA Molecular"), width="stretch")
        st.caption("Barreira intrínseca do alvo")

    with col2:
        st.plotly_chart(create_gauge_chart(ida_reg, "IDA Regional"), width="stretch")
        st.caption(f"Ajustado para {selected_region}")

    with col3:
        st.plotly_chart(create_gauge_chart(ida_prog, "IDA Programa"), width="stretch")
        st.caption("Melhor cenário de manejo")

    rec = get_action_recommendation(ida_reg)

    st.markdown(
        f"""
        <div style="background: {rec['cor']}; padding: 20px; border-radius: 8px; margin: 16px 0;">
            <div style="font-size: 18px; font-weight: 600; color: #fff; margin-bottom: 8px;">
                {rec['nivel']}
            </div>
            <div style="font-size: 14px; color: #fff; margin-bottom: 12px;">
                <b>Ação recomendada:</b> {rec['acao']}
            </div>
            <ul style="margin: 0; padding-left: 20px; color: rgba(255,255,255,0.9); font-size: 13px;">
                {"".join(f"<li>{d}</li>" for d in rec["detalhes"])}
            </ul>
        </div>
        """,
        unsafe_allow_html=True
    )

    with st.expander("Como interpretar os IDAs"):
        st.markdown("""
        | IDA | Significado |
        |-----|-------------|
        | **Molecular** | Barreira intrínseca do alvo contra mutações de escape |
        | **Regional** | Ajuste pela pressão de seleção local (uso, mutações detectadas) |
        | **Programa** | Efeito do programa de manejo na durabilidade |

        **Escala:** 0.7-1.0 = Baixo risco · 0.5-0.7 = Moderado · 0.3-0.5 = Elevado · <0.3 = Alto
        """)

elif page == "🔬 Análise Molecular":
    st.markdown("## Análise Molecular")
    show_context_header(mol_data, use_demo, selected_region)

    agg = mol_data.get("aggregate_score", {})
    comp = mol_data.get("component_summary", {})
    conf = mol_data.get("confidence", {})
    residues = mol_data.get("residue_details", [])

    col_main, col_side = st.columns([2, 1])

    with col_main:
        show_score_card("IDA Molecular", ida_mol, "Barreira evolutiva intrínseca")

    with col_side:
        plddt = conf.get("avg_plddt", 0)
        qual = "Muito alta" if plddt >= 90 else "Alta" if plddt >= 70 else "Moderada" if plddt >= 50 else "Baixa"

        st.markdown("**Qualidade dos dados**")
        st.markdown(f"Modelo estrutural: **{qual}** (pLDDT {plddt:.0f})")
        st.markdown(f"Posições avaliadas: **{conf.get('n_residues', 0)}**")

    st.divider()

    st.markdown("#### Composição do Score")
    st.caption(f"IDA = SRS × {WEIGHTS['srs']:.0%} + SCE × {WEIGHTS['sce']:.0%} + Robustez × {WEIGHTS['robustness']:.0%}")

    col1, col2, col3 = st.columns(3)

    with col1:
        srs = comp.get("avg_srs", 0.5)
        st.metric("Restrição Estrutural (SRS)", f"{srs:.2f}")
        st.progress(srs)
        st.caption("Espaço físico para mutação")

    with col2:
        sce = comp.get("avg_sce", 0.5)
        st.metric("Conservação Evolutiva (SCE)", f"{sce:.2f}")
        st.progress(sce)
        st.caption("Importância funcional")

    with col3:
        rob = comp.get("avg_robustness", 0.5)
        st.metric("Robustez do Pocket", f"{rob:.2f}")
        st.progress(rob)
        st.caption("Tolerância a mudanças")

    limiting = identify_limiting_factors(comp)
    if limiting:
        factors_text = ", ".join([f"{label} ({value:.2f})" for label, value, _ in limiting])
        st.warning(f"**Fatores limitantes:** {factors_text}")

    if residues:
        st.divider()
        st.markdown("#### Detalhes por Posição")

        with st.expander("Ver tabela de posições"):
            df = pd.DataFrame(residues)
            df["Função"] = df["functional_role"].map(FUNCTIONAL_ROLES_PT)
            df["FRAC"] = df["frac_documented"].apply(lambda x: "Sim" if x else "—")

            st.dataframe(
                df[["res_seq", "res_name", "srs", "sce", "robustness", "ida_molecular", "Função", "FRAC"]],
                hide_index=True,
                column_config={
                    "res_seq": "Posição",
                    "res_name": "AA",
                    "srs": st.column_config.NumberColumn("SRS", format="%.2f"),
                    "sce": st.column_config.NumberColumn("SCE", format="%.2f"),
                    "robustness": st.column_config.NumberColumn("Robustez", format="%.2f"),
                    "ida_molecular": st.column_config.NumberColumn("IDA", format="%.2f"),
                }
            )

        frac_count = sum(1 for r in residues if r.get("frac_documented", False))
        if frac_count > 0:
            st.error(f"**{frac_count} posição(ões) com mutações FRAC documentadas** — requer atenção especial no manejo.")

elif page == "🗺️ Análise Regional":
    st.markdown("## Análise Regional")
    show_context_header(mol_data, use_demo, selected_region)

    ida_mol_base = reg_data.get("ida_molecular_base", 0.5)
    ida_reg = reg_data.get("ida_regional", ida_mol_base)

    col1, col2 = st.columns(2)

    with col1:
        show_score_card("IDA Molecular (Base)", ida_mol_base)

    with col2:
        show_score_card("IDA Regional", ida_reg, f"Região: {selected_region}")

    delta = ida_reg - ida_mol_base
    delta_pct = (delta / ida_mol_base) * 100 if ida_mol_base > 0 else 0

    if delta < 0:
        st.markdown(f"**Redução de {abs(delta_pct):.0f}%** devido à pressão regional")

    st.divider()

    status = reg_data.get("status", "—")
    reason = reg_data.get("status_reason", "")

    status_colors = {"Atenção": "#fff3cd", "Alerta": "#f8d7da", "Monitorar": "#e3f2fd", "Favorável": "#e8f5e9"}

    st.markdown(
        f"""
        <div style="background: {status_colors.get(status, '#f5f5f5')}; padding: 16px; border-radius: 6px;">
            <div style="font-size: 16px; font-weight: 600;">Status: {status}</div>
            <div style="font-size: 13px; color: #555; margin-top: 4px;">{reason}</div>
        </div>
        """,
        unsafe_allow_html=True
    )

    st.markdown("#### Fatores de Pressão")

    factors = reg_data.get("factors", {})
    labels = {
        "selection_pressure": "Pressão de Seleção",
        "mutations_detected": "Mutações Detectadas",
        "moa_saturation": "Saturação do MoA",
    }

    for key, data in factors.items():
        col_f1, col_f2 = st.columns([3, 1])
        with col_f1:
            st.markdown(f"**{labels.get(key, key)}**")
            st.caption(data.get("reason", ""))
        with col_f2:
            score = data.get("score", 1.0)
            color = "#dc3545" if score < 0.8 else "#ffc107" if score < 0.95 else "#28a745"
            st.markdown(f'<span style="color:{color}; font-weight:600;">●</span> {score:.0%}', unsafe_allow_html=True)

elif page == "📋 Comparar Programas":
    st.markdown("## Comparar Programas")
    show_context_header(mol_data, use_demo, selected_region)

    ida_base = prog_data.get("ida_regional_base", 0.35)
    programs = prog_data.get("programs", [])

    st.markdown(f"**IDA Regional (base):** {ida_base:.2f}")
    st.divider()

    cols = st.columns(len(programs))

    effect_styles = {
        "protetor": {"bg": "#e8f5e9", "border": "#4caf50", "label": "Protetor"},
        "neutro": {"bg": "#fff8e1", "border": "#ffc107", "label": "Neutro"},
        "acelerador": {"bg": "#ffebee", "border": "#f44336", "label": "Acelerador"},
    }

    for i, prog in enumerate(programs):
        with cols[i]:
            effect = prog.get("effect", "neutro")
            style = effect_styles.get(effect, effect_styles["neutro"])
            ida_prog = prog.get("ida_programa", ida_base)
            bonus = prog.get("bonus", 0)

            st.markdown(
                f"""
                <div style="background: {style['bg']}; padding: 16px; border-radius: 8px;
                            border-top: 3px solid {style['border']}; min-height: 240px;">
                    <div style="font-size: 14px; font-weight: 600; margin-bottom: 8px; color: #000;">
                        {prog.get('name', f'Programa {i+1}')}
                    </div>
                    <div style="font-size: 28px; font-weight: 700; color: #000;">{ida_prog:.2f}</div>
                    <div style="font-size: 12px; color: #333; margin-bottom: 12px;">
                        {style['label']} ({bonus:+.0%})
                    </div>
                    <div style="font-size: 11px; color: #000;">
                        <b>Sequência:</b><br>
                        {' → '.join(prog.get('sequence', []))}
                    </div>
                </div>
                """,
                unsafe_allow_html=True
            )

            for issue in prog.get("issues", []):
                st.caption(f"▲ {issue}")

    with st.expander("Ver detalhamento de bônus/penalidades"):
        for prog in programs:
            st.markdown(f"**{prog.get('name')}**")
            for key, val in prog.get("details", {}).items():
                color = "green" if val > 0 else "red" if val < 0 else "gray"
                st.markdown(f"- {key}: <span style='color:{color}'>{val:+.0%}</span>", unsafe_allow_html=True)
            st.markdown("---")

elif page == "📖 Documentação":
    st.markdown("## Documentação")
    show_context_header(mol_data, use_demo, selected_region)

    st.markdown("### O que é o IDA?")
    st.markdown("""
    O Índice de Durabilidade Agroquímica (IDA) estima o risco de perda de eficácia
    de um fungicida devido à evolução de resistência.

    - **Não é** uma previsão de vida útil em anos
    - **É** uma estimativa comparativa para apoiar decisões de manejo
    """)

    st.markdown("### Escala de Risco")
    st.markdown("""
    | Score | Barreira | Risco | Ação |
    |-------|----------|-------|------|
    | 0.7 - 1.0 | Alta | Baixo | Monitoramento padrão |
    | 0.5 - 0.7 | Moderada-Alta | Moderado | Atenção ao manejo |
    | 0.3 - 0.5 | Moderada | Elevado | Rotação obrigatória |
    | 0.0 - 0.3 | Baixa | Alto | Revisão urgente |
    """)

    st.markdown("### Fontes de Dados")
    st.markdown("""
    | Fonte | Informação |
    |-------|------------|
    | PDB / AlphaFold | Estrutura 3D do alvo |
    | UniProt | Sequências homólogas |
    | FRAC | Mutações de resistência |
    | Dados regionais | Pressão de seleção local |
    """)

    st.markdown("### Limitações")
    st.markdown("""
    - Baseado em modelos estruturais (AlphaFold/PDB)
    - Dados regionais podem ter cobertura limitada
    - Não substitui monitoramento de campo
    """)

    st.markdown("### Glossário")
    st.markdown("""
    | Termo | Significado |
    |-------|-------------|
    | MoA | Mecanismo de Ação |
    | SDHI | Inibidor da Succinato Desidrogenase |
    | DMI | Inibidor da Desmetilação |
    | QoI | Inibidor do citocromo bc1 |
    | FRAC | Fungicide Resistance Action Committee |
    | pLDDT | Confiança estrutural (0-100) |
    """)

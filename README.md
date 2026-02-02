![CI](https://github.com/bruno-portfolio/IDA-Indice-de-Durabilidade-Agroquimica/actions/workflows/ci.yml/badge.svg) ![Python](https://img.shields.io/badge/Python-3.10+-3776AB?logo=python&logoColor=white) ![Streamlit](https://img.shields.io/badge/Streamlit-1.29+-FF4B4B?logo=streamlit&logoColor=white) ![BioPython](https://img.shields.io/badge/BioPython-1.81+-3776AB?logo=python&logoColor=white) ![Plotly](https://img.shields.io/badge/Plotly-5.18+-3F4F75?logo=plotly&logoColor=white) ![Pandas](https://img.shields.io/badge/Pandas-2.0+-150458?logo=pandas&logoColor=white) ![NumPy](https://img.shields.io/badge/NumPy-1.24+-013243?logo=numpy&logoColor=white) ![License](https://img.shields.io/badge/License-MIT-green)
# IDA - Índice de Durabilidade Agroquímica

## O que é o IDA?

O IDA é uma ferramenta que estima o **risco de perda de eficácia** de fungicidas devido à evolução de resistência.

Quando um fungicida é usado repetidamente, os fungos podem desenvolver resistência. O IDA **não prevê quantos anos o produto vai durar**, mas mede a **barreira evolutiva** contra esse processo, ajudando você a:

- Escolher programas de aplicação mais duráveis
- Decidir quando rotacionar mecanismos de ação
- Priorizar fungicidas com menor risco de resistência

<img width="1402" height="890" alt="image" src="https://github.com/user-attachments/assets/b096eaac-93f6-46a9-a52e-00c5e1cf82cb" />
<img width="1400" height="687" alt="image" src="https://github.com/user-attachments/assets/a8445277-d635-432f-84e2-d42853b541d5" />
<img width="1349" height="674" alt="image" src="https://github.com/user-attachments/assets/0112dd63-0a0b-456a-bea6-649d68955a54" />
<img width="1263" height="769" alt="image" src="https://github.com/user-attachments/assets/d79cee84-c6e5-4aa6-b493-ae0f5b82dde8" />
<img width="1344" height="868" alt="image" src="https://github.com/user-attachments/assets/8671e4b2-b075-439a-bfe7-715069534833" />

## Instalação

### Pré-requisitos

- Python 3.10 ou superior
- pip (gerenciador de pacotes Python)

### Passo a passo

```bash
# 1. Clone ou baixe o projeto
cd pasta_do_projeto/ida_projeto

# 2. Crie um ambiente virtual (recomendado)
python -m venv venv

# 3. Ative o ambiente virtual
# Windows:
venv\Scripts\activate
# Linux/Mac:
source venv/bin/activate

# 4. Instale as dependências
pip install -r requirements.txt
```

### Dependências principais

- `biopython` - Análise de sequências e estruturas
- `numpy` - Cálculos numéricos
- `pandas` - Manipulação de dados
- `streamlit` - Interface web
- `plotly` - Gráficos interativos

---

## Como usar

### Opção 1: Interface Web (mais fácil)

```bash
cd ida_projeto
streamlit run app.py
```

Acesse no navegador: `http://localhost:8501`

A interface permite:
- Ver os resultados de forma visual
- Comparar programas de manejo
- Entender o significado de cada score

### Opção 2: Linha de comando

```bash
cd ida_projeto
python -m src.pipeline
```

---

## Como analisar um novo fungo/fungicida

### O que você precisa

Para analisar um novo alvo, você precisa de **3 arquivos**:

| Arquivo | O que é | Onde conseguir |
|---------|---------|----------------|
| **Estrutura 3D** (.pdb) | Estrutura da proteína alvo | [PDB](https://www.rcsb.org/) ou [AlphaFold](https://alphafold.ebi.ac.uk/) |
| **Sequências homólogas** (.fasta) | Proteínas similares de outros organismos | [UniProt](https://www.uniprot.org/) |
| **Confiança** (.json) | Qualidade da estrutura (opcional para PDB experimental) | AlphaFold ou criar manualmente |

### Passo 1: Obter a estrutura 3D

#### Opção A: Estrutura experimental (PDB) - RECOMENDADO

1. Acesse [RCSB PDB](https://www.rcsb.org/)
2. Busque pelo nome do alvo (ex: "succinate dehydrogenase fungicide")
3. Baixe o arquivo `.pdb`
4. Salve em `dados/brutos/estruturas/`

**Dica:** Estruturas experimentais são mais confiáveis. Procure estruturas com:
- Resolução < 2.5 Å
- Ligante (fungicida) presente na estrutura

#### Opção B: Estrutura predita (AlphaFold)

1. Acesse [AlphaFold DB](https://alphafold.ebi.ac.uk/)
2. Busque pelo código UniProt do organismo
3. Baixe o arquivo `.pdb` e o `.json` de confiança
4. Salve ambos em `dados/brutos/estruturas/`

### Passo 2: Obter sequências homólogas

1. Acesse [UniProt](https://www.uniprot.org/)
2. Busque pela proteína alvo (ex: "SDH fungi")
3. Use o filtro "Reviewed" para qualidade
4. Selecione 50-100 sequências de diferentes espécies
5. Download em formato FASTA
6. Salve em `dados/brutos/sequencias/`

**Exemplo de busca UniProt:**
```
(protein_name:"succinate dehydrogenase") AND (taxonomy_id:4751) AND (reviewed:true)
```
- `taxonomy_id:4751` = Fungi
- `reviewed:true` = Sequências curadas

### Passo 3: Criar arquivo de confiança (para estruturas PDB)

Se você usou uma estrutura do PDB (não do AlphaFold), crie um arquivo de confiança:

```python
# Script: criar_confianca.py
import json
from pathlib import Path

# Ler o PDB para contar resíduos
pdb_path = Path("dados/brutos/estruturas/SUA_ESTRUTURA.pdb")
residues = set()

with open(pdb_path) as f:
    for line in f:
        if line.startswith("ATOM"):
            res_num = int(line[22:26].strip())
            residues.add(res_num)

max_res = max(residues)

# Criar confiança alta (estruturas experimentais são confiáveis)
confidence = {
    "source": "Estrutura experimental PDB - alta confiança",
    "confidenceScore": [95.0 if i in residues else 0 for i in range(1, max_res + 1)]
}

# Salvar
out_path = Path("dados/brutos/estruturas/SUA_ESTRUTURA_confidence.json")
with open(out_path, "w") as f:
    json.dump(confidence, f, indent=2)

print(f"Arquivo de confiança criado: {out_path}")
```

### Passo 4: Identificar o sítio de ligação

O sítio de ligação são os resíduos da proteína que interagem com o fungicida.

#### Se a estrutura tem o ligante (fungicida):

Use um visualizador molecular (PyMOL, ChimeraX) para identificar os resíduos próximos ao ligante (< 5 Å).

#### Se não tem o ligante:

Consulte a literatura ou use ferramentas de predição de cavidades.

**Exemplo para SDH (SDHI):**
Os resíduos típicos do sítio de ligação são: 169, 170, 172, 173, 216, 218 (na subunidade B).

### Passo 5: Executar o pipeline

```python
from src.pipeline import run_pipeline_v2
from pathlib import Path

# Configurar caminhos
pdb_file = Path("dados/brutos/estruturas/SUA_ESTRUTURA.pdb")
fasta_file = Path("dados/brutos/sequencias/SEUS_HOMOLOGOS.fasta")
confidence_file = Path("dados/brutos/estruturas/SUA_ESTRUTURA_confidence.json")

# Resíduos do sítio de ligação (ajuste para seu alvo!)
binding_site = [169, 170, 172, 173, 216, 218]

# Executar
results = run_pipeline_v2(
    pdb_path=pdb_file,
    binding_site_residues=binding_site,
    sequences_fasta=fasta_file,
    confidence_json=confidence_file,
)

# Ver resultado
print(f"IDA Molecular: {results.ida_molecular_result.ida_molecular_aggregate:.2f}")
print(f"Confiança: {results.ida_molecular_result.confidence_level}")
```

### Passo 6: Visualizar resultados

```bash
streamlit run app.py
```

Selecione "Dados do Pipeline" na barra lateral para ver seus resultados.

---

## Alvos suportados

O IDA foi desenvolvido para os principais alvos de fungicidas:

| Alvo | Classe de Fungicida | Exemplos |
|------|---------------------|----------|
| **SDH** (Succinato desidrogenase) | SDHI / Carboxamidas | Boscalida, Fluxapiroxade |
| **CYP51** (Esterol 14α-demetilase) | DMI / Triazóis | Tebuconazol, Propiconazol |
| **Cyt b** (Citocromo bc1) | QoI / Estrobilurinas | Azoxistrobina, Piraclostrobina |

---

## Estrutura do projeto

```
ida_projeto/
├── app.py                 # Interface Streamlit
├── requirements.txt       # Dependências
├── README.md              # Este arquivo
│
├── src/                   # Código fonte
│   ├── pipeline.py        # Pipeline principal
│   ├── ida_molecular.py   # Cálculo IDA Molecular
│   ├── ida_regional.py    # Cálculo IDA Regional
│   ├── ida_programa.py    # Cálculo IDA Programa
│   ├── pdb_parser.py      # Leitura de estruturas PDB
│   ├── sasa_calculator.py # Cálculo de área superficial
│   └── ...
│
├── dados/                 # Dados de entrada
│   ├── brutos/
│   │   ├── estruturas/    # Arquivos PDB
│   │   └── sequencias/    # Arquivos FASTA
│   └── exemplo/           # Dados de demonstração
│
└── resultados/            # Resultados do pipeline
    ├── ida_molecular.json
    ├── ida_regional.json
    └── ida_programa.json
```

---

## Interpretação dos resultados

### Score IDA

| Score | Risco | Significado | O que fazer |
|-------|-------|-------------|-------------|
| **0.7 - 1.0** | Baixo | Barreira evolutiva ALTA | Monitorar normalmente |
| **0.5 - 0.7** | Moderado | Barreira MODERADA | Rotacionar MoA |
| **0.3 - 0.5** | Elevado | Barreira BAIXA | Rotação obrigatória |
| **0.0 - 0.3** | Alto | Barreira MUITO BAIXA | Revisar programa |

**IMPORTANTE:** IDA ALTO = BOM (maior proteção contra resistência)

### Níveis de confiança

| Nível | Significado |
|-------|-------------|
| **Muito alta** | Resultados muito confiáveis |
| **Alta** | Resultados confiáveis |
| **Média** | Resultados indicativos |
| **Baixa** | Interpretar com cautela |

---

## Perguntas frequentes

### O IDA prevê quantos anos o fungicida vai durar?

**Não.** O IDA mede o **risco relativo** de perda de eficácia, não o tempo exato. Use-o para comparar opções e tomar decisões de manejo.

### Posso usar para qualquer fungo?

Sim, desde que você tenha a estrutura 3D da proteína alvo e sequências homólogas. O sistema foi validado para fungos fitopatogênicos.

### Os dados de demonstração são reais?

Os dados de demonstração são **simulados** para ilustrar o funcionamento. Para análise real, execute o pipeline com seus próprios dados.

### Preciso instalar softwares adicionais?

O IDA funciona com as dependências Python básicas. Para melhores resultados, recomendamos:
- **MAFFT** - Alinhamento múltiplo (opcional, melhora SCE)
- **fpocket** - Detecção de cavidades (opcional)

---

## Exemplo completo: Analisando ferrugem asiática da soja

```bash
# 1. Estrutura: baixar SDH do PDB
# Usamos PDB 2FBW (SDH com carboxin ligado)

# 2. Sequências: baixar do UniProt
# Busca: "(protein_name:SDH) AND (taxonomy:Fungi)"

# 3. Executar pipeline
cd ida_projeto
python -c "
from src.pipeline import run_pipeline_v2
from pathlib import Path

results = run_pipeline_v2(
    pdb_path=Path('dados/brutos/estruturas/2FBW_chainB.pdb'),
    binding_site_residues=[169, 170, 172, 173, 216, 218],
    sequences_fasta=Path('dados/brutos/sequencias/sdh_b_fungi_homologs.fasta'),
    confidence_json=Path('dados/brutos/estruturas/2FBW_chainB_confidence.json'),
)

print(f'IDA: {results.ida_molecular_result.ida_molecular_aggregate:.2f}')
"

# 4. Visualizar
streamlit run app.py
```

## Suporte

Este é um projeto de pesquisa em desenvolvimento.
brunoescalhao@gmail.com

# [Licensa](https://github.com/bruno-portfolio/IDA-Indice-de-Durabilidade-Agroquimica/blob/main/LICENSE)

---

## Glossário

| Termo | Significado |
|-------|-------------|
| **MoA** | Mecanismo de Ação - como o fungicida mata o fungo |
| **SDHI** | Inibidor da Succinato Desidrogenase |
| **DMI** | Inibidor da Desmetilação (triazóis) |
| **QoI** | Inibidor do citocromo bc1 (estrobilurinas) |
| **FRAC** | Fungicide Resistance Action Committee |
| **Multissítio** | Fungicida que atua em múltiplos alvos (ex: mancozeb) |
| **pLDDT** | Métrica de confiança do AlphaFold (0-100) |
| **PDB** | Protein Data Bank - banco de estruturas 3D |
| **UniProt** | Banco de dados de sequências de proteínas |

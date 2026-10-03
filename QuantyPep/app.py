# =====================================
# Troplia — FastAPI (Site + API)
# Idiomas completados: de, it, ja, zh, ru
# =====================================

from fastapi import FastAPI, Request, Response
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from datetime import datetime
from io import StringIO
import csv
import os
import re
import urllib.parse
import secrets

app = FastAPI(title="Troplia")

# Monta /static
static_dir = os.path.join(os.getcwd(), "static")
os.makedirs(static_dir, exist_ok=True)
app.mount("/static", StaticFiles(directory=static_dir), name="static")

HISTORY: list[dict] = []  # (legacy; not used)
# Per-user in-memory history, keyed by a cookie 'troplia_id'
HISTORY_BY_USER: dict[str, list[dict]] = {}

def _get_sid_from_request(request: Request) -> str:
    return request.cookies.get('troplia_id') or ''

def _get_user_history(request: Request) -> list[dict]:
    sid = _get_sid_from_request(request)
    return HISTORY_BY_USER.get(sid, [])


# Aminoácidos válidos
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

# Massas monoisotópicas de RESÍDUOS (aa – H2O) em Da
AA_WEIGHTS = {
    'A': 71.037114, 'R': 156.101111, 'N': 114.042927, 'D': 115.026943,
    'C': 103.009185, 'E': 129.042593, 'Q': 128.058578, 'G': 57.021464,
    'H': 137.058912, 'I': 113.084064, 'L': 113.084064, 'K': 128.094963,
    'M': 131.040485, 'F': 147.068414, 'P': 97.052764, 'S': 87.032028,
    'T': 101.047678, 'V': 99.068414, 'W': 186.079313, 'Y': 163.063329
}

# Extinção molar (280 nm)
EXT_W = 5500  # Trp
EXT_Y = 1490  # Tyr

# ===== i18n (backend p/ exportação) =====
SUPPORTED_LANGS = ['pt','en']

def pick_lang_from_header(h: str | None) -> str:
    if not h:
        return 'en'
    parts = [p.split(';')[0].strip().lower() for p in h.split(',')]
    for p in parts:
        base = p.split('-')[0]
        if base in SUPPORTED_LANGS:
            return base
    return 'en'

# nomes localizados dos modos (para CSV/TXT)
CALC_NAME_I18N = {
    'pt': {'Com W':'Com W','Com YYY':'Com YYY','Sem W':'Sem W'},
    'en': {'Com W':'With W','Com YYY':'With YYY','Sem W':'No W'},}

OBS_I18N = {
    'pt': {'Com W':'Modo Com W (A280, W+Y) aplicado.',
           'Com YYY':'Modo Com YYY (A280, W+Y) aplicado.',
           'Sem W':'Modo Sem W (205/215/225) aplicado.'},
    'en': {'Com W':'Applied "With W" (A280, W+Y) mode.',
           'Com YYY':'Applied "With YYY" (A280, W+Y) mode.',
           'Sem W':'Applied "No W" (205/215/225) mode.'},}

# Para TXT, sempre com símbolo µM:
TXT_LABELS = {
    'pt': {'project':'Projeto','date':'Data','peptide':'Peptídeo','sequence':'Sequência','type':'Tipo',
           'result':'Concentração','obs':'Observação','unit':'µM'},
    'en': {'project':'Project','date':'Date','peptide':'Peptide','sequence':'Sequence','type':'Type',
           'result':'Concentration','obs':'Notes','unit':'µM'},}

CSV_HEADERS = {
    'pt': ["Projeto","Data (YYYY-MM-DD)","Peptídeo","Sequência","Tipo","A205","A215","A225","A280",
           "Diluição","#Y","#W","MM (Da)","µM","Obs"],
    'en': ["Project","Date (YYYY-MM-DD)","Peptide","Sequence","Type","A205","A215","A225","A280",
           "Dilution","#Y","#W","MM (Da)","µM","Notes"],
    'es': ["Proyecto","Fecha (YYYY-MM-DD)","Péptido","Secuencia","Tipo","A205","A215","A225","A280",
           "Dilución","#Y","#W","MM (Da)","µM","Notas"],
    'fr': ["Projet","Date (YYYY-MM-DD)","Peptide","Séquence","Type","A205","A215","A225","A280",
           "Dilution","#Y","#W","MM (Da)","µM","Remarque"],
    'de': ["Projekt","Datum (YYYY-MM-DD)","Peptid","Sequenz","Typ","A205","A215","A225","A280",
           "Verdünnung","#Y","#W","MM (Da)","µM","Hinweis"],
    'it': ["Progetto","Data (YYYY-MM-DD)","Peptide","Sequenza","Tipo","A205","A215","A225","A280",
           "Diluizione","#Y","#W","MM (Da)","µM","Nota"],
    'ja': ["プロジェクト","日付 (YYYY-MM-DD)","ペプチド","配列","タイプ","A205","A215","A225","A280",
           "希釈","#Y","#W","MM (Da)","µM","注記"],
    'zh': ["项目","日期 (YYYY-MM-DD)","肽","序列","类型","A205","A215","A225","A280",
           "稀释倍数","#Y","#W","MM (Da)","µM","备注"],
    'ru': ["Проект","Дата (YYYY-MM-DD)","Пептид","Последовательность","Тип","A205","A215","A225","A280",
           "Разведение","#Y","#W","MM (Да)","µM","Примечание"],
}

def clean_seq(s: str) -> str:
    if not s:
        return ""
    s = s.upper().replace(" ", "")
    return ''.join(ch for ch in s if ch in VALID_AA)

def mw_from_seq(seq: str) -> float:
    if not seq:
        return 0.0
    return sum(AA_WEIGHTS.get(ch, 0.0) for ch in seq) + 18.01056

def autodetect_type(seq: str) -> str:
    y = seq.count('Y') if seq else 0
    w = seq.count('W') if seq else 0
    if w >= 1:
        return 'Com W'
    if w == 0 and y >= 3:
        return 'Com YYY'
    return 'Sem W'

def autodetect_type_counts(qw: int, qy: int) -> str:
    if qw >= 1:
        return 'Com W'
    if qw == 0 and qy >= 3:
        return 'Com YYY'
    return 'Sem W'

def canonical_calc_type(raw: str | None, seq: str) -> str:
    s = (raw or '').strip().lower()
    if s in {'com w','w','a280 (w)','a280 w','with w'}:
        return 'Com W'
    if s in {'com yyy','yyy','a280 (yyy)','a280 yyy','with yyy'}:
        return 'Com YYY'
    if s in {'sem w','205','a205','205/215/225','no w','without w'}:
        return 'Sem W'
    return autodetect_type(seq)

def slugify(text: str) -> str:
    text = text.strip().replace(" ", "_")
    text = re.sub(r"[^A-Za-z0-9_\-\.]+", "_", text)
    return text or "Troplia"

def make_download_headers(project_name: str, ext: str) -> dict:
    now_stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    base = f"{slugify(project_name)}_{now_stamp}.{ext}"
    quoted = urllib.parse.quote(base)
    return {"Content-Disposition": f'attachment; filename="{base}"; filename*=UTF-8\'\'{quoted}'}

# ===== HTML =====
INDEX_HTML = r"""
<!doctype html>
<html lang="pt">
<head>

 <!-- Favicon -->
  <link rel="icon" type="image/png" href="/static/favicon.png">
  <link rel="shortcut icon" type="image/png" href="/static/favicon.png">
   
  <!-- Google AdSense -->
  <script async src="https://pagead2.googlesyndication.com/pagead/js/adsbygoogle.js?client=ca-pub-9055264960112923"
     crossorigin="anonymous"></script>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>TROPLIA</title>

  <script src="https://cdn.tailwindcss.com"></script>
  <link href="https://fonts.googleapis.com/css2?family=Space+Grotesk:wght@700;800&family=Poppins:wght@400;600;800&display=swap" rel="stylesheet">

  <style>
    :root{
      --brand:#83E509;
      --page:#e2f7f5;
      --ink:#0c2a27;
      --muted:#3a5b56;
      --card:#ffffff;
      --border:#cfe9e5;
    }
    html,body{ background-color:var(--page); color:var(--ink); }
    body{
      font-family:Poppins, system-ui, -apple-system, Segoe UI, Roboto, Ubuntu, Cantarell, Noto Sans, Arial;
      text-align: justify;
    }
    .title-brand{
      font-family:"Space Grotesk", Poppins, system-ui, -apple-system, Segoe UI, Roboto, Ubuntu, Cantarell, Noto Sans, Arial, sans-serif;
      font-weight:800; color:#6cd4c4; text-transform:uppercase; letter-spacing:.08em;
    }
    .card{ background-color:var(--card); border:1px solid var(--border); }
    .input,.select,.textarea{ background-color:#fff; color:var(--ink); border:1px solid var(--border); }
    .input::placeholder{ color:#708c88; }
    
    .link{ color:#0a7a5d; } .link:hover{ text-decoration:underline; }
    .warn{ color:#a05a00; }

    .logo-wrap{ background-color:var(--page); border-radius:1.25rem; padding:.5rem; }
    .logo-img{ display:block; height:11.5rem; width:auto; background-color:var(--page); }

    .troplia-title-size{ font-size:2.75rem; }
    @media (min-width:768px){ .troplia-title-size{ font-size:3.75rem; } }

    nav.tabs { position:relative; z-index:5; }
    nav.tabs button { cursor:pointer; }
    .clickable a { pointer-events:auto; }
  </style>

<style>
  /* Destaque visível para a aba ativa do topo */
  
</style>


<style>
/* Tab active: single thick underline (no double underline) */
.tab-active{
  font-weight: 700 !important;
  border-bottom: 3px solid var(--ink) !important;
  text-decoration: none !important;
}
</style>

</head>
<  <link rel="stylesheet" href="/static/helipep.css">
  <script src="/static/helipep.js" defer></script>
<body>

  <div class="mx-auto max-w-7xl px-3 sm:px-6 py-6">

    <header class="max-w-6xl mx-auto flex items-center justify-between gap-4">
      <div class="flex items-center gap-4">
        <div class="logo-wrap"><img src="/static/logo.png" alt="Logo Troplia" class="logo-img" id="logo-img" /></div>
        <div><h1 class="troplia-title-size title-brand">TROPLIA</h1></div>
      </div>
      <div class="flex items-center gap-2 text-sm">
        <span aria-hidden="true">🌐</span>
        <label for="lang" class="sr-only">Idioma</label>
        <select id="lang" class="select rounded px-2 py-1">
          <option value="pt">🇧🇷 Português</option>
          <option value="en">🇺🇸 English</option>
        </select>
      </div>
    </header>

    <nav class="tabs max-w-6xl mx-auto mt-6 border-b border-[color:var(--border)] flex flex-wrap gap-6 text-sm sm:text-base">
      <button type="button" id="tab-home"    data-tab="home"    class="py-2 tab-active" data-i18n="tab_home">Início</button>
      <button type="button" id="tab-news"    data-tab="news"    class="py-2" data-i18n="tab_news">Novidades</button>
      <button type="button" id="tab-about"   data-tab="about"   class="py-2" data-i18n="tab_about">Sobre Nós</button>
      <button type="button" id="tab-contact" data-tab="contact" class="py-2" data-i18n="tab_contact">Contato</button>
    </nav>

    <main class="mt-6 grid grid-cols-1 lg:grid-cols-[1fr_minmax(240px,280px)] gap-6">
      <div class="max-w-6xl mx-auto w-full">

        <!-- INÍCIO -->
        <section id="pane-home" class="block">
          <div class="card rounded-2xl shadow p-6">
            <h2 class="text-2xl font-bold mb-2" data-i18n="home_title">Bem-vindo à Troplia</h2>
            <p class="mb-4 text-[color:var(--muted)]" data-i18n="home_intro">
              Esta é a plataforma da Troplia para acesso a diversas ferramentas utilizadas em laboratório. O site será atualizado com novas ferramentas; no momento, conta com a <b>QuantyPep</b> para cálculo e registro de concentração de peptídeos por absorbância — simples, rápida e precisa.
            </p>
          </div>

          <div class="card rounded-2xl shadow p-6 mt-6 clickable">
            <h3 class="text-xl font-semibold mb-3" data-i18n="tools_title">Ferramentas</h3>
            <div class="grid sm:grid-cols-2 lg:grid-cols-3 gap-4">
              <a href="#" id="go-calc" class="block rounded-xl border border-[color:var(--border)] p-4 hover:shadow transition bg-white">
                <div class="font-semibold mb-1" data-i18n="tool_calc_title">QuantyPep - Calculadora de Concentração de Peptídeos (Absorbância)</div>
                <p class="text-sm text-[color:var(--muted)]" data-i18n="tool_calc_desc" data-i18n="tool_calc_desc">
                  Calcule concentrações por A280 (Com W ou Com YYY) ou por A205/215/225 (Sem W). Resultados sempre em micromolar (µM).
                </p>
              </a>
<!-- Card: HeliPep (iguais em estrutura/estilo ao QuantyPep) -->
<a id="go-helipep" href="#helipep" class="tool-card">
  <h3><strong>HeliPep — Helical Wheel Plot para Peptídeos</strong></h3>
  <p data-i18n="tool_helipep_desc">Gere helical wheel plots (preview + download), lote via .txt e personalização de cores/formas.</p>
</a>
                </p>
              </a>
            </div>
          </div>
        </section>

        <!-- CALCULADORA -->
        <section id="pane-calc" class="hidden" aria-hidden="true">
          <div class="card rounded-2xl shadow p-6 mb-6">
            <h2 class="text-2xl font-semibold" data-i18n="calc_title">QuantyPep — Calculadora de Concentração de Peptídeos (Absorbância)</h2>
          </div>

          <div class="mt-1 grid md:grid-cols-2 gap-6">
            <div class="card rounded-2xl shadow p-6">
              <h3 class="text-lg font-semibold mb-3" data-i18n="proj_info">Informações do Projeto</h3>
              <div class="space-y-3">
                <div><label class="block text-sm mb-1" data-i18n="lbl_project">Nome do Projeto</label><input id="project" class="input w-full rounded px-3 py-2" placeholder="Meu experimento"/></div>
                <div><label class="block text-sm mb-1" data-i18n="lbl_date">Data</label><input id="date" type="date" class="input w-full rounded px-3 py-2"/></div>
                <div><label class="block text-sm mb-1" data-i18n="lbl_peptide_name">Nome do Peptídeo</label><input id="peptide_name" class="input w-full rounded px-3 py-2" placeholder="Ex.: Pep-1"/></div>
                <div>
                  <label class="block text-sm mb-1"><span data-i18n="lbl_sequence">Sequência do Peptídeo</span> <span class="opacity-70" data-i18n="opt_field">(opcional)</span></label>
                  <input id="sequence" class="input w-full rounded px-3 py-2" placeholder="Ex.: ACDEFGHIKLMNPQRSTVWY"/>
                  <p id="seqWarn" class="text-xs warn mt-1 hidden"></p>
                </div>
              </div>
            </div>

            <div class="card rounded-2xl shadow p-6">
              <div class="flex items-center justify-between mb-1">
                <h3 class="text-lg font-semibold" data-i18n="calc_type">Tipo de Cálculo</h3>
                <small id="calcTip" class="opacity-80"></small>
              </div>
              <div class="flex gap-2 mb-3">
                <select id="calc_type" class="select rounded px-3 py-2">
                  <option value="Auto" data-i18n-opt="calc_opt_auto">Auto</option>
                  <option value="Com W" data-i18n-opt="calc_opt_w">Com W (A280)</option>
                  <option value="Com YYY" data-i18n-opt="calc_opt_yyy">Com YYY (A280)</option>
                  <option value="Sem W" data-i18n-opt="calc_opt_now">Sem W (205/215/225)</option>
                </select>
                <button type="button" id="btnSuggest" class="border border-[color:var(--border)] rounded px-3 py-2" title="Sugerir pelo peptídeo">🔎</button>
              </div>

              <div id="fields" class="grid grid-cols-1 sm:grid-cols-2 gap-3"></div>

              <div class="flex flex-wrap gap-3 mt-4">
                <button type="button" id="btnCompute" class="bg-[color:var(--brand)] text-black rounded px-4 py-2 font-semibold" data-i18n="btn_compute">Computar</button>
                <button type="button" id="btnClear" class="bg-white border border-[color:var(--border)] rounded px-4 py-2" data-i18n="btn_clear">Limpar</button>
                <a id="export_txt_main" href="/api/export/txt" class="bg-white border border-[color:var(--border)] rounded px-4 py-2" data-i18n="btn_save_txt">Salvar TXT</a>
                <a id="export_csv_main" href="/api/export/csv" class="bg-white border border-[color:var(--border)] rounded px-4 py-2" data-i18n="btn_export_csv">Exportar CSV</a>
              </div>
            </div>
          </div>

          <div class="card rounded-2xl shadow p-6 mt-6">
            <h3 class="text-lg font-semibold mb-2" data-i18n="result">Resultado</h3>
            <div id="results" class="text-sm sm:text-base whitespace-pre-line"></div>
          </div>
        
    <!-- === Contêiner interno QuantyPep para receber as abas realocadas === -->
    <div id="qp-internal-after" class="mt-6">
      <div class="flex gap-3 mb-3">
        <button type="button" id="qp-open-hist"
                class="px-4 py-2 rounded-xl font-semibold shadow hover:opacity-90 transition bg-[#00b894] text-white"
                data-i18n="tab_hist">📜 Histórico</button>
        <button type="button" id="qp-open-help"
                class="px-4 py-2 rounded-xl font-semibold shadow hover:opacity-90 transition bg-[#0984e3] text-white"
                data-i18n="tab_help">❓ Instruções</button>
      </div>
      <!-- As seções pane-hist e pane-help serão movidas para cá quando a QuantyPep abrir -->
    </div>
    </section>
        <!-- HISTÓRICO -->
        <section id="pane-hist" class="hidden" aria-hidden="true">
          <div class="card rounded-2xl shadow p-6">
            <div class="flex items-center justify-between gap-3 mb-4">
              <h2 class="text-lg font-semibold" data-i18n="hist_title">Histórico</h2>
              <div class="flex gap-2">
                <a id="export_txt_hist" href="/api/export/txt" class="bg-white border border-[color:var(--border)] rounded px-3 py-2" data-i18n="btn_save_txt">Salvar TXT</a>
                <a id="export_csv_hist" href="/api/export/csv" class="bg-white border border-[color:var(--border)] rounded px-3 py-2" data-i18n="btn_export_csv">Exportar CSV</a>
                <button type="button" id="btnClearHistory" class="bg-white border border-[color:var(--border)] rounded px-3 py-2" data-i18n="btn_clear_history">Limpar Histórico</button>
              </div>
            </div>
            <div class="overflow-x-auto">
              <table class="table-auto w-full text-sm border border-[color:var(--border)]">
                <thead class="bg-[#f3fbf9]">
                  <tr>
                    <th class="p-2" data-i18n="th_project">Projeto</th>
                    <th class="p-2" data-i18n="th_date">Data (YYYY-MM-DD)</th>
                    <th class="p-2" data-i18n="th_peptide">Peptídeo</th>
                    <th class="p-2" data-i18n="th_type">Tipo</th>
                    <th class="p-2">A205</th><th class="p-2">A215</th><th class="p-2">A225</th><th class="p-2">A280</th>
                    <th class="p-2" data-i18n="th_dil">Diluição</th>
                    <th class="p-2">#Y</th>
                    <th class="p-2">#W</th>
                    <th class="p-2">MM (Da)</th>
                    <th class="p-2">µM</th>
                    <th class="p-2" data-i18n="th_obs">Obs</th>
                  </tr>
                </thead>
                <tbody id="history" class="align-top"></tbody>
              </table>
            </div>
          </div>
        </section>

        <!-- INSTRUÇÕES -->
        <section id="pane-help" class="hidden" aria-hidden="true">
          <div class="card rounded-2xl shadow p-6">
            <h2 class="text-lg font-semibold mb-3" data-i18n="help_title">Instruções</h2>
            <div class="text-sm sm:text-base space-y-3">
              <p data-i18n="help_step"></p>
              <p data-i18n="help_auto"></p>
              <p data-i18n="help_w"></p>
              <p data-i18n="help_yyy"></p>
              <p data-i18n="help_now"></p>
              <p data-i18n="help_res"></p>
              <p data-i18n="help_mm"></p>
              <p data-i18n="help_eps"></p>
            </div>
          </div>
        </section>
                <!-- NOVIDADES -->
        <section id="pane-news" class="hidden" aria-hidden="true">
          <div class="card rounded-2xl shadow p-6">
            <h2 class="text-2xl font-bold mb-3" data-i18n="news_title">Novidades</h2>

            <h3 class="font-semibold mb-2" data-i18n="news_helipep_heading">Lançamento da Ferramenta HeliPep</h3>
            <p class="text-[color:var(--muted)]" data-i18n="news_p1">
              Lançamos a <b>HeliPep — Helical Wheel Plot para Peptídeos</b>, a mais nova ferramenta integrada à plataforma Troplia.
            </p>

            <p class="mt-3 text-[color:var(--muted)]" data-i18n="news_p2">
              O HeliPep foi desenvolvido para atender pesquisadores que buscam visualizar e analisar propriedades estruturais de
              peptídeos por meio de <i>helical wheel plots</i> (roda helicoidal). Com ele, você pode:
            </p>

            <ul class="list-disc pl-6 mt-2 text-[color:var(--muted)]">
              <li data-i18n="news_b1"><b>Gerar visualizações instantâneas</b> de peptídeos a partir de sequências individuais ou múltiplas em lote (.txt).</li>
              <li data-i18n="news_b2"><b>Exportar imagens em alta qualidade (PNG)</b>, prontas para artigos, relatórios e apresentações.</li>
              <li data-i18n="news_b3"><b>Personalizar cores e formas</b> de resíduos hidrofóbicos, básicos, ácidos e polares.</li>
              <li data-i18n="news_b4"><b>Explorar opções avançadas</b>, como guias circulares e momento hidrofóbico.</li>
            </ul>

            <p class="mt-3 text-[color:var(--muted)]" data-i18n="news_p3">
              Assim como no QuantyPep, a interface foi projetada para ser <b>simples, intuitiva e responsiva</b>.
              O HeliPep chega para complementar a Troplia e oferecer uma solução poderosa para análise visual de peptídeos —
              alinhada às exigências de publicações de alto impacto.
            </p>

            <p class="mt-3 font-semibold" data-i18n="news_cta">Experimente agora a ferramenta no menu principal da plataforma! 🚀</p>
          </div>
</section>

         <!-- HELIPEP -->
<section id="pane-helipep" class="hidden" aria-hidden="true">
  <div class="card rounded-2xl shadow p-6 mb-6" style="margin:16px auto; max-width:1100px;">
    <h2 class="text-2xl font-semibold" data-i18n="helipep_title">HeliPep — Helical Wheel Plot para peptídeos</h2>
    <p class="muted" data-i18n="helipep_subtitle" data-i18n="tool_helipep_desc">Gere helical wheel plots a partir de uma sequência única ou arquivo .txt com várias sequências.</p>

    <div class="grid-2" style="display:grid; grid-template-columns:1fr 1fr; gap:16px;">
      <div>
        <h3 class="section-title"><strong data-i18n="helipep_single_title">Sequência única</strong></h3>
        <label for="helipep-seq" data-i18n="helipep_seq_label">Sequência (apenas letras):</label>
        <textarea id="helipep-seq" rows="4" placeholder="ACDEFGHIKLMNPQRSTVWY"></textarea>

        <label for="helipep-name" data-i18n="helipep_name_label">Nome do peptídeo (opcional; aparece no centro):</label>
        <input id="helipep-name" type="text" placeholder="Pep-1">

        <div class="row" style="display:flex; gap:10px; align-items:center; flex-wrap:wrap; margin:8px 0;">
          <button id="helipep-preview" type="button" class="btn" data-i18n="helipep_preview_btn">Pré-visualizar</button>
          <a id="helipep-preview-dl" class="btn" href="#" download="helipep_preview.png" data-i18n="helipep_download_btn">Baixar PNG</a>
        </div>

        <div class="img-wrap" style="margin-top:10px;">
          <img id="helipep-img" alt="Helical wheel preview" style="max-width:420px;">
        </div>
      </div>

      <div>
        <h3 class="section-title"><strong data-i18n="helipep_bulk_title">Várias sequências (.txt)</strong></h3>
        <p class="muted" style="margin-top:4px;" data-i18n="helipep_bulk_hint">Exemplo de arquivo .txt (uma por linha):</p>
        <div class="file-pt-wrapper">
          <label for="helipep-file" class="btn" data-i18n="helipep_select_file_btn">Selecionar arquivo</label>
          <input id="helipep-file" type="file" accept=".txt" class="file-pt-input" />
          <span id="helipep-file-name" class="file-pt-name" data-i18n="helipep_file_none">Nenhum arquivo selecionado</span>
        </div>
        <div class="row" style="display:flex; gap:10px; align-items:center; flex-wrap:wrap; margin:8px 0;">
          <button id="helipep-bulk" class="btn" data-i18n="helipep_bulk_btn">Baixar resultado</button> 
        </div>
      </div>
    </div>

    <h3 class="section-title" style="margin-top:18px;"><strong data-i18n="helipep_custom_title">Personalização</strong></h3>
    <div class="grid-3" style="display:grid; grid-template-columns:repeat(3, minmax(180px, 1fr)); gap:12px;">
      <div><label data-i18n="helipep_color_hydrophobic">Cor hidrofóbicos</label><input type="color" id="color_hydrophobic" value="#008000"></div>
      <div><label data-i18n="helipep_color_basic">Cor básicos</label><input type="color" id="color_basic"       value="#FF0000"></div>
      <div><label data-i18n="helipep_color_acidic">Cor ácidos</label><input type="color" id="color_acidic"      value="#0000FF"></div>
      <div><label data-i18n="helipep_color_polar">Cor polares</label><input type="color" id="color_polar"      value="#ADD8E6"></div>

      <div>
        <label data-i18n="helipep_shape_hydrophobic">Forma hidrofóbicos</label>
        <select id="shape_hydrophobic">
          <option value="">●</option><option value="s">■</option><option value="^">▲</option><option value="D">◆</option>
        </select>
      </div>
      <div>
        <label data-i18n="helipep_shape_basic">Forma básicos</label>
        <select id="shape_basic">
          <option value="">●</option><option value="s">■</option><option value="^">▲</option><option value="D">◆</option>
        </select>
      </div>
      <div>
        <label data-i18n="helipep_shape_acidic">Forma ácidos</label>
        <select id="shape_acidic">
          <option value="">●</option><option value="s">■</option><option value="^">▲</option><option value="D">◆</option>
        </select>
      </div>
      <div>
        <label data-i18n="helipep_shape_polar">Forma polares</label>
        <select id="shape_polar">
          <option value="">●</option><option value="s">■</option><option value="^">▲</option><option value="D">◆</option>
        </select>
      </div>
    </div>

    <details style="margin-top:10px;">
      <summary class="section-title"><strong data-i18n="helipep_advanced_title">Opções avançadas</strong></summary>
      <div class="grid-3" style="display:grid; grid-template-columns:repeat(3, minmax(180px, 1fr)); gap:12px; margin-top:8px;">
        <div><label data-i18n="helipep_opt_angle">Δ°/resíduo</label><input id="opt_angle" type="number" step="0.1" value="100"></div>
        <div><label data-i18n="helipep_opt_start">Ângulo inicial</label><input id="opt_start" type="number" step="0.1" value="0"></div>
        <div><label data-i18n="helipep_opt_clockwise">Clockwise</label><input id="opt_clockwise" type="checkbox" checked></div>
        <div><label data-i18n="helipep_opt_layer_n">Camada a cada N resíduos</label><input id="opt_layer_n" type="number" value="18"></div>
        <div><label data-i18n="helipep_opt_links">Ligar resíduos</label><input id="opt_links" type="checkbox" checked></div>
        <div><label data-i18n="helipep_opt_guide">Círculo guia</label><input id="opt_guide" type="checkbox"></div>
        <div><label data-i18n="helipep_opt_nums">Números</label><input id="opt_nums" type="checkbox" checked></div>
        <div><label data-i18n="helipep_opt_letters">Letras</label><input id="opt_letters" type="checkbox" checked></div>
        <div><label data-i18n="helipep_opt_moment">Momento hidrofóbico</label><input id="opt_moment" type="checkbox"></div>
      </div>
    </details>
    <div class="refs" style="margin-top:10px;">
      <h4 class="section-title"><strong data-i18n="helipep_refs_title">Referências</strong></h4>
      <small>
        Schiffer, M., & Edmundson, A. B. (1967). Use of helical wheels to represent the structures of proteins and to identify periodicities of physicochemical properties. <em>Biophysical Journal</em>.<br>
        Eisenberg, D., Weiss, R. M., & Terwilliger, T. C. (1982). The hydrophobic moment: a measure of the amphipilicity of a helix. <em>Nature</em>.
      </small>
    </div>
  </div>
</section>
        <!-- SOBRE NÓS -->
        <section id="pane-about" class="hidden" aria-hidden="true">
          <div class="card rounded-2xl shadow p-6">
            <h2 class="text-2xl font-bold mb-3" data-i18n="about_title">Sobre Nós</h2>
            <h3 class="font-semibold mb-2" data-i18n="about_sub">Sobre a Troplia</h3>
            <p class="text-[color:var(--muted)]" data-i18n="about_body"></p>
            <h3 class="font-semibold mt-4 mb-2" data-i18n="about_maint_title">Manutenção do Site</h3>
            <p class="text-[color:var(--muted)]" data-i18n="about_maint_body"></p>
          </div>
        </section>

        <!-- CONTATO -->
        <section id="pane-contact" class="hidden" aria-hidden="true">
          <div class="card rounded-2xl shadow p-6">
            <h2 class="text-2xl font-bold mb-3" data-i18n="contact_title">Contato</h2>
            <p class="mb-3" data-i18n="contact_intro"></p>
            <p data-i18n="contact_list"></p>
          </div>
        </section>

      </div>
      <aside class="w-full lg:w-auto"><div class="sticky top-4 space-y-4"></div></aside>
    </main>
  </div>

<script>
/* ===== i18n ===== */
let CURRENT_LANG = 'pt';

const I18N = {
  pt:{
    /* Tabs / Navegação */
    tab_home:'Início', tab_hist:'Histórico', tab_help:'Instruções', tab_news:'Novidades', tab_about:'Sobre Nós', tab_contact:'Contato',

    /* Home */
    home_title:'Bem-vindo à Troplia',
    home_intro:'A Troplia se propõe a desenvolver soluções inteligentes que unem ciência e tecnologia para impactar positivamente a sociedade e o setor produtivo. Nosso site está em constante construção para auxiliar em diferentes frentes de inovação, oferecendo futuramente novas funcionalidades e integrações. No momento, conta com ferramentas de análise laboratorial voltadas para apoiar pesquisadores de diversas áreas no estudo de peptídeos, fornecendo cálculos e visualizações que tornam a pesquisa mais simples, rápida e precisa.',
    tools_title:'Ferramentas',

    /* Cards de ferramentas (existente + novo) */
    tool_calc_title:'QuantyPep - Calculadora de Concentração de Peptídeos (Absorbância)',
    tool_calc_desc:'Calcule concentrações por A280 (Com W ou Com YYY) ou por A205/215/225 (Sem W). Resultados sempre em micromolar (µM).',
    tool_helipep_title:'HeliPep — Helical Wheel Plot para Peptídeos',
    tool_helipep_desc:'Gere helical wheel plots instantaneamente (preview + download), suporte a .txt (várias sequências) e personalização de cores/formatos.',

    /* QuantyPep (Calculadora) */
    calc_title:'QuantyPep — Calculadora de Concentração de Peptídeos (Absorbância)',
    proj_info:'Informações do Projeto', lbl_project:'Nome do Projeto', lbl_date:'Data', lbl_peptide_name:'Nome do Peptídeo', lbl_sequence:'Sequência do Peptídeo', opt_field:'(opcional)',
    ph_project:'Meu experimento', ph_peptide_name:'Ex.: Pep-1', ph_sequence:'Ex.: ACDEFGHIKLMNPQRSTVWY',
    calc_type:'Tipo de Cálculo', btn_compute:'Computar', btn_clear:'Limpar', btn_export_csv:'Exportar CSV', btn_save_txt:'Salvar TXT', result:'Resultado',
    calc_opt_auto:'Auto', calc_opt_w:'Com W (A280)', calc_opt_yyy:'Com YYY (A280)', calc_opt_now:'Sem W (205/215/225)',
    calc_name_w:'Com W', calc_name_yyy:'Com YYY', calc_name_now:'Sem W',
    tip_w:'A280 (W+Y): ε_total = (#W×5500 + #Y×1490).', tip_yyy:'A280 (W+Y): ε_total = (#W×5500 + #Y×1490).', tip_now:'Usa 205/215/225 nm.',
    obs_w:'Modo Com W (A280, W+Y) aplicado.', obs_yyy:'Modo Com YYY (A280, W+Y) aplicado.', obs_now:'Modo Sem W (205/215/225) aplicado.',
    hist_title:'Histórico', th_project:'Projeto', th_date:'Data (YYYY-MM-DD)', th_peptide:'Peptídeo', th_type:'Tipo', th_dil:'Diluição', th_obs:'Obs',
    btn_clear_history:'Limpar Histórico',
    help_title:'Instruções',
    help_step:'<b>Passo a passo:</b> Informe <i>Nome do Projeto</i>, <i>Data</i> (YYYY-MM-DD), <i>Nome do Peptídeo</i> e, opcionalmente, a <i>Sequência</i>.',
    help_auto:'<b>Auto:</b> se a sequência tiver <b>W ≥ 1</b> → <i>Com W</i>; se <b>W = 0</b> e <b>Y ≥ 3</b> → <i>Com YYY</i>; caso contrário → <i>Sem W</i>.<br><b>Disclaimer:</b> <i>Com W</i> e <i>Com YYY</i> usam a <u>mesma calculadora A280</u>. Para A280, usamos <b>ε_total = (#W×5500 + #Y×1490)</b>. Se #W ou #Y for 0, esse termo não contribui. Peptídeos com <b>1–2 Y</b> e <b>nenhum W</b> devem usar <i>Sem W</i>.',
    help_w:'<b>Com W / Com YYY (A280):</b> <code>ε_total = #W×5500 + #Y×1490</code> (M⁻¹·cm⁻¹).<br><code>µM = (A280 / ε_total) × Dil × 1.000.000</code>.',
    help_yyy:'<b>Observação:</b> Peptídeos com ≥3 tirosinas e sem W usam a mesma fórmula A280; se houver W e Y, ambos contribuem via ε_total.',
    help_now:'<b>Sem W (205/215/225):</b> <code>X=(A215−A225)×144</code>; <code>Y=A205×31</code>; <code>mg/mL=((X+Y)/2)×Dil</code>; <code>µM=(mg/mL ÷ MM)×1000</code>.',
    help_res:'<b>Resultado:</b> sempre em <b>µM</b>.',
    help_mm:'<b>Disclaimer (MM):</b> A massa molar é calculada automaticamente pela sequência, mas você pode informar qualquer valor.',
    help_eps:'<b>Coeficientes (280 nm):</b> W = <b>5500</b>; Y = <b>1490</b> M⁻¹·cm⁻¹.',

    /* News */
    news_title:'Novidades',
    news_body:'Estamos trabalhando em novas ferramentas, mais idiomas e integrações. Em breve publicaremos atualizações aqui.',
    news_helipep_title:'Novo: HeliPep',
    news_helipep_body:'Apresentamos o HeliPep — gere helical wheel plots com rapidez, personalize cores/formas e exporte um .zip para múltiplas sequências.',

    /* About & Contact */
    about_title:'Sobre Nós', about_sub:'Sobre a Troplia',
    about_body:'<b>O que é a Troplia</b><br>A Troplia visa desenvolver soluções inteligentes para a sociedade e melhorar processos de produção, por exemplo na área alimentar, fornecendo ferramentas para uma melhor tomada de decisões por produtores de todos os portes, bem como desenvolvendo produtos capazes de controlar pragas agrícolas, que sejam amigáveis ao meio ambiente e mais seguros ao consumidor.<br><br><b>Nossa Missão</b><br>Nossa missão é transformar ciência e tecnologia em soluções práticas que impactem positivamente a vida das pessoas e o futuro do planeta. Buscamos unir biotecnologia, inteligência artificial e inovação sustentável para desenvolver produtos e serviços que melhorem processos produtivos, promovam maior segurança alimentar e ambiental e facilitem a tomada de decisões estratégicas em setores essenciais como a agropecuária e a indústria. Ao mesmo tempo, temos como objetivo democratizar o acesso a ferramentas analíticas e educacionais, capacitando produtores, pesquisadores e estudantes a explorarem novas possibilidades e contribuírem para os Objetivos de Desenvolvimento Sustentável da Agenda 2030 da ONU. Acreditamos que ciência e sociedade caminham juntas, e trabalhamos para ser um elo forte entre inovação, responsabilidade social e preservação ambiental.<br><br><b>Áreas de atuação</b><br>Biotecnologia, inteligência artificial, agropecuária e ensino.<br><br><b>Serviços e produtos</b><ul><li>Soluções biotecnológicas para o mercado agropecuário;</li><li>Desenvolvimento de produtos agropecuários mais seguros e amigáveis ao meio ambiente, em consonância com a Agenda 2030 da ONU;</li><li>Desenvolvimento de ferramentas sob demanda para ambiente acadêmico e industrial utilizando redes neurais;</li><li>Análises cromatográficas e de espectrometria de massa;</li><li>Consultoria em análise de dados e desenvolvimento de projetos;</li><li>Desenvolvimento de cursos de capacitação em prospecção e abordagens analíticas.</li></ul><br><b>Equipe</b><ul><li>Dr. Marcelo Henrique Soller Ramada;</li><li>Dr. Stephan Machado Dohms;</li><li>M.Sc. Lorena Ferreira Peixoto;</li><li>M.Sc. Daniel Gusmão de Morais.</li></ul>',
    about_maint_title:'Manutenção do Site',
    about_maint_body:'Equipe responsável:<br><br>Gabriel Iudy Yamaguchi Rocha<br>Marcelo Henrique Soller Ramada<br>Rivaldo Varejão Pasqual Saraiva',
    contact_title:'Contato', contact_intro:'Entre em contato para informações ou parcerias:',
    contact_list:'Gabriel Iudy Yamaguchi Rocha: <a class="link" href="mailto:yamaguchiiudy@gmail.com">yamaguchiiudy@gmail.com</a><br>Marcelo Henrique Soller Ramada: <a class="link" href="mailto:marceloramada@gmail.com">marceloramada@gmail.com</a><br>Rivaldo Varejão Pasqual Saraiva: <a class="link" href="mailto:rivaldo.vps8@gmail.com">rivaldo.vps8@gmail.com</a>',

    /* HeliPep — NOVO */
    helipep_title:'HeliPep — Helical Wheel Plot para peptídeos',
    helipep_subtitle:'Gere helical wheel plots a partir de uma sequência única ou de um arquivo .txt com várias sequências.',
    helipep_single_title:'Sequência única',
    helipep_seq_label:'Sequência (apenas letras):',
    helipep_name_label:'Nome do peptídeo (opcional; aparece no centro):',
    helipep_preview_btn:'Pré-visualizar',
    helipep_download_btn:'Baixar PNG',
    helipep_bulk_title:'Várias sequências (.txt)',
    helipep_bulk_hint:'Exemplo de arquivo .txt (uma por linha):',
    helipep_bulk_btn:'Baixar resultado',
    helipep_custom_title:'Personalização',
    helipep_color_hydrophobic:'Cor hidrofóbicos',
    helipep_color_basic:'Cor básicos',
    helipep_color_acidic:'Cor ácidos',
    helipep_color_polar:'Cor polares',
    helipep_shape_hydrophobic:'Forma hidrofóbicos',
    helipep_shape_basic:'Forma básicos',
    helipep_shape_acidic:'Forma ácidos',
    helipep_shape_polar:'Forma polares',
    helipep_advanced_title:'Opções avançadas',
    helipep_refs_title:'Referências'
  , label_conc:'Concentração', unit_micromolar:'µM'},

  en:{
    tab_home:'Home', tab_hist:'History', tab_help:'Instructions', tab_news:'News', tab_about:'About Us', tab_contact:'Contact',
    home_title:'Welcome to Troplia',
    home_intro:'Troplia is committed to developing intelligent solutions that bring together science and technology to positively impact society and the productive sector. Our website is constantly evolving to support different areas of innovation, and will soon offer new features and integrations. At the moment, it provides laboratory analysis tools designed to support researchers from various fields in the study of peptides, offering calculations and visualizations that make research simpler, faster, and more accurate.',
    tools_title:'Tools',
    tool_calc_title:'QuantyPep — Peptide Concentration Calculator (Absorbance)',
    tool_calc_desc:'Calculate concentrations using A280 (With W or With YYY) or A205/215/225 (Without W). Results are always given in micromolar (µM).',
    tool_helipep_title:'HeliPep — Helical Wheel Plot for Peptides',
    tool_helipep_desc:'Generate helical wheel plots (preview + download), batch processing via .txt, and color/shape customization.',
    calc_title:'QuantyPep — Peptide Concentration Calculator (Absorbance)',
    proj_info:'Project Info', lbl_project:'Project Name', lbl_date:'Date', lbl_peptide_name:'Peptide Name', lbl_sequence:'Peptide Sequence', opt_field:'(optional)',
    ph_project:'My experiment', ph_peptide_name:'e.g., Pep-1', ph_sequence:'e.g., ACDEFGHIKLMNPQRSTVWY',
    calc_type:'Calculation Type', btn_compute:'Compute', btn_clear:'Clear', btn_export_csv:'Export CSV', btn_save_txt:'Save TXT', result:'Result',
    calc_opt_auto:'Auto', calc_opt_w:'With W (A280)', calc_opt_yyy:'With YYY (A280)', calc_opt_now:'No W (205/215/225)',
    calc_name_w:'With W', calc_name_yyy:'With YYY', calc_name_now:'No W',
    tip_w:'A280 (W+Y): ε_total = (#W×5500 + #Y×1490).', tip_yyy:'A280 (W+Y): ε_total = (#W×5500 + #Y×1490).', tip_now:'Uses 205/215/225 nm.',
    obs_w:'Applied "With W" (A280, W+Y) mode.', obs_yyy:'Applied "With YYY" (A280, W+Y) mode.', obs_now:'Applied "No W" (205/215/225) mode.',
    hist_title:'History', th_project:'Project', th_date:'Date (YYYY-MM-DD)', th_peptide:'Peptide', th_type:'Type', th_dil:'Dilution', th_obs:'Notes',
    btn_clear_history:'Clear History',
    help_title:'Instructions',
    help_step:'<b>Steps:</b> Provide <i>Project Name</i>, <i>Date</i> (YYYY-MM-DD), <i>Peptide Name</i> and optionally the <i>Sequence</i>.',
    help_auto:'<b>Auto:</b> if sequence has <b>W ≥ 1</b> → <i>With W</i>; if <b>W = 0</b> and <b>Y ≥ 3</b> → <i>With YYY</i>; else → <i>No W</i>.<br><b>Disclaimer:</b> <i>With W</i> and <i>With YYY</i> use the <u>same A280 calculator</u>. For A280 we use <b>ε_total = (#W×5500 + #Y×1490)</b>. If #W or #Y is 0, that term doesn’t contribute. Peptides with <b>1–2 Y</b> and <b>no W</b> should use <i>No W</i>.',
    help_w:'<b>With W / With YYY (A280):</b> <code>ε_total = #W×5500 + #Y×1490</code> (M⁻¹·cm⁻¹).<br><code>µM = (A280 / ε_total) × Dil × 1,000,000</code>.',
    help_yyy:'<b>Note:</b> Peptides with ≥3 tyrosines and no W use the same A280 formula; with both W and Y, both contribute via ε_total.',
    help_now:'<b>No W (205/215/225):</b> <code>X=(A215−A225)×144</code>; <code>Y=A205×31</code>; <code>mg/mL=((X+Y)/2)×Dil</code>; <code>µM=(mg/mL ÷ MW)×1000</code>.',
    help_res:'<b>Result:</b> always in <b>µM</b>.',
    help_mm:'<b>Disclaimer (MW):</b> Molecular mass is auto-computed from the sequence, but you can provide any value.',
    help_eps:'<b>Coefficients (280 nm):</b> W = <b>5500</b>; Y = <b>1490</b> M⁻¹·cm⁻¹.',
    news_title:'News',
    news_body:'We’re working on new tools, more languages and integrations. Updates will be posted here.',
    news_helipep_title:'New: HeliPep',
    news_helipep_body:'Introducing HeliPep — generate helical wheel plots fast, customize colors/shapes, and export a .zip for multiple sequences.',
    about_title:'About Us', about_sub:'About Troplia',
    about_body:'<b>What is Troplia</b><br>Troplia aims to develop intelligent solutions for society and improve production processes, for example in the food sector, by providing tools to enhance decision-making for producers of all scales, as well as developing products capable of controlling agricultural pests that are environmentally friendly and safer for consumers.<br><br><b>Our Mission</b><br>Our mission is to transform science and technology into practical solutions that positively impact people’s lives and the planet’s future. We strive to combine biotechnology, artificial intelligence, and sustainable innovation to develop products and services that improve production processes, promote greater food and environmental safety, and facilitate strategic decision-making in key sectors such as agriculture and industry. At the same time, we aim to democratize access to analytical and educational tools, empowering producers, researchers, and students to explore new possibilities and contribute to the United Nations’ 2030 Agenda Sustainable Development Goals. We believe science and society go hand in hand, and we work to be a strong link between innovation, social responsibility, and environmental preservation.<br><br><b>Areas of Focus</b><br>Biotechnology, artificial intelligence, agriculture, and education.<br><br><b>Services and Products</b><ul><li>Biotechnological solutions for the agricultural market;</li><li>Development of agricultural products that are safer and environmentally friendly, in line with the UN 2030 Agenda;</li><li>Development of custom tools for academic and industrial environments using neural networks;</li><li>Chromatographic and mass spectrometry analyses;</li><li>Consulting in data analysis and project development;</li><li>Development of training courses in bioprospecting and analytical approaches.</li></ul><br><b>Team</b><ul><li>Dr. Marcelo Henrique Soller Ramada;</li><li>Dr. Stephan Machado Dohms;</li><li>M.Sc. Lorena Ferreira Peixoto;</li><li>M.Sc. Daniel Gusmão de Morais.</li></ul>',
    about_maint_title:'Site Maintenance',
    about_maint_body:'Team in charge:<br><br>Gabriel Iudy Yamaguchi Rocha<br>Marcelo Henrique Soller Ramada<br>Rivaldo Varejão Pasqual Saraiva',
    contact_title:'Contact', contact_intro:'Get in touch for information or partnerships:',
    contact_list:'Gabriel Iudy Yamaguchi Rocha: <a class="link" href="mailto:yamaguchiiudy@gmail.com">yamaguchiiudy@gmail.com</a><br>Marcelo Henrique Soller Ramada: <a class="link" href="mailto:marceloramada@gmail.com">marceloramada@gmail.com</a><br>Rivaldo Varejão Pasqual Saraiva: <a class="link" href="mailto:rivaldo.vps8@gmail.com">rivaldo.vps8@gmail.com</a>',
    helipep_title:'HeliPep — Helical Wheel Plot for peptides',
    helipep_subtitle:'Generate helical wheel plots from a single sequence or a .txt file with multiple sequences.',
    helipep_single_title:'Single sequence',
    helipep_seq_label:'Sequence (letters only):',
    helipep_name_label:'Peptide name (optional; appears in the center):',
    helipep_preview_btn:'Preview',
    helipep_download_btn:'Download PNG',
    helipep_bulk_title:'Multiple sequences (.txt)',
    helipep_bulk_hint:'Example .txt file (one per line):',
    helipep_bulk_btn:'Download result',
    helipep_custom_title:'Customization',
    helipep_color_hydrophobic:'Color hydrophobic',
    helipep_color_basic:'Color basic',
    helipep_color_acidic:'Color acidic',
    helipep_color_polar:'Color polar',
    helipep_shape_hydrophobic:'Shape hydrophobic',
    helipep_shape_basic:'Shape basic',
    helipep_shape_acidic:'Shape acidic',
    helipep_shape_polar:'Shape polar',
    helipep_advanced_title:'Advanced options',
    helipep_refs_title:'References'
  , news_helipep_heading:'Launch of the HeliPep Tool', news_p1:'We are excited to introduce HeliPep — Helical Wheel Plot for Peptides, the newest tool integrated into the Troplia platform.', news_p2:'HeliPep was developed to support researchers seeking to visualize and analyze the structural properties of peptides through helical wheel plots. With it, you can:', news_b1:'<b>Generate instant visualizations</b> of peptides from individual sequences or batch inputs (.txt).', news_b2:'<b>Export high-quality images (PNG)</b>, ready for articles, reports, and presentations.', news_b3:'<b>Customize the colors and shapes</b> of hydrophobic, basic, acidic, and polar residues.', news_b4:'<b>Explore advanced options</b> such as circular guides and hydrophobic moment.', news_p3:'As with QuantyPep, the interface was designed to be <b>simple, intuitive, and responsive</b>. HeliPep complements Troplia by offering a powerful solution for peptide visual analysis — fully aligned with the requirements of high-impact publications.', news_cta:'Try the tool now in the main menu of the platform! 🚀', helipep_select_file_btn:'Select file', helipep_file_none:'No file selected', helipep_opt_angle:'Δ°/residue', helipep_opt_start:'Starting angle', helipep_opt_clockwise:'Clockwise', helipep_opt_layer_n:'Layer every N residues', helipep_opt_links:'Connect residues', helipep_opt_guide:'Guide circle', helipep_opt_nums:'Numbers', helipep_opt_letters:'Letters', helipep_opt_moment:'Hydrophobic moment', label_conc:'Concentration', unit_micromolar:'µM'},

  es:{
    tab_home:'Inicio', tab_hist:'Historial', tab_help:'Instrucciones', tab_news:'Novedades', tab_about:'Acerca de', tab_contact:'Contacto',
    home_title:'Bienvenido a Troplia',
    home_intro:'Esta es la plataforma de Troplia para acceder a diversas herramientas de laboratorio. El sitio se actualizará con nuevas herramientas; actualmente incluye <b>QuantyPep</b> para concentración de péptidos por absorbancia: simple, rápida y precisa.',
    tools_title:'Herramientas',
    tool_calc_title:'QuantyPep – Calculadora de Concentración de Péptidos (Absorbancia)',
    tool_calc_desc:'Calcule concentraciones por A280 (Con W o Con YYY) o por A205/215/225 (Sin W). Resultados siempre en micromolar (µM).',
    tool_helipep_title:'HeliPep — Helical Wheel Plot para Péptidos',
    tool_helipep_desc:'Helical wheel plots al instante (vista previa + descarga), soporte .txt y personalización de colores/formas.',
    calc_title:'QuantyPep — Calculadora de Concentración de Péptidos (Absorbancia)',
    proj_info:'Información del Proyecto', lbl_project:'Nombre del Proyecto', lbl_date:'Fecha', lbl_peptide_name:'Nombre del Péptido', lbl_sequence:'Secuencia del Péptido', opt_field:'(opcional)',
    ph_project:'Mi experimento', ph_peptide_name:'Ej.: Pep-1', ph_sequence:'Ej.: ACDEFGHIKLMNPQRSTVWY',
    calc_type:'Tipo de Cálculo', btn_compute:'Calcular', btn_clear:'Limpiar', btn_export_csv:'Exportar CSV', btn_save_txt:'Guardar TXT', result:'Resultado',
    calc_opt_auto:'Auto', calc_opt_w:'Con W (A280)', calc_opt_yyy:'Con YYY (A280)', calc_opt_now:'Sin W (205/215/225)',
    calc_name_w:'Con W', calc_name_yyy:'Con YYY', calc_name_now:'Sin W',
    tip_w:'A280 (W+Y): ε_total = (#W×5500 + #Y×1490).', tip_yyy:'A280 (W+Y): ε_total = (#W×5500 + #Y×1490).', tip_now:'Usa 205/215/225 nm.',
    obs_w:'Modo Con W (A280, W+Y) aplicado.', obs_yyy:'Modo Con YYY (A280, W+Y) aplicado.', obs_now:'Modo Sin W (205/215/225) aplicado.',
    hist_title:'Historial', th_project:'Proyecto', th_date:'Fecha (YYYY-MM-DD)', th_peptide:'Péptido', th_type:'Tipo', th_dil:'Dilución', th_obs:'Notas',
    btn_clear_history:'Borrar historial',
    help_title:'Instrucciones',
    help_step:'<b>Pasos:</b> Indique <i>Nombre del Proyecto</i>, <i>Fecha</i> (YYYY-MM-DD), <i>Nombre del Péptido</i> y, opcionalmente, la <i>Secuencia</i>.',
    help_auto:'<b>Auto:</b> si la secuencia tiene <b>W ≥ 1</b> → <i>Con W</i>; si <b>W = 0</b> y <b>Y ≥ 3</b> → <i>Con YYY</i>; de lo contrario → <i>Sin W</i>.<br><b>Aviso:</b> <i>Con W</i> y <i>Con YYY</i> usan la <u>misma calculadora A280</u>. Para A280: <b>ε_total = (#W×5500 + #Y×1490)</b>.',
    help_w:'<b>Con W / Con YYY (A280):</b> <code>ε_total = #W×5500 + #Y×1490</code> (M⁻¹·cm⁻¹).<br><code>µM = (A280 / ε_total) × Dil × 1.000.000</code>.',
    help_yyy:'<b>Nota:</b> Péptidos con ≥3 tirosinas y sin W usan la misma fórmula A280; con W y Y, ambos contribuyen vía ε_total.',
    help_now:'<b>Sin W (205/215/225):</b> <code>X=(A215−A225)×144</code>; <code>Y=A205×31</code>; <code>mg/mL=((X+Y)/2)×Dil</code>; <code>µM=(mg/mL ÷ MM)×1000</code>.',
    help_res:'<b>Resultado:</b> siempre en <b>µM</b>.',
    help_mm:'<b>Aviso (MM):</b> La masa molar se calcula automáticamente por la secuencia, pero puede indicarse cualquier valor.',
    help_eps:'<b>Coeficientes (280 nm):</b> W = <b>5500</b>; Y = <b>1490</b> M⁻¹·cm⁻¹.',
    news_title:'Novedades',
    news_body:'Estamos trabajando en nuevas herramientas, más idiomas e integraciones. Publicaremos novedades aquí.',
    news_helipep_title:'Nuevo: HeliPep',
    news_helipep_body:'Presentamos HeliPep — genere helical wheel plots con rapidez, personalice colores/formas y exporte un .zip para múltiples secuencias.',
    about_title:'Acerca de', about_sub:'Sobre Troplia',
    about_body:'Bajo la coordinación del Prof. Marcelo Henrique Soller Ramada, formamos parte del Programa de Posgrado en Ciencias Genómicas y Biotecnología de la Universidad Católica de Brasilia (UCB). Contamos con apoyo de CAPES, CNPq y FAPDF, utilizando la infraestructura de la UCB.',
    about_maint_title:'Mantenimiento del sitio',
    about_maint_body:'Equipo responsable:<br><br>Gabriel Iudy Yamaguchi Rocha<br>Marcelo Henrique Soller Ramada<br>Rivaldo Varejão Pasqual Saraiva',
    contact_title:'Contacto', contact_intro:'Contáctenos para información o propuestas de colaboración:',
    contact_list:'Gabriel Iudy Yamaguchi Rocha: <a class="link" href="mailto:yamaguchiiudy@gmail.com">yamaguchiiudy@gmail.com</a><br>Marcelo Henrique Soller Ramada: <a class="link" href="mailto:marceloramada@gmail.com">marceloramada@gmail.com</a><br>Rivaldo Varejão Pasqual Saraiva: <a class="link" href="mailto:rivaldo.vps8@gmail.com">rivaldo.vps8@gmail.com</a>',
    helipep_title:'HeliPep — Helical Wheel Plot para péptidos',
    helipep_subtitle:'Genere helical wheel plots a partir de una secuencia o de un archivo .txt con múltiples secuencias.',
    helipep_single_title:'Secuencia única',
    helipep_seq_label:'Secuencia (solo letras):',
    helipep_name_label:'Nombre del péptido (opcional; aparece en el centro):',
    helipep_preview_btn:'Vista previa',
    helipep_download_btn:'Descargar PNG',
    helipep_bulk_title:'Múltiples secuencias (.txt)',
    helipep_bulk_hint:'Ejemplo de archivo .txt (una por línea):',
    helipep_bulk_btn:'Generar .zip',
    helipep_custom_title:'Personalización',
    helipep_color_hydrophobic:'Color hidrofóbicos',
    helipep_color_basic:'Color básicos',
    helipep_color_acidic:'Color ácidos',
    helipep_color_polar:'Color polares',
    helipep_shape_hydrophobic:'Forma hidrofóbicos',
    helipep_shape_basic:'Forma básicos',
    helipep_shape_acidic:'Forma ácidos',
    helipep_shape_polar:'Forma polares',
    helipep_advanced_title:'Opciones avanzadas',
    helipep_refs_title:'Referencias'
  },

  fr:{
    tab_home:'Accueil', tab_hist:'Historique', tab_help:'Instructions', tab_news:'Actualités', tab_about:'À propos', tab_contact:'Contact',
    home_title:'Bienvenue sur Troplia',
    home_intro:'Voici la plateforme Troplia pour accéder à plusieurs outils de laboratoire. De nouveaux outils seront ajoutés ; pour l’instant, nous proposons <b>QuantyPep</b> pour la concentration de peptides par absorbance — simple, rapide et précise.',
    tools_title:'Outils',
    tool_calc_title:'QuantyPep – Calculatrice de Concentration de Peptides (Absorbance)',
    tool_calc_desc:'Calculez les concentrations par A280 (Avec W ou Avec YYY) ou par A205/215/225 (Sans W). Résultats toujours en micromolaire (µM).',
    tool_helipep_title:'HeliPep — Helical Wheel Plot pour Peptides',
    tool_helipep_desc:'Plots instantanés (aperçu + téléchargement), support .txt et personnalisation couleurs/formes.',
    calc_title:'QuantyPep — Calculatrice de Concentration de Peptides (Absorbance)',
    proj_info:'Informations du projet', lbl_project:'Nom du projet', lbl_date:'Date', lbl_peptide_name:'Nom du peptide', lbl_sequence:'Séquence du peptide', opt_field:'(optionnel)',
    ph_project:'Mon expérience', ph_peptide_name:'Ex. : Pep-1', ph_sequence:'Ex. : ACDEFGHIKLMNPQRSTVWY',
    calc_type:'Type de calcul', btn_compute:'Calculer', btn_clear:'Effacer', btn_export_csv:'Exporter CSV', btn_save_txt:'Enregistrer TXT', result:'Résultat',
    calc_opt_auto:'Auto', calc_opt_w:'Avec W (A280)', calc_opt_yyy:'Avec YYY (A280)', calc_opt_now:'Sans W (205/215/225)',
    calc_name_w:'Avec W', calc_name_yyy:'Avec YYY', calc_name_now:'Sans W',
    tip_w:'A280 (W+Y) : ε_total = (#W×5500 + #Y×1490).', tip_yyy:'A280 (W+Y) : ε_total = (#W×5500 + #Y×1490).', tip_now:'Utilise 205/215/225 nm.',
    obs_w:'Mode « Avec W » (A280, W+Y) appliqué.', obs_yyy:'Mode « Avec YYY » (A280, W+Y) appliqué.', obs_now:'Mode « Sans W » (205/215/225) appliqué.',
    hist_title:'Historique', th_project:'Projet', th_date:'Date (YYYY-MM-DD)', th_peptide:'Peptide', th_type:'Type', th_dil:'Dilution', th_obs:'Remarque',
    btn_clear_history:'Vider l’historique',
    help_title:'Instructions',
    help_step:'<b>Étapes :</b> Indiquez <i>Nom du projet</i>, <i>Date</i> (YYYY-MM-DD), <i>Nom du peptide</i> et éventuellement la <i>Séquence</i>.',
    help_auto:'<b>Auto :</b> si la séquence a <b>W ≥ 1</b> → <i>Avec W</i> ; si <b>W = 0</b> et <b>Y ≥ 3</b> → <i>Avec YYY</i> ; sinon → <i>Sans W</i>.<br><b>Note :</b> <i>Avec W</i> et <i>Avec YYY</i> utilisent la <u>même calculatrice A280</u>.',
    help_w:'<b>Avec W / Avec YYY (A280) :</b> <code>ε_total = #W×5500 + #Y×1490</code> (M⁻¹·cm⁻¹).<br><code>µM = (A280 / ε_total) × Dil × 1 000 000</code>.',
    help_yyy:'<b>Remarque :</b> ≥3 tyrosines sans W : même formule A280 ; avec W+Y, les deux contribuent via ε_total.',
    help_now:'<b>Sans W (205/215/225) :</b> <code>X=(A215−A225)×144</code> ; <code>Y=A205×31</code> ; <code>mg/mL=((X+Y)/2)×Dil</code> ; <code>µM=(mg/mL ÷ MM)×1000</code>.',
    help_res:'<b>Résultat :</b> toujours en <b>µM</b>.',
    help_mm:'<b>Note (MM) :</b> La masse molaire est calculée automatiquement ; vous pouvez saisir une valeur.',
    help_eps:'<b>Coefficients (280 nm) :</b> W = <b>5500</b> ; Y = <b>1490</b> M⁻¹·cm⁻¹.',
    news_title:'Actualités',
    news_body:'Nous travaillons sur de nouveaux outils, plus de langues et des intégrations. Des mises à jour seront publiées ici.',
    news_helipep_title:'Nouveau : HeliPep',
    news_helipep_body:'Découvrez HeliPep — wheel plots rapides, personnalisez couleurs/formes et exportez un .zip pour plusieurs séquences.',
    about_title:'À propos', about_sub:'À propos de Troplia',
    about_body:'Sous la direction du Prof. Marcelo Henrique Soller Ramada, nous faisons partie du programme de master en sciences génomiques et biotechnologie de l’UCB. Partenariats stratégiques et soutien CAPES, CNPq, FAPDF.',
    about_maint_title:'Maintenance du site',
    about_maint_body:'Équipe responsable :<br><br>Gabriel Iudy Yamaguchi Rocha<br>Marcelo Henrique Soller Ramada<br>Rivaldo Varejão Pasqual Saraiva',
    contact_title:'Contact', contact_intro:'Contactez-nous pour toute information ou partenariat :',
    contact_list:'Gabriel Iudy Yamaguchi Rocha : <a class="link" href="mailto:yamaguchiiudy@gmail.com">yamaguchiiudy@gmail.com</a><br>Marcelo Henrique Soller Ramada : <a class="link" href="mailto:marceloramada@gmail.com">marceloramada@gmail.com</a><br>Rivaldo Varejão Pasqual Saraiva : <a class="link" href="mailto:rivaldo.vps8@gmail.com">rivaldo.vps8@gmail.com</a>',
    helipep_title:'HeliPep — Helical Wheel Plot pour peptides',
    helipep_subtitle:'Générez des helical wheel plots à partir d’une séquence unique ou d’un fichier .txt contenant plusieurs séquences.',
    helipep_single_title:'Séquence unique',
    helipep_seq_label:'Séquence (lettres uniquement) :',
    helipep_name_label:'Nom du peptide (optionnel ; au centre) :',
    helipep_preview_btn:'Aperçu',
    helipep_download_btn:'Télécharger PNG',
    helipep_bulk_title:'Plusieurs séquences (.txt)',
    helipep_bulk_hint:'Exemple de fichier .txt (une par ligne) :',
    helipep_bulk_btn:'Générer .zip',
    helipep_custom_title:'Personnalisation',
    helipep_color_hydrophobic:'Couleur hydrophobes',
    helipep_color_basic:'Couleur basiques',
    helipep_color_acidic:'Couleur acides',
    helipep_color_polar:'Couleur polaires',
    helipep_shape_hydrophobic:'Forme hydrophobes',
    helipep_shape_basic:'Forme basiques',
    helipep_shape_acidic:'Forme acides',
    helipep_shape_polar:'Forme polaires',
    helipep_advanced_title:'Options avancées',
    helipep_refs_title:'Références'
  },

  de:{
    tab_home:'Start', tab_hist:'Verlauf', tab_help:'Anleitung', tab_news:'Neuigkeiten', tab_about:'Über uns', tab_contact:'Kontakt',
    home_title:'Willkommen bei Troplia',
    home_intro:'Dies ist die Troplia-Plattform für verschiedene Laborwerkzeuge. Neue Tools folgen; derzeit gibt es <b>QuantyPep</b> zur Peptidkonzentration per Absorption — einfach, schnell, präzise.',
    tools_title:'Werkzeuge',
    tool_calc_title:'QuantyPep – Peptid-Konzentrationsrechner (Absorption)',
    tool_calc_desc:'Konzentration per A280 (Mit W oder Mit YYY) oder A205/215/225 (Ohne W). Ergebnis stets in Mikromolar (µM).',
    tool_helipep_title:'HeliPep — Helical Wheel Plot für Peptide',
    tool_helipep_desc:'Sofortige Wheel-Plots (Vorschau + Download), .txt-Batch und Farb-/Formanpassung.',
    calc_title:'QuantyPep — Peptid-Konzentrationsrechner (Absorption)',
    proj_info:'Projektinfo', lbl_project:'Projektname', lbl_date:'Datum', lbl_peptide_name:'Peptidname', lbl_sequence:'Peptidsequenz', opt_field:'(optional)',
    ph_project:'Mein Experiment', ph_peptide_name:'z. B. Pep-1', ph_sequence:'z. B. ACDEFGHIKLMNPQRSTVWY',
    calc_type:'Rechentyp', btn_compute:'Berechnen', btn_clear:'Leeren', btn_export_csv:'CSV exportieren', btn_save_txt:'TXT speichern', result:'Ergebnis',
    calc_opt_auto:'Auto', calc_opt_w:'Mit W (A280)', calc_opt_yyy:'Mit YYY (A280)', calc_opt_now:'Ohne W (205/215/225)',
    calc_name_w:'Mit W', calc_name_yyy:'Mit YYY', calc_name_now:'Ohne W',
    tip_w:'A280 (W+Y): ε_total = (#W×5500 + #Y×1490).', tip_yyy:'A280 (W+Y): ε_total = (#W×5500 + #Y×1490).', tip_now:'Verwendet 205/215/225 nm.',
    obs_w:'Modus „Mit W“ (A280, W+Y) angewendet.', obs_yyy:'Modus „Mit YYY“ (A280, W+Y) angewendet.', obs_now:'Modus „Ohne W“ (205/215/225) angewendet.',
    hist_title:'Verlauf', th_project:'Projekt', th_date:'Datum (YYYY-MM-DD)', th_peptide:'Peptid', th_type:'Typ', th_dil:'Verdünnung', th_obs:'Hinweis',
    btn_clear_history:'Verlauf löschen',
    help_title:'Anleitung',
    help_step:'<b>Vorgehen:</b> <i>Projektname</i>, <i>Datum</i> (YYYY-MM-DD), <i>Peptidname</i> und optional die <i>Sequenz</i> angeben.',
    help_auto:'<b>Auto:</b> <b>W ≥ 1</b> → <i>Mit W</i>; <b>W = 0</b> & <b>Y ≥ 3</b> → <i>Mit YYY</i>; sonst → <i>Ohne W</i>. <u>Mit W</u> und <u>Mit YYY</u> nutzen denselben A280-Rechner.',
    help_w:'<b>Mit W / Mit YYY (A280):</b> <code>ε_total = #W×5500 + #Y×1490</code>.',
    help_yyy:'<b>Hinweis:</b> ≥3 Tyrosine ohne W: gleiche A280-Formel; mit W+Y tragen beide bei.',
    help_now:'<b>Ohne W (205/215/225):</b> <code>X=(A215−A225)×144</code>; <code>Y=A205×31</code>; <code>mg/mL=((X+Y)/2)×Verdünnung</code>; <code>µM=(mg/mL ÷ MM)×1000</code>.',
    help_res:'<b>Ergebnis:</b> immer in <b>µM</b>.',
    help_mm:'<b>Hinweis (MM):</b> Molmasse wird aus der Sequenz berechnet; Wert kann überschrieben werden.',
    help_eps:'<b>Koeffizienten (280 nm):</b> W = <b>5500</b>; Y = <b>1490</b>.',
    news_title:'Neuigkeiten',
    news_body:'Wir arbeiten an neuen Tools, mehr Sprachen und Integrationen. Updates folgen hier.',
    news_helipep_title:'Neu: HeliPep',
    news_helipep_body:'Neu: HeliPep — Wheel-Plots schnell erzeugen, Farben/Formen anpassen und ZIP für mehrere Sequenzen exportieren.',
    about_title:'Über uns', about_sub:'Über Troplia',
    about_body:'Unter Leitung von Prof. Marcelo Henrique Soller Ramada am Graduiertenprogramm der UCB (Genomik & Biotechnologie). Unterstützt von CAPES, CNPq, FAPDF.',
    about_maint_title:'Wartung der Website',
    about_maint_body:'Verantwortliches Team:<br><br>Gabriel Iudy Yamaguchi Rocha<br>Marcelo Henrique Soller Ramada<br>Rivaldo Varejão Pasqual Saraiva',
    contact_title:'Kontakt', contact_intro:'Kontakt für Infos oder Partnerschaften:',
    contact_list:'Gabriel Iudy Yamaguchi Rocha: <a class="link" href="mailto:yamaguchiiudy@gmail.com">yamaguchiiudy@gmail.com</a><br>Marcelo Henrique Soller Ramada: <a class="link" href="mailto:marceloramada@gmail.com">marceloramada@gmail.com</a><br>Rivaldo Varejão Pasqual Saraiva: <a class="link" href="mailto:rivaldo.vps8@gmail.com">rivaldo.vps8@gmail.com</a>',
    helipep_title:'HeliPep — Helical Wheel Plot für Peptide',
    helipep_subtitle:'Erstellen Sie Helical-Wheel-Plots aus einer einzelnen Sequenz oder aus einer .txt-Datei mit mehreren Sequenzen.',
    helipep_single_title:'Einzelne Sequenz',
    helipep_seq_label:'Sequenz (nur Buchstaben):',
    helipep_name_label:'Peptidname (optional; in der Mitte):',
    helipep_preview_btn:'Vorschau',
    helipep_download_btn:'PNG herunterladen',
    helipep_bulk_title:'Mehrere Sequenzen (.txt)',
    helipep_bulk_hint:'Beispiel .txt-Datei (eine pro Zeile):',
    helipep_bulk_btn:'ZIP erzeugen',
    helipep_custom_title:'Anpassung',
    helipep_color_hydrophobic:'Farbe hydrophob',
    helipep_color_basic:'Farbe basisch',
    helipep_color_acidic:'Farbe sauer',
    helipep_color_polar:'Farbe polar',
    helipep_shape_hydrophobic:'Form hydrophob',
    helipep_shape_basic:'Form basisch',
    helipep_shape_acidic:'Form sauer',
    helipep_shape_polar:'Form polar',
    helipep_advanced_title:'Erweiterte Optionen',
    helipep_refs_title:'Referenzen'
  },

  it:{
    tab_home:'Home', tab_hist:'Storico', tab_help:'Istruzioni', tab_news:'Novità', tab_about:'Chi siamo', tab_contact:'Contatti',
    home_title:'Benvenuto su Troplia',
    home_intro:'Questa è la piattaforma Troplia per vari strumenti di laboratorio. Il sito verrà aggiornato con nuovi strumenti; al momento include <b>QuantyPep</b> per la concentrazione dei peptidi per assorbanza — semplice, veloce e precisa.',
    tools_title:'Strumenti',
    tool_calc_title:'QuantyPep – Calcolatore concentrazione peptidi (Assorbanza)',
    tool_calc_desc:'Calcola per A280 (Con W o Con YYY) oppure A205/215/225 (Senza W). Risultati sempre in micromolare (µM).',
    tool_helipep_title:'HeliPep — Helical Wheel Plot per Peptidi',
    tool_helipep_desc:'Wheel plot istantanei (anteprima + download), supporto .txt e personalizzazione colori/forme.',
    calc_title:'QuantyPep — Calcolatore concentrazione peptidi (Assorbanza)',
    proj_info:'Info progetto', lbl_project:'Nome progetto', lbl_date:'Data', lbl_peptide_name:'Nome peptide', lbl_sequence:'Sequenza peptide', opt_field:'(opzionale)',
    ph_project:'Il mio esperimento', ph_peptide_name:'Es.: Pep-1', ph_sequence:'Es.: ACDEFGHIKLMNPQRSTVWY',
    calc_type:'Tipo di calcolo', btn_compute:'Calcola', btn_clear:'Pulisci', btn_export_csv:'Esporta CSV', btn_save_txt:'Salva TXT', result:'Risultato',
    calc_opt_auto:'Auto', calc_opt_w:'Con W (A280)', calc_opt_yyy:'Con YYY (A280)', calc_opt_now:'Senza W (205/215/225)',
    calc_name_w:'Con W', calc_name_yyy:'Con YYY', calc_name_now:'Senza W',
    tip_w:'A280 (W+Y): ε_total = (#W×5500 + #Y×1490).',
    tip_yyy:'A280 (W+Y): ε_total = (#W×5500 + #Y×1490).',
    tip_now:'Usa 205/215/225 nm.',
    obs_w:'Modalità "Con W" (A280, W+Y) applicata.', obs_yyy:'Modalità "Con YYY" (A280, W+Y) applicata.', obs_now:'Modalità "Senza W" (205/215/225) applicata.',
    hist_title:'Storico', th_project:'Progetto', th_date:'Data (YYYY-MM-DD)', th_peptide:'Peptide', th_type:'Tipo', th_dil:'Diluizione', th_obs:'Nota',
    btn_clear_history:'Svuota storico',
    help_title:'Istruzioni',
    help_step:'<b>Passi:</b> Inserisci <i>Nome progetto</i>, <i>Data</i> (YYYY-MM-DD), <i>Nome peptide</i> e, opzionalmente, la <i>Sequenza</i>.',
    help_auto:'<b>Auto:</b> se <b>W ≥ 1</b> → <i>Con W</i>; se <b>W = 0</b> e <b>Y ≥ 3</b> → <i>Con YYY</i>; altrimenti → <i>Senza W</i>. Stessa calcolatrice A280 per W/YYY.',
    help_w:'<b>Con W / Con YYY (A280):</b> <code>ε_total = #W×5500 + #Y×1490</code>.',
    help_yyy:'<b>Nota:</b> ≥3 tirosine senza W: stessa formula A280; con W+Y contribuiscono entrambi.',
    help_now:'<b>Senza W (205/215/225):</b> <code>X=(A215−A225)×144</code>; <code>Y=A205×31</code>; <code>mg/mL=((X+Y)/2)×Diluizione</code>; <code>µM=(mg/mL ÷ MM)×1000</code>.',
    help_res:'<b>Risultato:</b> sempre in <b>µM</b>.',
    help_mm:'<b>Nota (MM):</b> La massa molecolare è calcolata automaticamente, ma puoi inserirne una.',
    help_eps:'<b>Coefficienti (280 nm):</b> W = <b>5500</b>; Y = <b>1490</b>.',
    news_title:'Novità',
    news_body:'Stiamo lavorando su nuovi strumenti, più lingue e integrazioni. Gli aggiornamenti appariranno qui.',
    news_helipep_title:'Novità: HeliPep',
    news_helipep_body:'Presentiamo HeliPep — wheel plot rapidi, personalizzazione colori/forme, export .zip per più sequenze.',
    about_title:'Chi siamo', about_sub:'Su Troplia',
    about_body:'Sotto la guida del Prof. Marcelo Henrique Soller Ramada (UCB), con supporto CAPES, CNPq, FAPDF.',
    about_maint_title:'Manutenzione del sito',
    about_maint_body:'Team responsabile:<br><br>Gabriel Iudy Yamaguchi Rocha<br>Marcelo Henrique Soller Ramada<br>Rivaldo Varejão Pasqual Saraiva',
    contact_title:'Contatti', contact_intro:'Contattaci per informazioni o collaborazioni:',
    contact_list:'Gabriel Iudy Yamaguchi Rocha: <a class="link" href="mailto:yamaguchiiudy@gmail.com">yamaguchiiudy@gmail.com</a><br>Marcelo Henrique Soller Ramada: <a class="link" href="mailto:marceloramada@gmail.com">marceloramada@gmail.com</a><br>Rivaldo Varejão Pasqual Saraiva: <a class="link" href="mailto:rivaldo.vps8@gmail.com">rivaldo.vps8@gmail.com</a>',
    helipep_title:'HeliPep — Helical Wheel Plot per peptidi',
    helipep_subtitle:'Genera helical wheel plot da una singola sequenza o da un file .txt con più sequenze.',
    helipep_single_title:'Sequenza singola',
    helipep_seq_label:'Sequenza (solo lettere):',
    helipep_name_label:'Nome del peptide (opzionale; al centro):',
    helipep_preview_btn:'Anteprima',
    helipep_download_btn:'Scarica PNG',
    helipep_bulk_title:'Sequenze multiple (.txt)',
    helipep_bulk_hint:'Esempio di file .txt (una per riga):',
    helipep_bulk_btn:'Genera .zip',
    helipep_custom_title:'Personalizzazione',
    helipep_color_hydrophobic:'Colore idrofobici',
    helipep_color_basic:'Colore basici',
    helipep_color_acidic:'Colore acidi',
    helipep_color_polar:'Colore polari',
    helipep_shape_hydrophobic:'Forma idrofobici',
    helipep_shape_basic:'Forma basici',
    helipep_shape_acidic:'Forma acidi',
    helipep_shape_polar:'Forma polari',
    helipep_advanced_title:'Opzioni avanzate',
    helipep_refs_title:'Riferimenti'
  },

  ja:{
    tab_home:'ホーム', tab_hist:'履歴', tab_help:'手順', tab_news:'ニュース', tab_about:'概要', tab_contact:'連絡先',
    home_title:'Troplia へようこそ',
    home_intro:'本プラットフォームでは研究室向けツールを提供します。順次追加予定。現在は<b>QuantyPep</b>（吸光度によるペプチド濃度計算）を提供しています。',
    tools_title:'ツール',
    tool_calc_title:'QuantyPep – ペプチド濃度計算（吸光度）',
    tool_calc_desc:'A280（W あり / YYY あり）または A205/215/225（W なし）で計算。結果は常に µM。',
    tool_helipep_title:'HeliPep — ペプチド用 Helical Wheel Plot',
    tool_helipep_desc:'即時プロット（プレビュー＋ダウンロード）、.txt バッチ、色/形のカスタマイズ対応。',
    calc_title:'QuantyPep — ペプチド濃度計算（吸光度）',
    proj_info:'プロジェクト情報', lbl_project:'プロジェクト名', lbl_date:'日付', lbl_peptide_name:'ペプチド名', lbl_sequence:'配列', opt_field:'（任意）',
    ph_project:'私の実験', ph_peptide_name:'例: Pep-1', ph_sequence:'例: ACDEFGHIKLMNPQRSTVWY',
    calc_type:'計算タイプ', btn_compute:'計算', btn_clear:'クリア', btn_export_csv:'CSV エクスポート', btn_save_txt:'TXT を保存', result:'結果',
    calc_opt_auto:'自動', calc_opt_w:'W あり (A280)', calc_opt_yyy:'YYY あり (A280)', calc_opt_now:'W なし (205/215/225)',
    calc_name_w:'W あり', calc_name_yyy:'YYY あり', calc_name_now:'W なし',
    tip_w:'A280 (W+Y): ε_total = (#W×5500 + #Y×1490).',
    tip_yyy:'A280 (W+Y): ε_total = (#W×5500 + #Y×1490).',
    tip_now:'205/215/225 nm を使用。',
    obs_w:'「W あり」（A280, W+Y）モードを適用。', obs_yyy:'「YYY あり」（A280, W+Y）モードを適用。', obs_now:'「W なし」（205/215/225）モードを適用。',
    hist_title:'履歴', th_project:'プロジェクト', th_date:'日付 (YYYY-MM-DD)', th_peptide:'ペプチド', th_type:'タイプ', th_dil:'希釈', th_obs:'注記',
    btn_clear_history:'履歴を消去',
    help_title:'手順',
    help_step:'<b>手順:</b> <i>プロジェクト名</i>, <i>日付</i> (YYYY-MM-DD), <i>ペプチド名</i>、任意で<i>配列</i>を入力。',
    help_auto:'<b>自動:</b> <b>W ≥ 1</b> → <i>W あり</i>； <b>W = 0</b> かつ <b>Y ≥ 3</b> → <i>YYY あり</i>；それ以外 → <i>W なし</i>。W/YYY は同じ A280 計算を使用。',
    help_w:'<b>W あり / YYY あり (A280):</b> <code>ε_total = #W×5500 + #Y×1490</code>.',
    help_yyy:'<b>注:</b> W なしで Y≥3 の場合も A280 公式は同じ。W と Y がある場合は両方寄与。',
    help_now:'<b>W なし (205/215/225):</b> <code>X=(A215−A225)×144</code>、<code>Y=A205×31</code>、<code>mg/mL=((X+Y)/2)×希釈</code>、<code>µM=(mg/mL ÷ 分子量)×1000</code>。',
    help_res:'<b>結果:</b> 常に <b>µM</b>。',
    help_mm:'<b>注 (分子量):</b> 配列から自動計算しますが手動入力も可能です。',
    help_eps:'<b>係数 (280 nm):</b> W = <b>5500</b>; Y = <b>1490</b>。',
    news_title:'ニュース',
    news_body:'新しいツール、多言語、連携機能を順次リリースします。更新情報はこちらに掲載します。',
    news_helipep_title:'新機能: HeliPep',
    news_helipep_body:'HeliPep 公開 — Helical Wheel Plot を素早く生成、色/形をカスタマイズし、複数配列を .zip でエクスポート。',
    about_title:'概要', about_sub:'Troplia について',
    about_body:'UCB の Marcelo Henrique Soller Ramada 教授の指導の下、CAPES/CNPq/FAPDF の支援を受けています。',
    about_maint_title:'サイト運用',
    about_maint_body:'担当チーム：<br><br>Gabriel Iudy Yamaguchi Rocha<br>Marcelo Henrique Soller Ramada<br>Rivaldo Varejão Pasqual Saraiva',
    contact_title:'連絡先', contact_intro:'情報・連携のご相談は以下まで：',
    contact_list:'Gabriel Iudy Yamaguchi Rocha: <a class="link" href="mailto:yamaguchiiudy@gmail.com">yamaguchiiudy@gmail.com</a><br>Marcelo Henrique Soller Ramada: <a class="link" href="mailto:marceloramada@gmail.com">marceloramada@gmail.com</a><br>Rivaldo Varejão Pasqual Saraiva: <a class="link" href="mailto:rivaldo.vps8@gmail.com">rivaldo.vps8@gmail.com</a>',
    helipep_title:'HeliPep — ペプチド用 Helical Wheel Plot',
    helipep_subtitle:'単一配列、または複数配列を含む .txt ファイルから Helical Wheel Plot を作成します。',
    helipep_single_title:'単一配列',
    helipep_seq_label:'配列（英字のみ）:',
    helipep_name_label:'ペプチド名（任意；中央に表示）:',
    helipep_preview_btn:'プレビュー',
    helipep_download_btn:'PNG をダウンロード',
    helipep_bulk_title:'複数配列（.txt）',
    helipep_bulk_hint:'例（1 行につき 1 配列）:',
    helipep_bulk_btn:'.zip を作成',
    helipep_custom_title:'カスタマイズ',
    helipep_color_hydrophobic:'疎水性の色',
    helipep_color_basic:'塩基性の色',
    helipep_color_acidic:'酸性の色',
    helipep_color_polar:'極性の色',
    helipep_shape_hydrophobic:'疎水性の形',
    helipep_shape_basic:'塩基性の形',
    helipep_shape_acidic:'酸性の形',
    helipep_shape_polar:'極性の形',
    helipep_advanced_title:'詳細オプション',
    helipep_refs_title:'参考文献'
  },

  zh:{
    tab_home:'首页', tab_hist:'历史', tab_help:'说明', tab_news:'新闻', tab_about:'关于我们', tab_contact:'联系',
    home_title:'欢迎来到 Troplia',
    home_intro:'这是 Troplia 的实验室工具平台。我们将持续更新新工具；目前提供 <b>QuantyPep</b>（基于吸收的肽浓度计算）—— 简单、快速、准确。',
    tools_title:'工具',
    tool_calc_title:'QuantyPep – 肽浓度计算器（吸收）',
    tool_calc_desc:'可用 A280（含 W 或 含 YYY）或 A205/215/225（无 W）计算。结果为微摩尔（µM）。',
    tool_helipep_title:'HeliPep — 肽类 Helical Wheel Plot',
    tool_helipep_desc:'即时生成（预览+下载）、支持 .txt 批量、颜色/形状自定义。',
    calc_title:'QuantyPep — 肽浓度计算器（吸收）',
    proj_info:'项目信息', lbl_project:'项目名称', lbl_date:'日期', lbl_peptide_name:'肽名称', lbl_sequence:'肽序列', opt_field:'（可选）',
    ph_project:'我的实验', ph_peptide_name:'如：Pep-1', ph_sequence:'如：ACDEFGHIKLMNPQRSTVWY',
    calc_type:'计算类型', btn_compute:'计算', btn_clear:'清除', btn_export_csv:'导出 CSV', btn_save_txt:'保存 TXT', result:'结果',
    calc_opt_auto:'自动', calc_opt_w:'含 W (A280)', calc_opt_yyy:'含 YYY (A280)', calc_opt_now:'无 W (205/215/225)',
    calc_name_w:'含 W', calc_name_yyy:'含 YYY', calc_name_now:'无 W',
    tip_w:'A280 (W+Y): ε_total = (#W×5500 + #Y×1490).', tip_yyy:'A280 (W+Y): ε_total = (#W×5500 + #Y×1490).', tip_now:'使用 205/215/225 nm。',
    obs_w:'已应用“含 W”（A280，W+Y）模式。', obs_yyy:'已应用“含 YYY”（A280，W+Y）模式。', obs_now:'已应用“无 W”（205/215/225）模式。',
    hist_title:'历史', th_project:'项目', th_date:'日期 (YYYY-MM-DD)', th_peptide:'肽', th_type:'类型', th_dil:'稀释倍数', th_obs:'备注',
    btn_clear_history:'清空历史',
    help_title:'说明',
    help_step:'<b>步骤：</b> 填写<i>项目名称</i>、<i>日期</i> (YYYY-MM-DD)、<i>肽名称</i>，以及可选的<i>序列</i>。',
    help_auto:'<b>自动：</b> 若 <b>W ≥ 1</b> → <i>含 W</i>；若 <b>W = 0</b> 且 <b>Y ≥ 3</b> → <i>含 YYY</i>；否则 → <i>无 W</i>。W/YYY 使用同一 A280 计算。',
    help_w:'<b>含 W / 含 YYY (A280)：</b> <code>ε_total = #W×5500 + #Y×1490</code>。',
    help_yyy:'<b>注：</b> 无 W 且 Y≥3 时同一 A280 公式；有 W 与 Y 时二者共同贡献。',
    help_now:'<b>无 W (205/215/225)：</b> <code>X=(A215−A225)×144</code>；<code>Y=A205×31</code>；<code>mg/mL=((X+Y)/2)×稀释</code>；<code>µM=(mg/mL ÷ 相对分子质量)×1000</code>。',
    help_res:'<b>结果：</b> 始终为 <b>µM</b>。',
    help_mm:'<b>说明（分子质量）：</b> 可自动计算，也可手动输入。',
    help_eps:'<b>系数 (280 nm)：</b> W = <b>5500</b>；Y = <b>1490</b>。',
    news_title:'新闻',
    news_body:'我们正在开发新的工具、更多语言和集成。更新会发布在这里。',
    news_helipep_title:'新功能：HeliPep',
    news_helipep_body:'推出 HeliPep —— 快速生成 Helical Wheel Plot，支持颜色/形状自定义，并可将多条序列打包为 .zip 导出。',
    about_title:'关于我们', about_sub:'关于 Troplia',
    about_body:'在 Marcelo Henrique Soller Ramada 教授带领下，隶属巴西利亚天主教大学（UCB）相关项目，获 CAPES、CNPq、FAPDF 支持。',
    about_maint_title:'网站维护',
    about_maint_body:'负责团队：<br><br>Gabriel Iudy Yamaguchi Rocha<br>Marcelo Henrique Soller Ramada<br>Rivaldo Varejão Pasqual Saraiva',
    contact_title:'联系', contact_intro:'欢迎联系合作与咨询：',
    contact_list:'Gabriel Iudy Yamaguchi Rocha: <a class="link" href="mailto:yamaguchiiudy@gmail.com">yamaguchiiudy@gmail.com</a><br>Marcelo Henrique Soller Ramada: <a class="link" href="mailto:marceloramada@gmail.com">marceloramada@gmail.com</a><br>Rivaldo Varejão Pasqual Saraiva: <a class="link" href="mailto:rivaldo.vps8@gmail.com">rivaldo.vps8@gmail.com</a>',
    helipep_title:'HeliPep — 肽类 Helical Wheel Plot',
    helipep_subtitle:'可基于单条序列或包含多条序列的 .txt 文件生成 Helical Wheel Plot。',
    helipep_single_title:'单条序列',
    helipep_seq_label:'序列（仅字母）：',
    helipep_name_label:'肽名称（可选；显示在中心）：',
    helipep_preview_btn:'预览',
    helipep_download_btn:'下载 PNG',
    helipep_bulk_title:'多条序列（.txt）',
    helipep_bulk_hint:'示例 .txt 文件（每行一条序列）：',
    helipep_bulk_btn:'生成 .zip',
    helipep_custom_title:'自定义',
    helipep_color_hydrophobic:'疏水性颜色',
    helipep_color_basic:'碱性颜色',
    helipep_color_acidic:'酸性颜色',
    helipep_color_polar:'极性颜色',
    helipep_shape_hydrophobic:'疏水性形状',
    helipep_shape_basic:'碱性形状',
    helipep_shape_acidic:'酸性形状',
    helipep_shape_polar:'极性形状',
    helipep_advanced_title:'高级选项',
    helipep_refs_title:'参考文献'
  },

  ru:{
    tab_home:'Главная', tab_hist:'История', tab_help:'Инструкции', tab_news:'Новости', tab_about:'О нас', tab_contact:'Контакты',
    home_title:'Добро пожаловать в Troplia',
    home_intro:'Платформа Troplia для лабораторных инструментов. Сайт будет дополняться новыми сервисами; сейчас доступен <b>QuantyPep</b> для расчёта концентрации пептидов по абсорбции — просто, быстро и точно.',
    tools_title:'Инструменты',
    tool_calc_title:'QuantyPep – Калькулятор концентрации пептидов (абсорбция)',
    tool_calc_desc:'Расчёт по A280 (С W или С YYY) либо по A205/215/225 (Без W). Результат всегда в микромолях (µM).',
    tool_helipep_title:'HeliPep — Helical Wheel Plot для пептидов',
    tool_helipep_desc:'Мгновенные графики (предпросмотр + загрузка), пакет .txt и настройка цветов/форм.',
    calc_title:'QuantyPep — Калькулятор концентрации пептидов (абсорбция)',
    proj_info:'Информация о проекте', lbl_project:'Название проекта', lbl_date:'Дата', lbl_peptide_name:'Название пептида', lbl_sequence:'Последовательность пептида', opt_field:'(необязательно)',
    ph_project:'Мой эксперимент', ph_peptide_name:'Напр.: Pep-1', ph_sequence:'Напр.: ACDEFGHIKLMNPQRSTVWY',
    calc_type:'Тип расчёта', btn_compute:'Рассчитать', btn_clear:'Сброс', btn_export_csv:'Экспорт CSV', btn_save_txt:'Сохранить TXT', result:'Результат',
    calc_opt_auto:'Авто', calc_opt_w:'С W (A280)', calc_opt_yyy:'С YYY (A280)', calc_opt_now:'Без W (205/215/225)',
    calc_name_w:'С W', calc_name_yyy:'С YYY', calc_name_now:'Без W',
    tip_w:'A280 (W+Y): ε_total = (#W×5500 + #Y×1490).',
    tip_yyy:'A280 (W+Y): ε_total = (#W×5500 + #Y×1490).',
    tip_now:'Используются 205/215/225 нм.',
    obs_w:'Применён режим «С W» (A280, W+Y).', obs_yyy:'Применён режим «С YYY» (A280, W+Y).', obs_now:'Применён режим «Без W» (205/215/225).',
    hist_title:'История', th_project:'Проект', th_date:'Дата (YYYY-MM-DD)', th_peptide:'Пептид', th_type:'Тип', th_dil:'Разведение', th_obs:'Примечание',
    btn_clear_history:'Очистить историю',
    help_title:'Инструкции',
    help_step:'<b>Шаги:</b> укажите <i>Проект</i>, <i>Дату</i> (YYYY-MM-DD), <i>Название пептида</i> и, при желании, <i>Последовательность</i>.',
    help_auto:'<b>Авто:</b> если <b>W ≥ 1</b> → <i>С W</i>; если <b>W = 0</b> и <b>Y ≥ 3</b> → <i>С YYY</i>; иначе → <i>Без W</i>. Режимы W/YYY используют один A280 калькулятор.',
    help_w:'<b>С W / С YYY (A280):</b> <code>ε_total = #W×5500 + #Y×1490</code>.',
    help_yyy:'<b>Примечание:</b> ≥3 тирозина без W — та же формула A280; при наличии W и Y вносят вклад оба.',
    help_now:'<b>Без W (205/215/225):</b> <code>X=(A215−A225)×144</code>; <code>Y=A205×31</code>; <code>mg/mL=((X+Y)/2)×Разведение</code>; <code>µM=(mg/mL ÷ ММ)×1000</code>.',
    help_res:'<b>Результат:</b> всегда в <b>µM</b>.',
    help_mm:'<b>Примечание (ММ):</b> Молекулярная масса рассчитывается автоматически, но может быть указана вручную.',
    help_eps:'<b>Коэффициенты (280 нм):</b> W = <b>5500</b>; Y = <b>1490</b>.',
    news_title:'Новости',
    news_body:'Мы работаем над новыми инструментами, языками и интеграциями. Обновления будут публиковаться здесь.',
    news_helipep_title:'Новое: HeliPep',
    news_helipep_body:'Представляем HeliPep — быстрое создание графиков, настройка цветов/форм и экспорт нескольких последовательностей в .zip.',
    about_title:'О нас', about_sub:'О Troplia',
    about_body:'Под руководством проф. Марсело Энрике Соллер Рамада (UCB). Поддержка CAPES, CNPq, FAPDF; инфраструктура UCB.',
    about_maint_title:'Поддержка сайта',
    about_maint_body:'Команда:<br><br>Gabriel Iudy Yamaguchi Rocha<br>Marcelo Henrique Soller Ramada<br>Rivaldo Varejão Pasqual Saraiva',
    contact_title:'Контакты', contact_intro:'Свяжитесь с нами для информации и сотрудничества:',
    contact_list:'Gabriel Iudy Yamaguchi Rocha: <a class="link" href="mailto:yamaguchiiudy@gmail.com">yamaguchiiudy@gmail.com</a><br>Marcelo Henrique Soller Ramada: <a class="link" href="mailto:marceloramada@gmail.com">marceloramada@gmail.com</a><br>Rivaldo Varejão Pasqual Saraiva: <a class="link" href="mailto:rivaldo.vps8@gmail.com">rivaldo.vps8@gmail.com</a>',
    helipep_title:'HeliPep — Helical Wheel Plot для пептидов',
    helipep_subtitle:'Создавайте Helical Wheel Plot из одной последовательности или из .txt-файла с несколькими последовательностями.',
    helipep_single_title:'Одна последовательность',
    helipep_seq_label:'Последовательность (только буквы):',
    helipep_name_label:'Название пептида (необязательно; по центру):',
    helipep_preview_btn:'Предпросмотр',
    helipep_download_btn:'Скачать PNG',
    helipep_bulk_title:'Несколько последовательностей (.txt)',
    helipep_bulk_hint:'Пример .txt-файла (по одной в строке):',
    helipep_bulk_btn:'Сформировать .zip',
    helipep_custom_title:'Настройка',
    helipep_color_hydrophobic:'Цвет гидрофобных',
    helipep_color_basic:'Цвет основных',
    helipep_color_acidic:'Цвет кислых',
    helipep_color_polar:'Цвет полярных',
    helipep_shape_hydrophobic:'Форма гидрофобных',
    helipep_shape_basic:'Форма основных',
    helipep_shape_acidic:'Форма кислых',
    helipep_shape_polar:'Форма полярных',
    helipep_advanced_title:'Доп. параметры',
    helipep_refs_title:'Ссылки'
  }
};

const SUPPORTED = ['pt','en','es','fr','de','it','ja','zh','ru'];

function setHtmlLang(code){ document.documentElement.setAttribute('lang', code); }
function t(key){ const d = I18N[CURRENT_LANG] || I18N.en; return d[key] ?? key; }

function updatePlaceholders(){
  const proj = document.getElementById('project');
  const pepn = document.getElementById('peptide_name');
  const seq  = document.getElementById('sequence');
  if(proj) proj.placeholder = t('ph_project');
  if(pepn) pepn.placeholder = t('ph_peptide_name');
  if(seq)  seq.placeholder  = t('ph_sequence');
}

function setLang(code){
  CURRENT_LANG = I18N[code] ? code : 'en';
  localStorage.setItem('troplia_lang', CURRENT_LANG);
  setHtmlLang(CURRENT_LANG);

  document.querySelectorAll('[data-i18n]').forEach(el=>{
    const k=el.getAttribute('data-i18n');
    const v=(I18N[CURRENT_LANG]||{})[k];
    if(v!==undefined){
      if(typeof v==='string' && v.includes('<')) el.innerHTML=v; else el.textContent=v;
    }
  });

  const map = {'Auto':'calc_opt_auto','Com W':'calc_opt_w','Com YYY':'calc_opt_yyy','Sem W':'calc_opt_now'};
  const sel = document.getElementById('calc_type');
  if(sel){ [...sel.options].forEach(opt=>{ const key = map[opt.value]; if(key && t(key)) opt.textContent = t(key); }); }

  updatePlaceholders();
  updateCalcTipOnly();
  updateExportLinks();
  if(!document.getElementById('pane-hist').classList.contains('hidden')) updateHistory();
}

function chooseInitialLang(){
  const saved = localStorage.getItem('troplia_lang');
  if (saved && SUPPORTED.includes(saved)) return saved;
  const langs = (navigator.languages && navigator.languages.length ? navigator.languages : [navigator.language || 'en']).map(s => (s||'').toLowerCase());
  for (const l of langs){ const base = l.split('-')[0]; if (SUPPORTED.includes(base)) return base; }
  return 'en';
}

function updateExportLinks(){
  const lang = CURRENT_LANG;
  ['export_txt_main','export_csv_main','export_txt_hist','export_csv_hist'].forEach(id=>{
    const el = document.getElementById(id);
    if(!el) return;
    const base = el.getAttribute('href').split('?')[0];
    el.setAttribute('href', `${base}?lang=${encodeURIComponent(lang)}`);
  });
}

/* ===== Navegação ===== */
function switchTab(tgt){
  const panes=['home','calc','hist','help','news','about','contact'];
  panes.forEach(x=>{
    const el = document.getElementById('pane-'+x);
    const on = (x===tgt);
    if(el){ el.classList.toggle('hidden', !on); el.classList.toggle('block', on); el.setAttribute('aria-hidden', String(!on)); }
    const tb = document.getElementById('tab-'+x);
    if(tb){ tb.classList.toggle('tab-active', on); }
  });
  if(tgt==='hist') updateHistory();
  const ch = document.getElementById('btnClearHistory');
  if(ch && !ch._bound){ ch.addEventListener('click', clearHistory); ch._bound=true; }
}

/* ===== Calculadora ===== */
const fieldsDiv = document.getElementById('fields');
const calcSel   = document.getElementById('calc_type');
const seqInput  = document.getElementById('sequence');
const dateInput = document.getElementById('date');
const calcTip   = document.getElementById('calcTip');

if(dateInput){ dateInput.value = new Date().toISOString().slice(0,10); }

function sanitizeSeq(s){ return (s||'').toUpperCase().replace(/\s+/g,'').replace(/[^ACDEFGHIKLMNPQRSTVWY]/g,''); }
function countLetters(s,ch){ return (s.match(new RegExp(ch,'g'))||[]).length; }

function updateCalcTipOnly(){
  if(!calcSel || !calcTip) return;
  const mode = calcSel.value;
  if(mode==='Com W' || mode==='Com YYY') calcTip.textContent = t('tip_w');
  else if(mode==='Sem W') calcTip.textContent = t('tip_now');
  else calcTip.textContent = '';
}

function makeNumberInput(id,label,defv=''){
  const w = document.createElement('div');
  w.innerHTML = `<label class="block text-sm mb-1">${label}</label><input id="${id}" class="input w-full rounded px-3 py-2" value="${defv}" inputmode="decimal" placeholder="0.0"/>`;
  return w;
}
function makeIntInput(id,label,defv=''){
  const w = document.createElement('div');
  w.innerHTML = `<label class="block text-sm mb-1">${label}</label><input id="${id}" class="input w-full rounded px-3 py-2" value="${defv}" inputmode="numeric" placeholder="0"/>`;
  return w;
}

/* ===== Controle Auto x Manual ===== */
let userTouchedMode = false;
let programmaticSet = false;

function setMode(val){
  programmaticSet = true;
  if(calcSel){ calcSel.value = val; }
  renderFields(val);
  programmaticSet = false;
}

function attachQtyAutoSwitch(){
  const iw=document.getElementById('qty_w');
  const iy=document.getElementById('qty_y');
  const dil=document.getElementById('dilution');
  if(!iw && !iy) return;

  function onChange(){
    if (userTouchedMode) return; // auto-sugestão só quando não tocou no modo

    const w = parseInt(iw?.value||'0',10) || 0;
    const y = parseInt(iy?.value||'0',10) || 0;
    const prevDil = dil ? dil.value : '';

    if(w >= 1){ setMode('Com W'); }
    else if(w === 0 && y >= 3){ setMode('Com YYY'); }
    else { setMode('Sem W'); }

    const d2 = document.getElementById('dilution');
    if(d2 && prevDil) d2.value = prevDil;
  }

  if(iw) iw.addEventListener('input', onChange);
  if(iy) iy.addEventListener('input', onChange);
}

function renderFields(mode){
  if(!fieldsDiv) return;
  fieldsDiv.innerHTML = '';
  if(mode==='Com W' || mode==='Com YYY'){
    fieldsDiv.append(
      makeNumberInput('abs280', 'A280'),
      makeNumberInput('dilution', (I18N[CURRENT_LANG]||I18N.en)['lbl_dilution'] || 'Dilution','1'),
      makeIntInput('qty_w', '#W','0'),
      makeIntInput('qty_y', '#Y','0')
    );
    const s=sanitizeSeq(seqInput?.value||'');
    if(s){
      const w=countLetters(s,'W'), y=countLetters(s,'Y');
      const iw=document.getElementById('qty_w'); const iy=document.getElementById('qty_y');
      if(iw && (iw.value==='' || iw.value==='0') && w>0) iw.value=String(w);
      if(iy && (iy.value==='' || iy.value==='0') && y>0) iy.value=String(y);
    }
    attachQtyAutoSwitch();

  } else if(mode==='Sem W'){
    fieldsDiv.append(
      makeNumberInput('abs205','A205'),
      makeNumberInput('abs215','A215'),
      makeNumberInput('abs225','A225'),
      makeNumberInput('dilution', (I18N[CURRENT_LANG]||I18N.en)['lbl_dilution'] || 'Dilution','1'),
      makeNumberInput('mw','MM (Da)')
    );
    autofillMW();
  }
  updateCalcTipOnly();
}

function autofillMW(){
  const s = sanitizeSeq(seqInput?.value||'');
  const mwEl = document.getElementById('mw');
  if(!mwEl) return;
  if(!s){ mwEl.placeholder='0'; return; }
  const weights = {'A':71.037114,'R':156.101111,'N':114.042927,'D':115.026943,'C':103.009185,'E':129.042593,'Q':128.058578,'G':57.021464,'H':137.058912,'I':113.084064,'L':113.084064,'K':128.094963,'M':131.040485,'F':147.068414,'P':97.052764,'S':87.032028,'T':101.047678,'V':99.068414,'W':186.079313,'Y':163.063329};
  let total = 18.01056;
  for(const ch of s){ total += (weights[ch]||0); }
  mwEl.value = total.toFixed(3);
}

/* ===== Sugerir pelo texto da sequência ===== */
function suggestModeFromSequence(){
  const s = sanitizeSeq(seqInput?.value||'');
  const w = countLetters(s,'W');
  const y = countLetters(s,'Y');
  let suggested = 'Sem W';
  if(w >= 1) suggested = 'Com W';
  else if(y >= 3) suggested = 'Com YYY';
  setMode(suggested);
}

/* ===== Escutas ===== */
if(calcSel){
  calcSel.addEventListener('change', ()=>{
    if(programmaticSet) return;
    userTouchedMode = true;   // escolheu manualmente
    renderFields(calcSel.value);
  });
}

const btnSuggest = document.getElementById('btnSuggest');
if(btnSuggest){
  btnSuggest.addEventListener('click', ()=>{
    userTouchedMode = false;  // volta a permitir auto-sugestão
    suggestModeFromSequence();
  });
}

if(seqInput){
  seqInput.addEventListener('input', ()=>{
    const s = sanitizeSeq(seqInput.value);
    const cleaned = s;
    if(seqInput.value !== cleaned){
      const pos = seqInput.selectionStart;
      seqInput.value = cleaned;
      if(typeof pos === 'number') seqInput.setSelectionRange(Math.max(0,(pos-1)),Math.max(0,(pos-1)));
    }
    if(!userTouchedMode){
      suggestModeFromSequence();
    }else{
      if((calcSel&&calcSel.value)==='Com W' || (calcSel&&calcSel.value)==='Com YYY'){
        const iw=document.getElementById('qty_w'); const iy=document.getElementById('qty_y');
        const wseq = countLetters(cleaned,'W'), yseq = countLetters(cleaned,'Y');
        if(iw && (iw.value==='' || iw.value==='0') && wseq>0) iw.value=String(wseq);
        if(iy && (iy.value==='' || iy.value==='0') && yseq>0) iy.value=String(yseq);
      }else if((calcSel&&calcSel.value)==='Sem W'){
        autofillMW();
      }
    }
  });
}

/* ===== Ações ===== */
const tabsNav = document.querySelector('nav.tabs');
if(tabsNav){
  tabsNav.addEventListener('click', (ev)=>{
    const btn = ev.target.closest('button[data-tab]');
    if(!btn) return;
    switchTab(btn.dataset.tab);
  });
}

const goCalc = document.getElementById('go-calc');
if(goCalc){
  goCalc.addEventListener('click', (e)=>{ e.preventDefault(); switchTab('calc'); });
}

function getPayload(){
  const project = (document.getElementById('project')?.value) || 'Sem Nome';
  const date    = (document.getElementById('date')?.value) || new Date().toISOString().slice(0,10);
  const pepname = (document.getElementById('peptide_name')?.value) || 'Desconhecido';
  const seq     = (document.getElementById('sequence')?.value || '').toUpperCase().replace(/\s+/g,'').replace(/[^ACDEFGHIKLMNPQRSTVWY]/g,'');

  let calcType = (calcSel?.value) || 'Auto';
  if(calcType==='Auto'){
    const w = (seq.match(/W/g)||[]).length;
    const y = (seq.match(/Y/g)||[]).length;
    if(w>=1) calcType='Com W';
    else if(y>=3) calcType='Com YYY';
    else calcType='Sem W';
    setMode(calcType);
  }

  const payload = { project, date, peptide_name: pepname, sequence: seq, calc_type: calcType };

  if(calcType==='Com W' || calcType==='Com YYY'){
    payload.abs280   = document.getElementById('abs280')?.value || '0';
    payload.dilution = document.getElementById('dilution')?.value || '1';
    payload.qty_w    = document.getElementById('qty_w')?.value || '0';
    payload.qty_y    = document.getElementById('qty_y')?.value || '0';
  } else {
    payload.abs205  = document.getElementById('abs205')?.value || '0';
    payload.abs215  = document.getElementById('abs215')?.value || '0';
    payload.abs225  = document.getElementById('abs225')?.value || '0';
    payload.dilution= document.getElementById('dilution')?.value || '1';
    payload.mw      = document.getElementById('mw')?.value || '';
  }
  return payload;
}

function renderResult(entry){
  const box = document.getElementById('results');
  if(!box){ return; }
  if(!entry){ box.textContent=''; return; }
  if(entry.error){
    box.innerHTML = `<span class="warn">${entry.error}</span>`;
    return;
  }
  const label = t('label_conc') || 'Concentração';
  const unit  = t('unit_micromolar') || 'µM';
  const obs = (function(){
    if(entry.calcType==='Com W') return t('obs_w');
    if(entry.calcType==='Com YYY') return t('obs_yyy');
    return t('obs_now');
  })();

  const line = `${label} = ${Number(entry.result_uM).toFixed(2)} ${unit}.`;
  box.innerHTML = line + `\n${obs}`;
}

async function computeNow(){
  const payload = getPayload();

  if(payload.calc_type==='Com W' || payload.calc_type==='Com YYY'){
    const w = parseInt(document.getElementById('qty_w')?.value||'0',10) || 0;
    const y = parseInt(document.getElementById('qty_y')?.value||'0',10) || 0;

    if(payload.calc_type==='Com W' && w<1){
      const msg = (CURRENT_LANG==='en') ? 'For "With W", please set #W ≥ 1 (sequence optional).' :
                  (CURRENT_LANG==='es' ? 'Para "Con W", establezca #W ≥ 1 (secuencia opcional).' :
                  (CURRENT_LANG==='fr' ? 'Pour « Avec W », définissez #W ≥ 1 (séquence facultative).' :
                  (CURRENT_LANG==='de' ? 'Für „Mit W“ bitte #W ≥ 1 angeben (Sequenz optional).' :
                  (CURRENT_LANG==='it' ? 'Per "Con W", imposta #W ≥ 1 (sequenza opzionale).' :
                  (CURRENT_LANG==='ja' ? '「W あり」を使用するには #W ≥ 1 を設定してください（配列は任意）。' :
                  (CURRENT_LANG==='zh' ? '若选择“含 W”，请将 #W 设为 ≥ 1（序列可选）。' :
                  'Para "Com W", defina #W ≥ 1 (sequência opcional).'))))));
      renderResult({error: msg});
      return;
    }
    if(payload.calc_type==='Com YYY' && !(w===0 && y>=3)){
      const msg = (CURRENT_LANG==='en') ? 'Use "With YYY" when #W=0 and #Y≥3 (sequence optional).' :
                  (CURRENT_LANG==='es' ? 'Use "Con YYY" cuando #W=0 y #Y≥3 (secuencia opcional).' :
                  (CURRENT_LANG==='fr' ? 'Utilisez « Avec YYY » quand #W=0 et #Y≥3 (séquence facultative).' :
                  (CURRENT_LANG==='de' ? '„Mit YYY“ verwenden, wenn #W=0 und #Y≥3 (Sequenz optional).' :
                  (CURRENT_LANG==='it' ? 'Usa "Con YYY" quando #W=0 e #Y≥3 (sequenza opzionale).' :
                  (CURRENT_LANG==='ja' ? '「YYY あり」は #W=0 かつ #Y≥3 の場合に使用します（配列は任意）。' :
                  (CURRENT_LANG==='zh' ? '当 #W=0 且 #Y≥3 时请选择“含 YYY”（序列可选）。' :
                  'Use "Com YYY" quando #W=0 e #Y≥3 (sequência opcional).'))))));
      renderResult({error: msg});
      return;
    }
  }

  const resp = await fetch('/api/compute', { method:'POST', headers:{'Content-Type':'application/json'}, body: JSON.stringify(payload) });
  const data = await resp.json();
  renderResult(data);
  if(!data.error) updateHistory();
}

const btnCompute = document.getElementById('btnCompute');
if(btnCompute){ btnCompute.addEventListener('click', computeNow); }

const btnClear = document.getElementById('btnClear');
if(btnClear){
  btnClear.addEventListener('click', ()=>{
    ['abs205','abs215','abs225','abs280','dilution','mw','qty_w','qty_y'].forEach(id=>{ const el=document.getElementById(id); if(el) el.value=''; });
    const resEl = document.getElementById('results');
    if(resEl) resEl.textContent='';
  });
}

/* ===== Histórico ===== */
async function clearHistory(){
  try{
    await fetch('/api/history/clear', {method:'POST'});
  }catch(e){}
  await updateHistory();
}

async function updateHistory(){
  const tBody = document.getElementById('history');
  if(!tBody) return;
  const res = await fetch('/api/history');
  const data = await res.json();
  const items = data.items || [];
  tBody.innerHTML = '';
  for(const e of items){
    const tr = document.createElement('tr');
    tr.innerHTML = `
      <td class="p-2">${e.project||''}</td>
      <td class="p-2">${e.date||''}</td>
      <td class="p-2">${e.peptide||''}</td>
      <td class="p-2">${(e.calcType==='Com W')?t('calc_name_w'):(e.calcType==='Com YYY')?t('calc_name_yyy'):t('calc_name_now')}</td>
      <td class="p-2">${e.a205??''}</td>
      <td class="p-2">${e.a215??''}</td>
      <td class="p-2">${e.a225??''}</td>
      <td class="p-2">${e.abs280??''}</td>
      <td class="p-2">${e.dilution??''}</td>
      <td class="p-2">${e.qty_y??''}</td>
      <td class="p-2">${e.qty_w??''}</td>
      <td class="p-2">${e.mw??''}</td>
      <td class="p-2">${e.result_uM??''}</td>
      <td class="p-2">${(e.calcType==='Com W')?t('obs_w'):(e.calcType==='Com YYY')?t('obs_yyy'):t('obs_now')}</td>
    `;
    tBody.appendChild(tr);
  }
}

/* ===== Inicialização ===== */
function init(){
  const initLang = (function(){
    const saved = localStorage.getItem('troplia_lang');
    if (saved && SUPPORTED.includes(saved)) return saved;
    const langs = (navigator.languages && navigator.languages.length ? navigator.languages : [navigator.language || 'en']).map(s => (s||'').toLowerCase());
    for (const l of langs){ const base = l.split('-')[0]; if (SUPPORTED.includes(base)) return base; }
    return 'en';
  })();
  const sel = document.getElementById('lang');
  if(sel){
    sel.value = initLang;
    sel.addEventListener('change', ()=> setLang(sel.value));
  }
  setLang(initLang);

  switchTab('home');
  suggestModeFromSequence();
  updatePlaceholders();
}
document.addEventListener('DOMContentLoaded', init);
</script>
<script src="/static/helipep.js"></script>
<script>

// ==== QuantyPep Relocate - v6 (destaque de botões internos + highlight de abas do topo) ====
(function(){
  if (window.__tropliaQPRelocateV6) return;
  window.__tropliaQPRelocateV6 = true;

  function $(s){ return document.querySelector(s); }
  function $all(s){ return Array.prototype.slice.call(document.querySelectorAll(s)); }
  function on(el, ev, fn){ if (el) el.addEventListener(ev, fn, false); }

  // ---- Realocar as sections originais (sem alterar conteúdo) quando a QuantyPep abrir ----
  var moved = false;
  function relocateOnce(){
    if (moved) return;
    var host = $('#qp-internal-after');
    if (!host){ return; }
    var hist = $('#pane-hist');
    var help = $('#pane-help');

    if (hist){
      host.appendChild(hist);
      hist.style.display = 'none';
      hist.classList.remove('hidden');
      hist.removeAttribute('aria-hidden');
    }
    if (help){
      host.appendChild(help);
      help.style.display = 'none';
      help.classList.remove('hidden');
      help.removeAttribute('aria-hidden');
    }
    moved = true;
    // binds locais para abrir/fechar
    var btnHist = $('#qp-open-hist');
    var btnHelp = $('#qp-open-help');
    on(btnHist, 'click', function(){
      if (!hist) return;
      var vis = (hist.style.display !== 'none');
      hist.style.display = vis ? 'none' : 'block';
      if (!vis){ document.dispatchEvent(new CustomEvent('troplia:showHistory')); }
    });
    on(btnHelp, 'click', function(){
      if (!help) return;
      var vis = (help.style.display !== 'none');
      help.style.display = vis ? 'none' : 'block';
    });
  }

  // Detecta quando a aba calc abre sem interferir no sistema existente
  function watchCalcVisibility(){
    var calc = $('#pane-calc');
    if (!calc){ relocateOnce(); return; }
    if (!calc.classList.contains('hidden')){
      relocateOnce();
      return;
    }
    var obs = new MutationObserver(function(muts){
      for (var i=0;i<muts.length;i++){
        if (muts[i].attributeName === 'class'){
          var nowVisible = !calc.classList.contains('hidden');
          if (nowVisible){
            relocateOnce();
            obs.disconnect();
            break;
          }
        }
      }
    });
    obs.observe(calc, { attributes: true, attributeFilter: ['class'] });
  }

  // ---- Destaque da aba ativa no topo ----
  function setupTopTabHighlight(){
    var tabs = $all('[data-tab]');
    function applyActiveByPane(){
      // tenta achar qual pane está visível e ativar o tab correspondente
      var visible = null;
      var panes = $all('section[id^="pane-"]');
      for (var i=0;i<panes.length;i++){
        if (!panes[i].classList.contains('hidden')){
          visible = panes[i].id.replace('pane-','');
          break;
        }
      }
      if (!visible) return;
      tabs.forEach(function(b){ b.classList.remove('tab-active'); });
      var btn = document.querySelector('[data-tab="'+visible+'"]');
      if (btn){ btn.classList.add('tab-active'); }
    }

    // 1) clique normal
    tabs.forEach(function(b){
      on(b, 'click', function(){
        // deixa o código original lidar com navegação, só marcamos visualmente
        tabs.forEach(function(o){ o.classList.remove('tab-active'); });
        b.classList.add('tab-active');
        // fallback: em caso de navegação assíncrona, revalidar após um tick
        setTimeout(applyActiveByPane, 50);
      });
    });

    // 2) quando scripts de abas internos trocarem de aba, refletir highlight
    document.addEventListener('troplia:tabActivated', function(){
      setTimeout(applyActiveByPane, 0);
    }, false);

    // 3) estado inicial
    setTimeout(applyActiveByPane, 0);
  }

  if (document.readyState === 'loading'){
    document.addEventListener('DOMContentLoaded', function(){
      watchCalcVisibility();
      setupTopTabHighlight();
    }, false);
  } else {
    watchCalcVisibility();
    setupTopTabHighlight();
  }
})();

</script>

<!-- === Active Tab Highlighter (border-only, no interference) === -->
<script>
(function(){
  if (window.__tabHighlighterV2) return;
  window.__tabHighlighterV2 = true;

  function refreshActive(){
    var panes = document.querySelectorAll('section[id^="pane-"]');
    var active = null;
    panes.forEach(function(p){
      if (!p.classList.contains('hidden')){
        active = p.id.replace('pane-','');
      }
    });
    var tabs = document.querySelectorAll('nav.tabs button[data-tab]');
    tabs.forEach(function(btn){
      if (btn.getAttribute('data-tab') === active){
        btn.classList.add('tab-active');
        btn.setAttribute('aria-selected','true');
      } else {
        btn.classList.remove('tab-active');
        btn.removeAttribute('aria-selected');
      }
    });
  }

  function bindClicks(){
    var tabs = document.querySelectorAll('nav.tabs button[data-tab]');
    tabs.forEach(function(btn){
      btn.addEventListener('click', function(){
        // Wait a tick so the site's native tab logic runs first
        setTimeout(refreshActive, 40);
      }, false);
    });
  }

  // Update highlight when other scripts switch tabs
  document.addEventListener('troplia:tabActivated', function(){
    setTimeout(refreshActive, 0);
  }, false);

  // Observe visibility changes to panes (class 'hidden') to keep highlight in sync
  var mo = new MutationObserver(function(muts){
    for (var i=0;i<muts.length;i++){
      if (muts[i].type === 'attributes' && muts[i].attributeName === 'class'){
        refreshActive();
        break;
      }
    }
  });
  function observePanes(){
    document.querySelectorAll('section[id^="pane-"]').forEach(function(p){
      mo.observe(p, { attributes:true, attributeFilter:['class'] });
    });
  }

  if (document.readyState === 'loading'){
    document.addEventListener('DOMContentLoaded', function(){
      bindClicks(); observePanes(); refreshActive();
    }, false);
  } else {
    bindClicks(); observePanes(); refreshActive();
  }
})();
</script>

</body>
</html>
"""

@app.get("/", response_class=HTMLResponse)
async def index(request: Request):
    resp = HTMLResponse(INDEX_HTML)
    sid = request.cookies.get('troplia_id')
    if not sid:
        sid = secrets.token_urlsafe(16)
        resp.set_cookie('troplia_id', sid, max_age=60*60*24*365*2, httponly=True, samesite='Lax')
    return resp

# ------------------ ROTA PARA ADS.TXT ------------------
from fastapi.responses import FileResponse
from pathlib import Path

@app.get("/ads.txt")
async def serve_ads_txt():
    file_path = Path(__file__).with_name("ads.txt")
    return FileResponse(file_path, media_type="text/plain; charset=utf-8")

# ================== API ==================
def to_float(x, default=0.0):
    try:
        return float(str(x).replace(',','.'))
    except Exception:
        return default

def to_int(x, default=0):
    try:
        return int(x)
    except Exception:
        return default

@app.post("/api/compute/")
async def api_compute(payload: dict, request: Request):
    project = (payload.get("project") or "Sem Nome").strip()
    date_str = payload.get("date") or datetime.today().strftime('%Y-%m-%d')  # YYYY-MM-DD
    peptide_name = (payload.get("peptide_name") or "Desconhecido").strip()
    sequence = clean_seq(payload.get("sequence") or "")
    raw_type = payload.get("calc_type") or "Auto"

    qty_y = to_int(payload.get("qty_y"), (sequence.count('Y') if sequence else 0))
    qty_w = to_int(payload.get("qty_w"), (sequence.count('W') if sequence else 0))

    sraw = (raw_type or '').strip().lower()
    if sraw == 'auto' or sraw == '':
        calc_type = autodetect_type_counts(qty_w, qty_y) if (qty_w or qty_y) else autodetect_type(sequence)
    else:
        calc_type = canonical_calc_type(raw_type, sequence)

    if calc_type == 'Com W' and qty_w < 1:
        return JSONResponse({"error":"Para 'Com W', defina #W ≥ 1 (sequência opcional)."}, status_code=400)
    if calc_type == 'Com YYY' and not (qty_w == 0 and qty_y >= 3):
        return JSONResponse({"error":"Use 'Com YYY' quando #W=0 e #Y≥3 (sequência opcional)."}, status_code=400)

    result_uM = None
    obs = ""

    try:
        if calc_type in ("Com W", "Com YYY"):
            abs280 = to_float(payload.get("abs280"), 0.0)
            dil    = to_float(payload.get("dilution"), 1.0)
            denom = (qty_w * EXT_W) + (qty_y * EXT_Y)  # ε_total
            if denom <= 0:
                return JSONResponse({"error":"Para A280, informe #W/#Y (>0) ou altere para 'Sem W'."}, status_code=400)
            result_uM = (abs280/denom) * dil * 1_000_000
            obs = "Modo Com W (A280, W+Y) aplicado." if calc_type=="Com W" else "Modo Com YYY (A280, W+Y) aplicado."
        else:  # Sem W
            a205 = to_float(payload.get("abs205"), 0.0)
            a215 = to_float(payload.get("abs215"), 0.0)
            a225 = to_float(payload.get("abs225"), 0.0)
            dil  = to_float(payload.get("dilution"), 1.0)
            mw   = to_float(payload.get("mw"), mw_from_seq(sequence) if sequence else 0.0)
            x = (a215 - a225) * 144
            yv = a205 * 31
            mgml = ((x + yv) / 2) * dil
            if mw <= 0:
                return JSONResponse({"error":"Informe a massa molecular (Da) ou forneça sequência válida."}, status_code=400)
            result_uM = (mgml / mw) * 1000
            obs = "Modo Sem W (205/215/225) aplicado."

        entry = {
            "project": project,
            "date": date_str,
            "peptide": peptide_name,
            "sequence": sequence,
            "calcType": calc_type,
            "abs280": payload.get("abs280"),
            "a205": payload.get("abs205"),
            "a215": payload.get("abs215"),
            "a225": payload.get("abs225"),
            "dilution": payload.get("dilution"),
            "qty_y": qty_y,
            "qty_w": qty_w,
            "mw": payload.get("mw"),
            "result_uM": round(result_uM, 2) if result_uM is not None else None,
            "obs": obs,
        }
        sid = request.cookies.get('troplia_id')
        if not sid:
            sid = secrets.token_urlsafe(16)
        HISTORY_BY_USER.setdefault(sid, []).append(entry)
        return entry

    except Exception:
        return JSONResponse({"error":"Erro no cálculo. Verifique os dados."}, status_code=400)

@app.get("/api/history")
async def api_history(request: Request):
    return {"items": _get_user_history(request)}

def resolve_lang(request: Request) -> str:
    ql = request.query_params.get('lang')
    if ql and ql in SUPPORTED_LANGS:
        return ql
    return pick_lang_from_header(request.headers.get('accept-language'))

def localized_mode(calc_type: str, lang: str) -> str:
    return CALC_NAME_I18N.get(lang, CALC_NAME_I18N['en']).get(calc_type, calc_type)

def localized_obs(calc_type: str, lang: str) -> str:
    return OBS_I18N.get(lang, OBS_I18N['en']).get(calc_type, "")

@app.get("/api/export/csv")
async def api_export_csv(request: Request):
    hist = _get_user_history(request)
    if not hist:
        return JSONResponse({"error":"Sem dados para exportar."}, status_code=400)

    lang = resolve_lang(request)
    headers_disp = make_download_headers(hist[-1].get("project") or "Troplia", "csv")
    header_names = CSV_HEADERS.get(lang, CSV_HEADERS['en'])

    output = StringIO()
    writer = csv.DictWriter(output, fieldnames=header_names)
    writer.writeheader()

    for e in hist:
        row = {
            header_names[0]:  e['project'],
            header_names[1]:  e['date'],
            header_names[2]:  e['peptide'],
            header_names[3]:  e['sequence'],
            header_names[4]:  localized_mode(e['calcType'], lang),
            header_names[5]:  e.get('a205'),
            header_names[6]:  e.get('a215'),
            header_names[7]:  e.get('a225'),
            header_names[8]:  e.get('abs280'),
            header_names[9]:  e.get('dilution'),
            header_names[10]: e.get('qty_y'),
            header_names[11]: e.get('qty_w'),
            header_names[12]: e.get('mw'),
            header_names[13]: e.get('result_uM'),
            header_names[14]: localized_obs(e['calcType'], lang),
        }
        writer.writerow(row)

    output.seek(0)
    return HTMLResponse(output.getvalue(), headers=headers_disp, media_type="text/csv")

@app.get("/api/export/txt")
async def api_export_txt(request: Request):
    hist = _get_user_history(request)
    if not hist:
        return JSONResponse({"error":"Sem dados para exportar."}, status_code=400)

    lang = resolve_lang(request)
    labels = TXT_LABELS.get(lang, TXT_LABELS['en'])
    headers_disp = make_download_headers(hist[-1].get("project") or "Troplia", "txt")

    output = StringIO()
    for e in hist:
        output.write(f"{labels['project']}: {e['project']}\n")
        output.write(f"{labels['date']}: {e['date']}\n")
        output.write(f"{labels['peptide']}: {e['peptide']}\n")
        output.write(f"{labels['sequence']}: {e['sequence']}\n")
        output.write(f"{labels['type']}: {localized_mode(e['calcType'], lang)}\n")
        if e.get('result_uM') is not None:
            output.write(f"{labels['result']}: {e['result_uM']} {labels['unit']}\n")
        obs_line = localized_obs(e['calcType'], lang)
        if obs_line:
            output.write(f"{labels['obs']}: {obs_line}\n")
        output.write("═"*50 + "\n")

    output.seek(0)
    return HTMLResponse(output.getvalue(), headers=headers_disp, media_type="text/plain")


@app.post("/api/history/clear")
async def api_history_clear(request: Request):
    sid = request.cookies.get('troplia_id')
    if sid and sid in HISTORY_BY_USER:
        HISTORY_BY_USER[sid].clear()
    return {"ok": True}
# ========= HeliPep API =========
from fastapi import UploadFile, File, Form
from fastapi.responses import StreamingResponse
from troplia_helipep import bulk_generate_zip, save_wheel_png, png_bytes_for_sequence
from pathlib import Path
import io
from fastapi.responses import PlainTextResponse, StreamingResponse

from fastapi.responses import PlainTextResponse, StreamingResponse

@app.post("/api/helipep/preview")
async def api_helipep_preview(
    sequence: str = Form(...),
    peptide_name: str = Form(""),
    angle_per_residue_deg: float = Form(100.0),
    start_angle_deg: float = Form(0.0),
    clockwise: bool = Form(True),
    layer_every_n_residues: int = Form(18),
    draw_links: bool = Form(True),
    draw_circle_guide: bool = Form(False),         # <- agora False por padrão
    annotate_numbers: bool = Form(True),
    annotate_letters: bool = Form(True),
    show_hydrophobic_moment: bool = Form(False),   # <- agora False por padrão
    color_hydrophobic: str = Form(""),
    color_basic: str = Form(""),
    color_acidic: str = Form(""),
    color_polar: str = Form(""),
    shape_hydrophobic: str = Form(""),
    shape_basic: str = Form(""),
    shape_acidic: str = Form(""),
    shape_polar: str = Form("")
):
    try:
        print(f"[HeliPep] preview len={len(sequence)} name='{peptide_name}'")
        color_cfg = {}
        if color_hydrophobic: color_cfg["hydrophobic"] = color_hydrophobic
        if color_basic:       color_cfg["basic"]       = color_basic
        if color_acidic:      color_cfg["acidic"]      = color_acidic
        if color_polar:       color_cfg["polar"]       = color_polar

        shape_cfg = {}
        if shape_hydrophobic: shape_cfg["hydrophobic"] = shape_hydrophobic
        if shape_basic:       shape_cfg["basic"]       = shape_basic
        if shape_acidic:      shape_cfg["acidic"]      = shape_acidic
        if shape_polar:       shape_cfg["polar"]       = shape_polar

        from troplia_helipep import png_bytes_for_sequence
        png = png_bytes_for_sequence(
            sequence,
            angle_per_residue_deg=angle_per_residue_deg,
            start_angle_deg=start_angle_deg,
            clockwise=clockwise,
            layer_every_n_residues=layer_every_n_residues,
            draw_links=draw_links,
            draw_circle_guide=draw_circle_guide,
            annotate_numbers=annotate_numbers,
            annotate_letters=annotate_letters,
            show_hydrophobic_moment=show_hydrophobic_moment,
            color_cfg=color_cfg,
            shape_cfg=shape_cfg,
            center_name=(peptide_name.strip() or None),
            top_title=None
        )
        import io
        return StreamingResponse(io.BytesIO(png), media_type="image/png")
    except Exception as e:
        msg = f"Erro ao gerar preview: {type(e).__name__}: {e}"
        print("[HeliPep] ", msg)
        return PlainTextResponse(msg, status_code=400)

@app.post("/api/helipep/bulk")
async def api_helipep_bulk(
    file: UploadFile = File(...),
    peptide_name: str = Form(""),  # opcional
    angle_per_residue_deg: float = Form(100.0),
    start_angle_deg: float = Form(0.0),
    clockwise: bool = Form(True),
    layer_every_n_residues: int = Form(18),
    draw_links: bool = Form(True),
    draw_circle_guide: bool = Form(True),
    annotate_numbers: bool = Form(True),
    annotate_letters: bool = Form(True),
    show_hydrophobic_moment: bool = Form(True),
    color_hydrophobic: str = Form(""),
    color_basic: str = Form(""),
    color_acidic: str = Form(""),
    color_polar: str = Form(""),
    shape_hydrophobic: str = Form(""),
    shape_basic: str = Form(""),
    shape_acidic: str = Form(""),
    shape_polar: str = Form("")
):
    text = (await file.read()).decode("utf-8", errors="ignore")
    seq_list = [ln for ln in text.splitlines() if ln.strip()]

    color_cfg = {}
    if color_hydrophobic: color_cfg["hydrophobic"] = color_hydrophobic
    if color_basic: color_cfg["basic"] = color_basic
    if color_acidic: color_cfg["acidic"] = color_acidic
    if color_polar: color_cfg["polar"] = color_polar

    shape_cfg = {}
    if shape_hydrophobic: shape_cfg["hydrophobic"] = shape_hydrophobic
    if shape_basic: shape_cfg["basic"] = shape_basic
    if shape_acidic: shape_cfg["acidic"] = shape_acidic
    if shape_polar: shape_cfg["polar"] = shape_polar

    data = bulk_generate_zip(
        seq_list,
        angle_per_residue_deg=angle_per_residue_deg,
        start_angle_deg=start_angle_deg,
        clockwise=clockwise,
        layer_every_n_residues=layer_every_n_residues,
        draw_links=draw_links,
        draw_circle_guide=draw_circle_guide,
        annotate_numbers=annotate_numbers,
        annotate_letters=annotate_letters,
        show_hydrophobic_moment=show_hydrophobic_moment,
        color_cfg=color_cfg,
        shape_cfg=shape_cfg,
        center_name=(peptide_name.strip() or None),
        top_title=None
    )
    buf = io.BytesIO(data)
    return StreamingResponse(
        buf,
        media_type="application/zip",
        headers={"Content-Disposition": "attachment; filename=helipep_results.zip"}
    )


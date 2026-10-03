// static/helipep.js â€” roteador fechado (1 painel visÃ­vel), nav por IDs e hash; HeliPep em PT-BR
(function () {
  const $ = (s) => document.querySelector(s);

  // ---------- Cores padrÃ£o por classe ----------
  const DEFAULT_COLORS = {
    hydrophobic: "#008000", // verde
    basic: "#FF0000",       // vermelho
    acidic: "#0000FF",      // azul
    polar: "#ADD8E6"        // azul claro
  };

  function ensureDefaultColorPickers() {
    const map = {
      "#color_hydrophobic": DEFAULT_COLORS.hydrophobic,
      "#color_basic": DEFAULT_COLORS.basic,
      "#color_acidic": DEFAULT_COLORS.acidic,
      "#color_polar": DEFAULT_COLORS.polar
    };
    Object.entries(map).forEach(([sel, hex]) => {
      const el = document.querySelector(sel);
      if (el && (!el.value || el.value.toUpperCase() === "#000000")) el.value = hex;
    });
  }

  // ---------- Roteador: mostra SOMENTE 1 painel ----------
  const PANES = [
    "pane-home", "pane-calc", "pane-hist", "pane-help",
    "pane-news", "pane-about", "pane-contact", "pane-helipep"
  ];

  function showPane(idToShow) {
    PANES.forEach(id => {
      const el = document.getElementById(id);
      if (!el) return;
      const isTarget = id === idToShow;
      el.classList.toggle("hidden", !isTarget);
      el.setAttribute("aria-hidden", String(!isTarget));
    });
    window.scrollTo({ top: 0, behavior: "smooth" });
  }

  // Mapeia hash -> pane
  function paneFromHash(hash) {
    const h = (hash || "#home").toLowerCase();
    const ROUTES = {
      "#home": "pane-home", "#inicio": "pane-home",
      "#quantypep": "pane-calc", "#calc": "pane-calc",
      "#hist": "pane-hist", "#historico": "pane-hist",
      "#instrucoes": "pane-help", "#ajuda": "pane-help",
      "#novidades": "pane-news",
      "#sobre": "pane-about", "#sobrenos": "pane-about",
      "#contato": "pane-contact",
      "#helipep": "pane-helipep"
    };
    return ROUTES[h] || "pane-home";
  }

  function routeFromHash() {
    const pane = paneFromHash(location.hash);
    showPane(pane);
    if (pane === "pane-helipep") ensureDefaultColorPickers();
  }

  // ---------- OpÃ§Ãµes do formulÃ¡rio HeliPep ----------
  function gatherOptions() {
    return {
      peptide_name: ($("#helipep-name")?.value || "").trim(),
      angle_per_residue_deg: +($("#opt_angle")?.value || 100),
      start_angle_deg: +($("#opt_start")?.value || 0),
      clockwise: $("#opt_clockwise")?.checked ?? true,
      layer_every_n_residues: +($("#opt_layer_n")?.value || 18),
      draw_links: $("#opt_links")?.checked ?? true,
      draw_circle_guide: $("#opt_guide")?.checked ?? false,   // padrÃ£o: desmarcado
      annotate_numbers: $("#opt_nums")?.checked ?? true,
      annotate_letters: $("#opt_letters")?.checked ?? true,
      show_hydrophobic_moment: $("#opt_moment")?.checked ?? false, // padrÃ£o: desmarcado
      color_hydrophobic: $("#color_hydrophobic")?.value || "",
      color_basic: $("#color_basic")?.value || "",
      color_acidic: $("#color_acidic")?.value || "",
      color_polar: $("#color_polar")?.value || "",
      shape_hydrophobic: $("#shape_hydrophobic")?.value || "",
      shape_basic: $("#shape_basic")?.value || "",
      shape_acidic: $("#shape_acidic")?.value || "",
      shape_polar: $("#shape_polar")?.value || ""
    };
  }

  // ---------- Preview / PNG ----------
  async function preview() {
    const seq = ($("#helipep-seq")?.value || "").trim();
    if (!seq) { alert("Informe a sequÃªncia"); return; }

    const opts = gatherOptions();
    const fd = new FormData();
    fd.append("sequence", seq);
    Object.entries(opts).forEach(([k, v]) => fd.append(k, v));

    try {
      const res = await fetch("/api/helipep/preview", { method: "POST", body: fd });
      if (!res.ok) {
        const txt = await res.text();
        console.error("Preview falhou:", res.status, txt);
        alert("Erro ao gerar preview: " + res.status);
        return;
      }
      const blob = await res.blob();
      const url = URL.createObjectURL(blob);
      const img = $("#helipep-img");
      if (img) img.src = url;
      const a = $("#helipep-preview-dl");
      if (a) {
        a.href = url;
        a.download = (opts.peptide_name ? `${opts.peptide_name}.png` : "helipep_preview.png");
      }
    } catch (err) {
      console.error(err);
      alert("Falha de comunicaÃ§Ã£o com o servidor.");
    }
  }

  // ---------- ZIP em lote ----------
  async function bulk() {
    const fileInput = $("#helipep-file");
    if (!fileInput?.files?.length) { alert("Selecione um arquivo .txt"); return; }

    const opts = gatherOptions();
    const fd = new FormData();
    fd.append("file", fileInput.files[0]);
    Object.entries(opts).forEach(([k, v]) => fd.append(k, v));

    try {
      const res = await fetch("/api/helipep/bulk", { method: "POST", body: fd });
      if (!res.ok) {
        const txt = await res.text();
        console.error("ZIP falhou:", res.status, txt);
        alert("Erro ao gerar ZIP: " + res.status);
        return;
      }
      const blob = await res.blob();
      const url = URL.createObjectURL(blob);
      const a = document.createElement("a");
      a.href = url; a.download = "helipep_results.zip";
      document.body.appendChild(a); a.click(); a.remove();
      URL.revokeObjectURL(url);
    } catch (err) {
      console.error(err);
      alert("Falha de comunicaÃ§Ã£o com o servidor.");
    }
  }

  // ---------- Abrir diretamente o HeliPep ----------
  function openHeliPep(ev) {
    ev?.preventDefault();
    ev?.stopPropagation();
    location.hash = "#helipep";
    showPane("pane-helipep");
    ensureDefaultColorPickers();
  }

  // ---------- Seletor de arquivo em PT-BR ----------
  function setupFileInputPT() {
    const input = document.getElementById("helipep-file");
    const nameEl = document.getElementById("helipep-file-name");
    if (!input || !nameEl) return;
    nameEl.textContent = "Nenhum arquivo selecionado";
    input.addEventListener("change", () => {
      nameEl.textContent = (input.files && input.files.length > 0)
        ? input.files[0].name
        : "Nenhum arquivo selecionado";
    });
  }

  // ---------- Boot ----------
  document.addEventListener("DOMContentLoaded", () => {
    // 1) Roteamento por hash
    routeFromHash();
    window.addEventListener("hashchange", routeFromHash, true);

    // 2) NavegaÃ§Ã£o por IDs do menu (fecha HeliPep ao trocar de aba)
    const NAV_ID_TO_HASH = {
      "tab-home": "#home",
      "tab-hist": "#historico",
      "tab-instr": "#instrucoes",
      "tab-news": "#novidades",
      "tab-about": "#sobre",
      "tab-contact": "#contato"
    };
    Object.entries(NAV_ID_TO_HASH).forEach(([id, hash]) => {
      const el = document.getElementById(id);
      if (!el) return;
      el.addEventListener("click", (ev) => {
        ev.preventDefault(); ev.stopPropagation();
        location.hash = hash; // dispara routeFromHash -> showPane
      }, true);
    });

    // 3) Qualquer <a href="#..."> na barra de navegaÃ§Ã£o tambÃ©m fecha o HeliPep
    document.addEventListener("click", (ev) => {
      const t = ev.target;
      const a = t && t.closest ? t.closest('a[href^="#"]') : null;
      if (!a) return;
      const hash = a.getAttribute("href");
      const pane = paneFromHash(hash);
      if (pane) { ev.preventDefault(); location.hash = hash; } // routeFromHash farÃ¡ o resto
    }, true);

    // 4) AÃ§Ãµes do HeliPep
    $("#helipep-preview")?.addEventListener("click", preview);
    $("#helipep-bulk")?.addEventListener("click", bulk);

    // 5) Atalhos para abrir o HeliPep
    $("#go-helipep")?.addEventListener("click", openHeliPep, true);
    document.addEventListener("click", (ev) => {
      const t = ev.target;
      const hit = t && t.closest ? t.closest('#go-helipep, .go-helipep, [data-open="helipep"], a[href*="helipep"]') : null;
      if (hit) openHeliPep(ev);
    }, true);

    // 6) Defaults
    ensureDefaultColorPickers();
    setupFileInputPT();
  });
})();
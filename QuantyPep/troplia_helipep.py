# troplia_helipep.py — publication-ready (mantém API; adiciona rim-light e refinamentos)
from pathlib import Path
import io
import zipfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import FancyArrowPatch, Circle, Ellipse
import numpy as np

# ---------------- Helpers (cores) ----------------
def _hex_to_rgb(h: str):
    h = h.lstrip("#")
    return tuple(int(h[i:i+2], 16)/255.0 for i in (0, 2, 4))

def _rgb_to_hex(rgb):
    return "#" + "".join(f"{int(max(0,min(1,c))*255):02X}" for c in rgb)

def _mix(color_hex: str, with_hex: str, t: float) -> str:
    a = _hex_to_rgb(color_hex); b = _hex_to_rgb(with_hex)
    m = tuple((1-t)*av + t*bv for av,bv in zip(a,b))
    return _rgb_to_hex(m)

# ---------------- Estilo global ------------------
plt.rcParams.update({
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.transparent": True,
    "axes.linewidth": 0.7,
    "path.simplify": False,
    "text.antialiased": True,
    "lines.antialiased": True,
    "patch.antialiased": True,
    "font.size": 12,
    "font.family": "DejaVu Sans",
})

# Paleta padrão por resíduo (mantida)
AA_COLORS_COMPLETE = {
    'A':'#288f3a','V':'#288f3a','I':'#288f3a','L':'#288f3a','M':'#288f3a','F':'#288f3a','W':'#288f3a','Y':'#288f3a','G':'#288f3a',
    'R':'#e74c3c','K':'#e74c3c','H':'#e74c3c',
    'D':'#2978ff','E':'#2978ff',
    'S':'#8ecae6','T':'#8ecae6','N':'#8ecae6','Q':'#8ecae6','C':'#8ecae6'
}
DEFAULT_COLOR = '#2d2d2d'

CLASS_MAP = {
    "hydrophobic": set(list("AVILMFWYG")),
    "basic": set(list("RKH")),
    "acidic": set(list("DE")),
    "polar": set(list("STNQC")),
}

def _aa_class(aa: str) -> str:
    aa = aa.upper()
    if aa in CLASS_MAP["basic"]: return "basic"
    if aa in CLASS_MAP["acidic"]: return "acidic"
    if aa in CLASS_MAP["polar"]: return "polar"
    if aa in CLASS_MAP["hydrophobic"]: return "hydrophobic"
    return "other"

def _shape_for_class(cls: str, shape_cfg: dict) -> str:
    # 'o' círculo (esfera), 's' quadrado, '^' triângulo, 'D' losango
    return shape_cfg.get(cls, "o")

def _color_for_residue(residue: str, color_cfg: dict) -> str:
    cls = _aa_class(residue)
    if cls in color_cfg:
        return color_cfg[cls]
    return AA_COLORS_COMPLETE.get(residue.upper(), DEFAULT_COLOR)

# ---------------- Plot principal -----------------
def helical_wheel_plot(
    sequence: str,
    angle_per_residue_deg: float = 100.0,
    layer_every_n_residues: int = 18,
    start_angle_deg: float = 0.0,
    clockwise: bool = True,
    draw_links: bool = True,
    draw_circle_guide: bool = False,        # desativado por padrão (mantido)
    color_cfg: dict | None = None,
    shape_cfg: dict | None = None,
    annotate_numbers: bool = True,
    annotate_letters: bool = True,
    show_hydrophobic_moment: bool = False,  # desativado por padrão (mantido)
    center_name: str | None = None,
    top_title: str | None = None
):
    """Desenha o helical wheel (esferas glossy + rim-light; linhas nítidas; μH opcional)."""
    color_cfg = color_cfg or {}
    shape_cfg = shape_cfg or {}

    seq = "".join([c for c in sequence.strip().upper() if c.isalpha()])
    n = len(seq)
    if n == 0:
        raise ValueError("Empty sequence.")

    # Geometria (mantida)
    step = np.deg2rad(angle_per_residue_deg) * (1 if clockwise else -1)
    base = np.deg2rad(start_angle_deg)
    angles = (base + np.arange(n) * step) % (2*np.pi)
    radii  = 1.0 + (np.arange(n) // layer_every_n_residues) * 0.22

    # 0° no topo
    x = radii * np.sin(angles)
    y = radii * np.cos(angles)

    # Figura
    fig, ax = plt.subplots(figsize=(6.2, 6.2))
    ax.set_aspect("equal", adjustable="box")
    lim = max(1.65, 1.05 + (n // layer_every_n_residues)*0.28)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.axis("off")

    # Círculo-guia (opcional)
    if draw_circle_guide:
        circ = plt.Circle((0,0), 1.05, fill=False, ls=(0,(3,3)), lw=0.8, color="#9aa0a6", alpha=0.85)
        ax.add_artist(circ)

    # Star polygon (detalhe estético)
    try:
        k = 2; m = 10; rad = 0.80
        th = np.linspace(0, 2*np.pi, m, endpoint=False) + np.deg2rad(90)
        pts = [(rad*np.cos(t), rad*np.sin(t)) for t in th]
        order = [(i*k) % m for i in range(m+1)]
        xs = [pts[j][0] for j in order]; ys = [pts[j][1] for j in order]
        ax.plot(xs, ys, color="#1f1f1f", lw=1.05, alpha=0.95,
                path_effects=[pe.withStroke(linewidth=2.2, foreground="white", alpha=0.9)], zorder=0)
    except Exception:
        pass

    # Ligações nítidas (stroke branco por baixo)
    if draw_links and n > 1:
        r_dot = 0.085
        for i in range(n-1):
            x0 = x[i] - r_dot*np.sin(angles[i])
            y0 = y[i] - r_dot*np.cos(angles[i])
            x1 = x[i+1] - r_dot*np.sin(angles[i+1])
            y1 = y[i+1] - r_dot*np.cos(angles[i+1])
            ax.plot([x0, x1], [y0, y1],
                    color="#1f1f1f", lw=1.25, solid_capstyle="round",
                    path_effects=[pe.withStroke(linewidth=2.2, foreground="white", alpha=0.95)])

    # Resíduos: esferas com sombra + brilho + rim-light; número ao lado
    for idx, (aa, xi, yi) in enumerate(zip(seq, x, y), start=1):
        cls = _aa_class(aa)
        base_color = _color_for_residue(aa, color_cfg)
        marker = _shape_for_class(cls, shape_cfg)

        r = 0.185                          # raio da esfera
        shadow_off = (0.040, -0.040)       # sombra “abaixo/direita”
        gloss_off  = (-0.055, 0.055)       # highlight “topo/esquerda”
        rim_off    = (0.060, -0.060)       # rim-light “contraluz”
        ring_color = _mix(base_color, "#000000", 0.22)
        face_color = base_color

        if marker == "o":
            ax.add_patch(Circle((xi + shadow_off[0], yi + shadow_off[1]),
                                r*0.90, color="black", alpha=0.17, lw=0, zorder=1))
            ax.add_patch(Circle((xi, yi), r, color=face_color, ec=ring_color, lw=1.2, zorder=2))
            ax.add_patch(Circle((xi + gloss_off[0], yi + gloss_off[1]),
                                r*0.56, color="white", alpha=0.28, lw=0, zorder=3))
            ax.add_patch(Ellipse((xi + rim_off[0], yi + rim_off[1]),
                                 width=r*1.45, height=r*0.55, angle=20,
                                 color="white", alpha=0.10, lw=0, zorder=3))
        else:
            ax.plot([xi], [yi], marker=marker, ms=19, mfc=face_color, mew=1.2, mec=ring_color, zorder=2)

        # Letra no centro com halo branco
        if annotate_letters:
            ax.text(xi, yi, aa, ha="center", va="center",
                    fontsize=11.2, color="#0e0e0e", weight="bold",
                    path_effects=[pe.withStroke(linewidth=2.8, foreground="white", alpha=0.96)],
                    zorder=4)

        # Número “subscrito” ao lado da letra
        if annotate_numbers:
            num_x = xi + r*0.60
            num_y = yi - r*0.60
            ax.text(num_x, num_y, str(idx), ha="center", va="center",
                    fontsize=8.2, color="#2a2a2a",
                    path_effects=[pe.withStroke(linewidth=2.2, foreground="white", alpha=0.95)],
                    zorder=5)

    # Momento hidrofóbico (seta elegante; só se ativado)
    if show_hydrophobic_moment:
        eisenberg = {
            "A":1.8,"R":-4.5,"N":-3.5,"D":-3.5,"C":2.5,"Q":-3.5,"E":-3.5,"G":-0.4,"H":-3.2,"I":4.5,
            "L":3.8,"K":-3.9,"M":1.9,"F":2.8,"P":-1.6,"S":-0.8,"T":-0.7,"W":-0.9,"Y":-1.3,"V":4.2
        }
        vec = np.array([0.0, 0.0])
        for i, aa in enumerate(seq):
            h = eisenberg.get(aa, 0.0)
            theta = angles[i]
            vec += np.array([np.sin(theta), np.cos(theta)]) * h

        mag = float(np.linalg.norm(vec))
        if mag > 1e-6:
            u = vec / mag
            L = min(0.45 + 0.18*np.log1p(mag), 0.95*1.05)
            start = np.array([0.0, 0.0])
            end   = start + u * L
            arr = FancyArrowPatch(
                posA=(start[0], start[1]),
                posB=(end[0], end[1]),
                arrowstyle='-|>',
                mutation_scale=30,
                linewidth=2.9,
                color="#2c7bb6",
                zorder=6
            )
            arr.set_path_effects([pe.withStroke(linewidth=4.6, foreground="white", alpha=0.95)])
            ax.add_patch(arr)
            ax.text(end[0], end[1], "μH", fontsize=11, weight="bold", ha="left", va="bottom",
                    color="#2c7bb6",
                    path_effects=[pe.withStroke(linewidth=3, foreground="white", alpha=0.9)])

    # Nome no centro / título (mantido)
    if center_name:
        ax.text(0, 0, center_name, ha="center", va="center",
                fontsize=12.5, weight="bold",
                path_effects=[pe.withStroke(linewidth=3, foreground="white", alpha=0.9)])
    if top_title:
        ax.set_title(top_title, fontsize=16, weight="bold", pad=10)

    plt.tight_layout()
    return fig, ax

# ---------------- APIs utilitárias (inalteradas) ----------------
def png_bytes_for_sequence(sequence: str, **kwargs) -> bytes:
    buf = io.BytesIO()
    fig, _ = helical_wheel_plot(sequence, **kwargs)
    fig.savefig(buf, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return buf.read()

def save_wheel_png(sequence: str, path: Path, **kwargs):
    fig, _ = helical_wheel_plot(sequence, **kwargs)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)

def bulk_generate_zip(seq_list, **kwargs) -> bytes:
    mem = io.BytesIO()
    with zipfile.ZipFile(mem, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:
        for idx, raw in enumerate(seq_list, start=1):
            seq = "".join([c for c in raw.strip().upper() if c.isalpha()])
            if not seq:
                continue
            buf = io.BytesIO()
            fig, _ = helical_wheel_plot(seq, **kwargs)
            fig.savefig(buf, format="png", dpi=300, bbox_inches="tight")
            plt.close(fig)
            buf.seek(0)
            zf.writestr(f"helical_wheel_{idx:03d}.png", buf.read())
    mem.seek(0)
    return mem.read()

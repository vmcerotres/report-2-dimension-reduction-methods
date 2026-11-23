import csv
import math
import os
import shutil
import subprocess

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(
    BASE_DIR,
    "forty+soybean+cultivars+from+subsequent+harvests",
    "data.csv",
)

# Continuous feature columns (omitting target variable 'GY' and experimental identifiers)
FEATURE_COLS = ["PH", "IFP", "NLP", "NGP", "NGL", "NS", "MHG"]
TARGET_COL = "GY"


def load_data(data_path: str, feature_cols: list[str]):
    with open(data_path, mode="r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = []
        seasons = []
        targets = []
        for row in reader:
            rows.append([float(row[c]) for c in feature_cols])
            seasons.append(int(row["Season"]))
            targets.append(float(row[TARGET_COL]))
    return rows, seasons, targets


def compute_stats_and_matrices(rows: list[list[float]], feature_cols: list[str]):
    n = len(rows)
    p = len(feature_cols)
    means = [sum(rows[i][j] for i in range(n)) / n for j in range(p)]

    # Sample covariance matrix (ddof=1)
    cov_mat = [[0.0] * p for _ in range(p)]
    for j1 in range(p):
        for j2 in range(j1, p):
            cov = sum(
                (rows[i][j1] - means[j1]) * (rows[i][j2] - means[j2])
                for i in range(n)
            ) / (n - 1)
            cov_mat[j1][j2] = cov
            cov_mat[j2][j1] = cov

    stds = [math.sqrt(cov_mat[j][j]) for j in range(p)]

    # Pearson correlation matrix
    corr_mat = [[0.0] * p for _ in range(p)]
    for j1 in range(p):
        for j2 in range(p):
            denom = stds[j1] * stds[j2]
            corr_mat[j1][j2] = cov_mat[j1][j2] / denom if denom != 0 else 0.0

    # Centered data (for Covariance PCA) and Standardized data (for Correlation PCA)
    centered_rows = [
        [rows[i][j] - means[j] for j in range(p)]
        for i in range(n)
    ]
    standardized_rows = [
        [(rows[i][j] - means[j]) / stds[j] for j in range(p)]
        for i in range(n)
    ]

    return means, stds, cov_mat, corr_mat, centered_rows, standardized_rows


def jacobi_eigen(matrix: list[list[float]], tol: float = 1e-12, max_iter: int = 1000):
    """
    Computes eigenvalues and eigenvectors of a real symmetric matrix using the Jacobi method.
    Returns (eigenvalues, eigenvectors_columns) sorted in descending order of eigenvalues.
    Each eigenvector has its largest-magnitude element forced positive for deterministic signs.
    """
    p = len(matrix)
    A = [row[:] for row in matrix]
    V = [[1.0 if i == j else 0.0 for j in range(p)] for i in range(p)]

    for _ in range(max_iter):
        # Find largest off-diagonal element
        max_val = 0.0
        k, l = 0, 1
        for i in range(p):
            for j in range(i + 1, p):
                if abs(A[i][j]) > max_val:
                    max_val = abs(A[i][j])
                    k, l = i, j
        if max_val < tol:
            break

        diff = A[l][l] - A[k][k]
        if abs(A[k][l]) < abs(diff) * 1e-36:
            t = A[k][l] / diff
        else:
            phi = diff / (2.0 * A[k][l])
            t = 1.0 / (abs(phi) + math.sqrt(phi * phi + 1.0))
            if phi < 0.0:
                t = -t

        c = 1.0 / math.sqrt(t * t + 1.0)
        s = t * c
        tau = s / (1.0 + c)
        h = t * A[k][l]

        A[k][k] -= h
        A[l][l] += h
        A[k][l] = 0.0
        A[l][k] = 0.0

        for i in range(k):
            g = A[i][k]
            h_val = A[i][l]
            A[i][k] = g - s * (h_val + g * tau)
            A[i][l] = h_val + s * (g - h_val * tau)
            A[k][i] = A[i][k]
            A[l][i] = A[i][l]

        for i in range(k + 1, l):
            g = A[k][i]
            h_val = A[i][l]
            A[k][i] = g - s * (h_val + g * tau)
            A[i][l] = h_val + s * (g - h_val * tau)
            A[i][k] = A[k][i]
            A[l][i] = A[i][l]

        for i in range(l + 1, p):
            g = A[k][i]
            h_val = A[l][i]
            A[k][i] = g - s * (h_val + g * tau)
            A[l][i] = h_val + s * (g - h_val * tau)
            A[i][k] = A[k][i]
            A[i][l] = A[l][i]

        for i in range(p):
            g = V[i][k]
            h_val = V[i][l]
            V[i][k] = g - s * (h_val + g * tau)
            V[i][l] = h_val + s * (g - h_val * tau)

    eigenvalues = [A[i][i] for i in range(p)]
    # Extract columns of V as eigenvectors
    eigenvectors = [[V[row][col] for row in range(p)] for col in range(p)]

    # Sort descending by eigenvalue
    order = sorted(range(p), key=lambda idx: eigenvalues[idx], reverse=True)
    sorted_vals = [max(0.0, eigenvalues[idx]) for idx in order]
    sorted_vecs = []
    for idx in order:
        vec = eigenvectors[idx]
        # Deterministic sign: make the entry with largest absolute value positive
        max_comp = max(range(p), key=lambda r: abs(vec[r]))
        if vec[max_comp] < 0:
            vec = [-x for x in vec]
        sorted_vecs.append(vec)

    return sorted_vals, sorted_vecs


def perform_pca(
    matrix: list[list[float]],
    transformed_rows: list[list[float]],
    stds: list[float],
    is_corr: bool,
):
    """
    Runs PCA on either Covariance or Correlation matrix.
    Returns dict with:
      - eigenvalues
      - explained_variance_pct
      - cumulative_variance_pct
      - eigenvectors: list of p vectors (each of length p), representing the linear combination weights e_k
      - factor_loadings: p x p matrix (variable i, component k) = corr(X_i, PC_k)
          For correlation PCA: loading_{i,k} = e_{i,k} * sqrt(lambda_k)
          For covariance PCA: loading_{i,k} = e_{i,k} * sqrt(lambda_k) / std_i
      - scores: n x p matrix of principal component scores
    """
    p = len(matrix)
    eigenvalues, eigenvectors = jacobi_eigen(matrix)
    total_var = sum(eigenvalues)
    explained_pct = [(val / total_var) * 100.0 for val in eigenvalues]
    cum_pct = []
    running = 0.0
    for pct in explained_pct:
        running += pct
        cum_pct.append(running)

    # Factor loadings (correlations between original variables and PCs)
    # Also store raw eigenvectors (component weights / loadings in sklearn terminology)
    factor_loadings = [[0.0] * p for _ in range(p)]  # row = variable i, col = PC k
    for k in range(p):
        sqrt_lam = math.sqrt(eigenvalues[k])
        for i in range(p):
            if is_corr:
                factor_loadings[i][k] = eigenvectors[k][i] * sqrt_lam
            else:
                factor_loadings[i][k] = (eigenvectors[k][i] * sqrt_lam) / stds[i]

    # Project rows onto eigenvectors: scores[row][k]
    scores = []
    for row in transformed_rows:
        row_scores = [
            sum(row[j] * eigenvectors[k][j] for j in range(p))
            for k in range(p)
        ]
        scores.append(row_scores)

    return {
        "eigenvalues": eigenvalues,
        "explained_pct": explained_pct,
        "cumulative_pct": cum_pct,
        "eigenvectors": eigenvectors,  # eigenvectors[k][i] = weight of variable i in PC_{k+1}
        "factor_loadings": factor_loadings,  # factor_loadings[i][k] = correlation of var i with PC_{k+1}
        "scores": scores,
    }


def _svg_to_png(svg_lines: list[str], output_png: str, size_px: int = 1200):
    svg_path = output_png.replace(".png", ".svg")
    with open(svg_path, "w", encoding="utf-8") as f:
        f.write("\n".join(svg_lines))

    out_dir = os.path.dirname(os.path.abspath(output_png))
    subprocess.run(
        ["qlmanage", "-t", "-s", str(size_px), "-o", out_dir, svg_path],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    ql_png = svg_path + ".png"
    if os.path.exists(ql_png):
        shutil.move(ql_png, output_png)
        os.remove(svg_path)


def _lerp_color(c1: tuple[int, int, int], c2: tuple[int, int, int], t: float) -> tuple[int, int, int]:
    t = max(0.0, min(1.0, t))
    return (
        int(round(c1[0] + (c2[0] - c1[0]) * t)),
        int(round(c1[1] + (c2[1] - c1[1]) * t)),
        int(round(c1[2] + (c2[2] - c1[2]) * t)),
    )


def _diverging_color(val: float, vmin: float, vmax: float) -> tuple[str, str]:
    blue = (43, 92, 170)
    light_blue = (158, 202, 225)
    white = (247, 247, 247)
    light_red = (252, 174, 145)
    red = (185, 29, 45)

    if val < 0 and vmin < 0:
        ratio = min(1.0, abs(val) / abs(vmin))
        rgb = _lerp_color(white, light_blue, ratio * 2.0) if ratio < 0.5 else _lerp_color(light_blue, blue, (ratio - 0.5) * 2.0)
    elif val > 0 and vmax > 0:
        ratio = min(1.0, val / vmax)
        rgb = _lerp_color(white, light_red, ratio * 2.0) if ratio < 0.5 else _lerp_color(light_red, red, (ratio - 0.5) * 2.0)
    else:
        rgb = white

    bg_hex = f"#{rgb[0]:02x}{rgb[1]:02x}{rgb[2]:02x}"
    luminance = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]
    text_hex = "#ffffff" if luminance < 145 else "#1e293b"
    return bg_hex, text_hex


def save_heatmap_image(
    matrix: list[list[float]],
    cols: list[str],
    title: str,
    subtitle: str,
    output_png: str,
    is_corr: bool = True,
):
    size = 900
    n = len(cols)
    grid_x = 120
    grid_y = 140
    cell_size = 92
    grid_size = n * cell_size

    all_vals = [matrix[r][c] for r in range(n) for c in range(n)]
    if is_corr:
        vmin, vmax = -1.0, 1.0
    else:
        vmin = -math.log1p(max(abs(v) for v in all_vals))
        vmax = math.log1p(max(abs(v) for v in all_vals))

    svg_lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{size}" height="{size}" viewBox="0 0 {size} {size}">',
        f'  <rect width="{size}" height="{size}" fill="#ffffff"/>',
        f'  <text x="{size // 2}" y="58" font-family="Helvetica, Arial, sans-serif" font-size="26" font-weight="bold" fill="#0f172a" text-anchor="middle">{title}</text>',
        f'  <text x="{size // 2}" y="88" font-family="Helvetica, Arial, sans-serif" font-size="16" fill="#475569" text-anchor="middle">{subtitle}</text>',
    ]

    for idx, col in enumerate(cols):
        cx = grid_x + idx * cell_size + cell_size / 2
        cy = grid_y + idx * cell_size + cell_size / 2
        svg_lines.append(
            f'  <text x="{cx}" y="{grid_y - 16}" font-family="Helvetica, Arial, sans-serif" font-size="16" font-weight="bold" fill="#1e293b" text-anchor="middle">{col}</text>'
        )
        svg_lines.append(
            f'  <text x="{cx}" y="{grid_y + grid_size + 28}" font-family="Helvetica, Arial, sans-serif" font-size="16" font-weight="bold" fill="#1e293b" text-anchor="middle">{col}</text>'
        )
        svg_lines.append(
            f'  <text x="{grid_x - 16}" y="{cy + 6}" font-family="Helvetica, Arial, sans-serif" font-size="16" font-weight="bold" fill="#1e293b" text-anchor="end">{col}</text>'
        )

    for i in range(n):
        for j in range(n):
            val = matrix[i][j]
            if is_corr:
                scaled_val = val
                val_str = f"{val:.2f}"
            else:
                scaled_val = math.copysign(math.log1p(abs(val)), val)
                val_str = f"{val:.1f}" if abs(val) >= 100 else f"{val:.2f}"

            bg_hex, text_hex = _diverging_color(scaled_val, vmin, vmax)
            x = grid_x + j * cell_size
            y = grid_y + i * cell_size
            svg_lines.append(
                f'  <rect x="{x}" y="{y}" width="{cell_size}" height="{cell_size}" fill="{bg_hex}" stroke="#e2e8f0" stroke-width="1.5"/>'
            )
            svg_lines.append(
                f'  <text x="{x + cell_size / 2}" y="{y + cell_size / 2 + 6}" font-family="Helvetica, Arial, sans-serif" font-size="16" font-weight="bold" fill="{text_hex}" text-anchor="middle">{val_str}</text>'
            )

    svg_lines.append(
        f'  <rect x="{grid_x}" y="{grid_y}" width="{grid_size}" height="{grid_size}" fill="none" stroke="#334155" stroke-width="2"/>'
    )

    cb_x = grid_x + grid_size + 38
    cb_y = grid_y
    cb_w = 24
    cb_h = grid_size
    steps = 100
    step_h = cb_h / steps
    for s in range(steps):
        frac = 1.0 - (s / (steps - 1))
        c_val = vmin + frac * (vmax - vmin)
        bg_hex, _ = _diverging_color(c_val, vmin, vmax)
        sy = cb_y + s * step_h
        svg_lines.append(
            f'  <rect x="{cb_x}" y="{sy:.2f}" width="{cb_w}" height="{step_h + 1:.2f}" fill="{bg_hex}" stroke="none"/>'
        )
    svg_lines.append(
        f'  <rect x="{cb_x}" y="{cb_y}" width="{cb_w}" height="{cb_h}" fill="none" stroke="#334155" stroke-width="1.5"/>'
    )

    if is_corr:
        ticks = [(0.0, "1.00"), (0.25, "0.50"), (0.50, "0.00"), (0.75, "-0.50"), (1.0, "-1.00")]
    else:
        max_raw = max(abs(v) for v in all_vals)
        ticks = [(0.0, f"+{max_raw:.0f}"), (0.25, "+100"), (0.50, "0"), (0.75, "-10"), (1.0, f"-{max_raw:.0f}")]

    for frac, label in ticks:
        ty = cb_y + frac * cb_h
        svg_lines.append(
            f'  <line x1="{cb_x + cb_w}" y1="{ty}" x2="{cb_x + cb_w + 6}" y2="{ty}" stroke="#334155" stroke-width="1.5"/>'
        )
        svg_lines.append(
            f'  <text x="{cb_x + cb_w + 10}" y="{ty + 5}" font-family="Helvetica, Arial, sans-serif" font-size="13" fill="#334155">{label}</text>'
        )

    svg_lines.append("</svg>")
    _svg_to_png(svg_lines, output_png)


def save_scree_plot(pca_res: dict, title: str, subtitle: str, output_png: str):
    width, height = 900, 900
    plot_x, plot_y = 110, 175
    plot_w, plot_h = 720, 560

    explained = pca_res["explained_pct"]
    cum_pct = pca_res["cumulative_pct"]
    eigenvalues = pca_res["eigenvalues"]
    p = len(explained)

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'  <rect width="{width}" height="{height}" fill="#ffffff"/>',
        f'  <text x="{width // 2}" y="68" font-family="Helvetica, Arial, sans-serif" font-size="26" font-weight="bold" fill="#0f172a" text-anchor="middle">{title}</text>',
        f'  <text x="{width // 2}" y="102" font-family="Helvetica, Arial, sans-serif" font-size="16" fill="#475569" text-anchor="middle">{subtitle}</text>',
    ]

    # Horizontal grid lines (0% to 100%)
    for pct_tick in [0, 20, 40, 60, 80, 90, 95, 100]:
        y = plot_y + plot_h - (pct_tick / 100.0) * plot_h
        is_thresh = pct_tick in (90, 95)
        stroke = "#dc2626" if is_thresh else "#e2e8f0"
        dash = ' stroke-dasharray="6,4"' if is_thresh else ""
        svg.append(
            f'  <line x1="{plot_x}" y1="{y:.1f}" x2="{plot_x + plot_w}" y2="{y:.1f}" stroke="{stroke}" stroke-width="1.5"{dash}/>'
        )
        label_color = "#dc2626" if is_thresh else "#475569"
        font_w = ' font-weight="bold"' if is_thresh else ""
        svg.append(
            f'  <text x="{plot_x - 14}" y="{y + 5:.1f}" font-family="Helvetica, Arial, sans-serif" font-size="13"{font_w} fill="{label_color}" text-anchor="end">{pct_tick}%</text>'
        )

    # Y-axis title
    svg.append(
        f'  <text x="38" y="{plot_y + plot_h / 2}" font-family="Helvetica, Arial, sans-serif" font-size="17" font-weight="bold" fill="#1e293b" text-anchor="middle" transform="rotate(-90 38 {plot_y + plot_h / 2})">Explained Variance (%)</text>'
    )
    # X-axis title
    svg.append(
        f'  <text x="{plot_x + plot_w / 2}" y="{plot_y + plot_h + 80}" font-family="Helvetica, Arial, sans-serif" font-size="17" font-weight="bold" fill="#1e293b" text-anchor="middle">Principal Component</text>'
    )

    col_w = plot_w / p
    bar_w = col_w * 0.48
    cum_points = []

    for i in range(p):
        cx = plot_x + (i + 0.5) * col_w
        bar_h = (explained[i] / 100.0) * plot_h
        bx = cx - bar_w / 2
        by = plot_y + plot_h - bar_h

        # Bar for individual explained variance
        svg.append(
            f'  <rect x="{bx:.1f}" y="{by:.1f}" width="{bar_w:.1f}" height="{max(1.5, bar_h):.1f}" fill="#3b82f6" fill-opacity="0.82" rx="3"/>'
        )
        # Value label: for PC1 where cumulative == individual, place bar label inside/below or only show cumulative on top
        if i > 0:
            svg.append(
                f'  <text x="{cx:.1f}" y="{by - 8:.1f}" font-family="Helvetica, Arial, sans-serif" font-size="12" font-weight="bold" fill="#1d4ed8" text-anchor="middle">{explained[i]:.1f}%</text>'
            )
        else:
            svg.append(
                f'  <text x="{cx:.1f}" y="{by + 24:.1f}" font-family="Helvetica, Arial, sans-serif" font-size="12" font-weight="bold" fill="#ffffff" text-anchor="middle">{explained[i]:.1f}%</text>'
            )
        # X-axis PC label & eigenvalue
        svg.append(
            f'  <text x="{cx:.1f}" y="{plot_y + plot_h + 26}" font-family="Helvetica, Arial, sans-serif" font-size="15" font-weight="bold" fill="#1e293b" text-anchor="middle">PC{i + 1}</text>'
        )
        lam_str = f"λ={eigenvalues[i]:.1f}" if eigenvalues[i] >= 10 else f"λ={eigenvalues[i]:.2f}"
        svg.append(
            f'  <text x="{cx:.1f}" y="{plot_y + plot_h + 46}" font-family="Helvetica, Arial, sans-serif" font-size="12" fill="#64748b" text-anchor="middle">({lam_str})</text>'
        )

        cy_cum = plot_y + plot_h - (cum_pct[i] / 100.0) * plot_h
        cum_points.append((cx, cy_cum, cum_pct[i]))

    # Cumulative variance polyline
    poly_str = " ".join(f"{x:.1f},{y:.1f}" for x, y, _ in cum_points)
    svg.append(
        f'  <polyline points="{poly_str}" fill="none" stroke="#ea580c" stroke-width="3.2"/>'
    )
    for x, y, pct in cum_points:
        svg.append(
            f'  <circle cx="{x:.1f}" cy="{y:.1f}" r="5.5" fill="#ea580c" stroke="#ffffff" stroke-width="2"/>'
        )
        svg.append(
            f'  <text x="{x:.1f}" y="{y - 12:.1f}" font-family="Helvetica, Arial, sans-serif" font-size="12" font-weight="bold" fill="#c2410c" text-anchor="middle">{pct:.1f}%</text>'
        )

    # Plot frame
    svg.append(
        f'  <rect x="{plot_x}" y="{plot_y}" width="{plot_w}" height="{plot_h}" fill="none" stroke="#334155" stroke-width="2"/>'
    )

    # Legend box placed in bottom-right of plot area (where bars are low and cumulative line is at top)
    leg_x, leg_y = plot_x + plot_w - 255, plot_y + plot_h - 140
    svg.append(
        f'  <rect x="{leg_x}" y="{leg_y}" width="235" height="82" fill="#f8fafc" stroke="#cbd5e1" stroke-width="1.5" rx="6"/>'
    )
    svg.append(
        f'  <rect x="{leg_x + 14}" y="{leg_y + 14}" width="18" height="14" fill="#3b82f6" fill-opacity="0.82" rx="2"/>'
    )
    svg.append(
        f'  <text x="{leg_x + 42}" y="{leg_y + 26}" font-family="Helvetica, Arial, sans-serif" font-size="13" fill="#1e293b">Individual Variance (%)</text>'
    )
    svg.append(
        f'  <line x1="{leg_x + 14}" y1="{leg_y + 45}" x2="{leg_x + 32}" y2="{leg_y + 45}" stroke="#ea580c" stroke-width="3"/>'
    )
    svg.append(
        f'  <circle cx="{leg_x + 23}" cy="{leg_y + 45}" r="4.5" fill="#ea580c"/>'
    )
    svg.append(
        f'  <text x="{leg_x + 42}" y="{leg_y + 50}" font-family="Helvetica, Arial, sans-serif" font-size="13" fill="#1e293b">Cumulative Variance (%)</text>'
    )
    svg.append(
        f'  <line x1="{leg_x + 14}" y1="{leg_y + 66}" x2="{leg_x + 32}" y2="{leg_y + 66}" stroke="#dc2626" stroke-width="1.5" stroke-dasharray="5,3"/>'
    )
    svg.append(
        f'  <text x="{leg_x + 42}" y="{leg_y + 70}" font-family="Helvetica, Arial, sans-serif" font-size="12" fill="#dc2626">90% / 95% Thresholds</text>'
    )

    svg.append("</svg>")
    _svg_to_png(svg, output_png)


def save_2d_pca_plot(
    pca_res: dict,
    seasons: list[int],
    feature_cols: list[str],
    title: str,
    subtitle: str,
    output_png: str,
):
    width, height = 900, 900
    plot_x, plot_y = 115, 160
    plot_w, plot_h = 720, 620

    scores = pca_res["scores"]
    explained = pca_res["explained_pct"]
    loadings = pca_res["factor_loadings"]

    xs = [s[0] for s in scores]
    ys = [s[1] for s in scores]

    max_x = max(max(abs(v) for v in xs) * 1.15, 1e-6)
    max_y = max(max(abs(v) for v in ys) * 1.15, 1e-6)

    def to_px(x_val: float, y_val: float) -> tuple[float, float]:
        px = plot_x + ((x_val + max_x) / (2.0 * max_x)) * plot_w
        py = plot_y + plot_h - ((y_val + max_y) / (2.0 * max_y)) * plot_h
        return px, py

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'  <rect width="{width}" height="{height}" fill="#ffffff"/>',
        f'  <text x="{width // 2}" y="66" font-family="Helvetica, Arial, sans-serif" font-size="25" font-weight="bold" fill="#0f172a" text-anchor="middle">{title}</text>',
        f'  <text x="{width // 2}" y="98" font-family="Helvetica, Arial, sans-serif" font-size="16" fill="#475569" text-anchor="middle">{subtitle}</text>',
    ]

    # Grid lines and ticks
    for i in range(-4, 5):
        frac = i / 4.0
        xv = frac * max_x
        yv = frac * max_y
        px, _ = to_px(xv, 0.0)
        _, py = to_px(0.0, yv)

        stroke = "#94a3b8" if i == 0 else "#e2e8f0"
        sw = "2" if i == 0 else "1.2"
        dash = ' stroke-dasharray="5,4"' if i == 0 else ""

        svg.append(
            f'  <line x1="{px:.1f}" y1="{plot_y}" x2="{px:.1f}" y2="{plot_y + plot_h}" stroke="{stroke}" stroke-width="{sw}"{dash}/>'
        )
        svg.append(
            f'  <line x1="{plot_x}" y1="{py:.1f}" x2="{plot_x + plot_w}" y2="{py:.1f}" stroke="{stroke}" stroke-width="{sw}"{dash}/>'
        )

        x_lbl = f"{xv:.0f}" if max_x >= 20 else f"{xv:.2f}"
        y_lbl = f"{yv:.0f}" if max_y >= 20 else f"{yv:.2f}"
        svg.append(
            f'  <text x="{px:.1f}" y="{plot_y + plot_h + 26}" font-family="Helvetica, Arial, sans-serif" font-size="13" fill="#475569" text-anchor="middle">{x_lbl}</text>'
        )
        svg.append(
            f'  <text x="{plot_x - 12}" y="{py + 4:.1f}" font-family="Helvetica, Arial, sans-serif" font-size="13" fill="#475569" text-anchor="end">{y_lbl}</text>'
        )

    # Scatter points colored by Season
    for (s, season) in zip(scores, seasons):
        px, py = to_px(s[0], s[1])
        color = "#2563eb" if season == 1 else "#16a34a"
        svg.append(
            f'  <circle cx="{px:.1f}" cy="{py:.1f}" r="5.0" fill="{color}" fill-opacity="0.72" stroke="#ffffff" stroke-width="0.9"/>'
        )

    # Overlay loading vectors (biplot arrows scaled to 75% of plot radius)
    cx0, cy0 = to_px(0.0, 0.0)
    svg.append(
        '  <defs><marker id="arrow" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse"><path d="M 0 1 L 10 5 L 0 9 z" fill="#dc2626"/></marker></defs>'
    )
    for j, col in enumerate(feature_cols):
        lx = loadings[j][0]
        ly = loadings[j][1]
        if math.hypot(lx, ly) < 0.12:
            continue
        ax, ay = to_px(lx * max_x * 0.78, ly * max_y * 0.78)
        svg.append(
            f'  <line x1="{cx0:.1f}" y1="{cy0:.1f}" x2="{ax:.1f}" y2="{ay:.1f}" stroke="#dc2626" stroke-width="2.2" marker-end="url(#arrow)"/>'
        )
        tx, ty = to_px(lx * max_x * 0.86, ly * max_y * 0.86)
        svg.append(
            f'  <text x="{tx:.1f}" y="{ty + 4:.1f}" font-family="Helvetica, Arial, sans-serif" font-size="14" font-weight="bold" fill="#b91c1c" text-anchor="middle">{col}</text>'
        )

    # Frame
    svg.append(
        f'  <rect x="{plot_x}" y="{plot_y}" width="{plot_w}" height="{plot_h}" fill="none" stroke="#334155" stroke-width="2"/>'
    )

    # Axis titles
    svg.append(
        f'  <text x="{plot_x + plot_w / 2}" y="{plot_y + plot_h + 68}" font-family="Helvetica, Arial, sans-serif" font-size="17" font-weight="bold" fill="#0f172a" text-anchor="middle">PC1 ({explained[0]:.2f}% Explained Variance)</text>'
    )
    svg.append(
        f'  <text x="36" y="{plot_y + plot_h / 2}" font-family="Helvetica, Arial, sans-serif" font-size="17" font-weight="bold" fill="#0f172a" text-anchor="middle" transform="rotate(-90 36 {plot_y + plot_h / 2})">PC2 ({explained[1]:.2f}% Explained Variance)</text>'
    )

    # Legend inside top-left of plot area
    leg_x, leg_y = plot_x + 18, plot_y + 18
    svg.append(
        f'  <rect x="{leg_x}" y="{leg_y}" width="180" height="86" fill="#f8fafc" stroke="#cbd5e1" stroke-width="1.5" rx="6"/>'
    )
    svg.append(f'  <circle cx="{leg_x + 20}" cy="{leg_y + 22}" r="5.5" fill="#2563eb" fill-opacity="0.8"/>')
    svg.append(f'  <text x="{leg_x + 36}" y="{leg_y + 26}" font-family="Helvetica, Arial, sans-serif" font-size="13" fill="#1e293b">Season 1 (n=160)</text>')
    svg.append(f'  <circle cx="{leg_x + 20}" cy="{leg_y + 46}" r="5.5" fill="#16a34a" fill-opacity="0.8"/>')
    svg.append(f'  <text x="{leg_x + 36}" y="{leg_y + 50}" font-family="Helvetica, Arial, sans-serif" font-size="13" fill="#1e293b">Season 2 (n=160)</text>')
    svg.append(f'  <line x1="{leg_x + 12}" y1="{leg_y + 68}" x2="{leg_x + 28}" y2="{leg_y + 68}" stroke="#dc2626" stroke-width="2.2"/>')
    svg.append(f'  <text x="{leg_x + 36}" y="{leg_y + 72}" font-family="Helvetica, Arial, sans-serif" font-size="13" fill="#b91c1c">Feature Loadings</text>')

    svg.append("</svg>")
    _svg_to_png(svg, output_png)


def save_3d_pca_plot(
    pca_res: dict,
    seasons: list[int],
    title: str,
    subtitle: str,
    output_png: str,
):
    width, height = 900, 900
    cx_center, cy_center = 450, 485
    scale = 215.0

    scores = pca_res["scores"]
    explained = pca_res["explained_pct"]

    xs = [s[0] for s in scores]
    ys = [s[1] for s in scores]
    zs = [s[2] for s in scores]

    max_x = max(max(abs(v) for v in xs) * 1.1, 1e-6)
    max_y = max(max(abs(v) for v in ys) * 1.1, 1e-6)
    max_z = max(max(abs(v) for v in zs) * 1.1, 1e-6)

    # Isometric /axonometric 3D projection:
    # u in [-1, 1] (PC1), v in [-1, 1] (PC2), w in [-1, 1] (PC3, vertical)
    az = math.radians(38.0)
    el = math.radians(24.0)

    def project_3d(u: float, v: float, w: float) -> tuple[float, float, float]:
        # Rotate around vertical w-axis by az, then tilt by el
        x_rot = u * math.cos(az) - v * math.sin(az)
        y_rot = u * math.sin(az) + v * math.cos(az)
        px = cx_center + scale * x_rot
        py = cy_center - scale * (w * math.cos(el) - y_rot * math.sin(el))
        depth = y_rot * math.cos(el) + w * math.sin(el)
        return px, py, depth

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'  <rect width="{width}" height="{height}" fill="#ffffff"/>',
        f'  <text x="{width // 2}" y="66" font-family="Helvetica, Arial, sans-serif" font-size="24" font-weight="bold" fill="#0f172a" text-anchor="middle">{title}</text>',
        f'  <text x="{width // 2}" y="98" font-family="Helvetica, Arial, sans-serif" font-size="16" fill="#475569" text-anchor="middle">{subtitle}</text>',
    ]

    # Draw 3D bounding cube and back grid planes
    # Back planes: u = -1, v = 1, w = -1
    grid_ticks = [-1.0, -0.5, 0.0, 0.5, 1.0]

    # Floor grid (w = -1)
    for t in grid_ticks:
        p1 = project_3d(t, -1.0, -1.0)
        p2 = project_3d(t, 1.0, -1.0)
        svg.append(f'  <line x1="{p1[0]:.1f}" y1="{p1[1]:.1f}" x2="{p2[0]:.1f}" y2="{p2[1]:.1f}" stroke="#e2e8f0" stroke-width="1.2"/>')
        p3 = project_3d(-1.0, t, -1.0)
        p4 = project_3d(1.0, t, -1.0)
        svg.append(f'  <line x1="{p3[0]:.1f}" y1="{p3[1]:.1f}" x2="{p4[0]:.1f}" y2="{p4[1]:.1f}" stroke="#e2e8f0" stroke-width="1.2"/>')

    # Back-left wall (u = -1)
    for t in grid_ticks:
        p1 = project_3d(-1.0, t, -1.0)
        p2 = project_3d(-1.0, t, 1.0)
        svg.append(f'  <line x1="{p1[0]:.1f}" y1="{p1[1]:.1f}" x2="{p2[0]:.1f}" y2="{p2[1]:.1f}" stroke="#f1f5f9" stroke-width="1.2"/>')
        p3 = project_3d(-1.0, -1.0, t)
        p4 = project_3d(-1.0, 1.0, t)
        svg.append(f'  <line x1="{p3[0]:.1f}" y1="{p3[1]:.1f}" x2="{p4[0]:.1f}" y2="{p4[1]:.1f}" stroke="#f1f5f9" stroke-width="1.2"/>')

    # Back-right wall (v = 1)
    for t in grid_ticks:
        p1 = project_3d(t, 1.0, -1.0)
        p2 = project_3d(t, 1.0, 1.0)
        svg.append(f'  <line x1="{p1[0]:.1f}" y1="{p1[1]:.1f}" x2="{p2[0]:.1f}" y2="{p2[1]:.1f}" stroke="#f1f5f9" stroke-width="1.2"/>')
        p3 = project_3d(-1.0, 1.0, t)
        p4 = project_3d(1.0, 1.0, t)
        svg.append(f'  <line x1="{p3[0]:.1f}" y1="{p3[1]:.1f}" x2="{p4[0]:.1f}" y2="{p4[1]:.1f}" stroke="#f1f5f9" stroke-width="1.2"/>')

    # Origin axes inside cube (u=0, v=0, w=0)
    o_u1, o_u2 = project_3d(-1.0, 0.0, 0.0), project_3d(1.05, 0.0, 0.0)
    o_v1, o_v2 = project_3d(0.0, -1.0, 0.0), project_3d(0.0, 1.05, 0.0)
    o_w1, o_w2 = project_3d(0.0, 0.0, -1.0), project_3d(0.0, 0.0, 1.05)

    svg.append(f'  <line x1="{o_u1[0]:.1f}" y1="{o_u1[1]:.1f}" x2="{o_u2[0]:.1f}" y2="{o_u2[1]:.1f}" stroke="#64748b" stroke-width="1.6" stroke-dasharray="5,4"/>')
    svg.append(f'  <line x1="{o_v1[0]:.1f}" y1="{o_v1[1]:.1f}" x2="{o_v2[0]:.1f}" y2="{o_v2[1]:.1f}" stroke="#64748b" stroke-width="1.6" stroke-dasharray="5,4"/>')
    svg.append(f'  <line x1="{o_w1[0]:.1f}" y1="{o_w1[1]:.1f}" x2="{o_w2[0]:.1f}" y2="{o_w2[1]:.1f}" stroke="#64748b" stroke-width="1.6" stroke-dasharray="5,4"/>')

    # Outer frame edges of the 3D box
    box_edges = [
        ((-1, -1, -1), (1, -1, -1)),
        ((1, -1, -1), (1, 1, -1)),
        ((-1, -1, -1), (-1, 1, -1)),
        ((-1, 1, -1), (1, 1, -1)),
        ((-1, -1, -1), (-1, -1, 1)),
        ((1, -1, -1), (1, -1, 1)),
        ((1, 1, -1), (1, 1, 1)),
        ((-1, 1, -1), (-1, 1, 1)),
        ((-1, -1, 1), (1, -1, 1)),
        ((1, -1, 1), (1, 1, 1)),
        ((-1, -1, 1), (-1, 1, 1)),
        ((-1, 1, 1), (1, 1, 1)),
    ]
    for c1, c2 in box_edges:
        p1 = project_3d(*c1)
        p2 = project_3d(*c2)
        svg.append(f'  <line x1="{p1[0]:.1f}" y1="{p1[1]:.1f}" x2="{p2[0]:.1f}" y2="{p2[1]:.1f}" stroke="#94a3b8" stroke-width="1.5"/>')

    # Axis ticks and labels on outer edges
    for frac in [-1.0, 0.0, 1.0]:
        # PC1 axis along (u, -1, -1)
        p_u = project_3d(frac, -1.0, -1.0)
        val_u = frac * max_x
        lbl_u = f"{val_u:.0f}" if max_x >= 20 else f"{val_u:.1f}"
        svg.append(f'  <text x="{p_u[0] - 8:.1f}" y="{p_u[1] + 20:.1f}" font-family="Helvetica, Arial, sans-serif" font-size="12" fill="#475569" text-anchor="middle">{lbl_u}</text>')

        # PC2 axis along (1, v, -1)
        p_v = project_3d(1.0, frac, -1.0)
        val_v = frac * max_y
        lbl_v = f"{val_v:.0f}" if max_y >= 20 else f"{val_v:.1f}"
        svg.append(f'  <text x="{p_v[0] + 18:.1f}" y="{p_v[1] + 14:.1f}" font-family="Helvetica, Arial, sans-serif" font-size="12" fill="#475569" text-anchor="start">{lbl_v}</text>')

        # PC3 axis along (-1, -1, w)
        p_w = project_3d(-1.0, -1.0, frac)
        val_w = frac * max_z
        lbl_w = f"{val_w:.0f}" if max_z >= 20 else f"{val_w:.1f}"
        svg.append(f'  <text x="{p_w[0] - 14:.1f}" y="{p_w[1] + 4:.1f}" font-family="Helvetica, Arial, sans-serif" font-size="12" fill="#475569" text-anchor="end">{lbl_w}</text>')

    # Project and depth-sort all data points (painter's algorithm: back to front)
    projected_pts = []
    for s, season in zip(scores, seasons):
        u = s[0] / max_x
        v = s[1] / max_y
        w = s[2] / max_z
        px, py, depth = project_3d(u, v, w)
        floor_px, floor_py, _ = project_3d(u, v, -1.0)
        projected_pts.append((depth, px, py, floor_px, floor_py, season))

    projected_pts.sort(key=lambda item: item[0], reverse=True)

    # Subtle drop shadows on the w=-1 floor for 3D depth perception
    for _, _, _, fpx, fpy, _ in projected_pts:
        svg.append(
            f'  <circle cx="{fpx:.1f}" cy="{fpy:.1f}" r="2.2" fill="#cbd5e1" fill-opacity="0.55"/>'
        )

    # 3D points
    for _, px, py, _, _, season in projected_pts:
        color = "#2563eb" if season == 1 else "#16a34a"
        svg.append(
            f'  <circle cx="{px:.1f}" cy="{py:.1f}" r="5.0" fill="{color}" fill-opacity="0.78" stroke="#ffffff" stroke-width="0.9"/>'
        )

    # Axis titles
    lbl_pc1 = project_3d(0.0, -1.32, -1.15)
    lbl_pc2 = project_3d(1.32, 0.0, -1.15)
    lbl_pc3 = project_3d(-1.32, -1.15, 0.0)

    svg.append(
        f'  <text x="{lbl_pc1[0]:.1f}" y="{lbl_pc1[1]:.1f}" font-family="Helvetica, Arial, sans-serif" font-size="16" font-weight="bold" fill="#0f172a" text-anchor="middle">PC1 ({explained[0]:.2f}%)</text>'
    )
    svg.append(
        f'  <text x="{lbl_pc2[0]:.1f}" y="{lbl_pc2[1]:.1f}" font-family="Helvetica, Arial, sans-serif" font-size="16" font-weight="bold" fill="#0f172a" text-anchor="middle">PC2 ({explained[1]:.2f}%)</text>'
    )
    svg.append(
        f'  <text x="{lbl_pc3[0]:.1f}" y="{lbl_pc3[1]:.1f}" font-family="Helvetica, Arial, sans-serif" font-size="16" font-weight="bold" fill="#0f172a" text-anchor="end">PC3 ({explained[2]:.2f}%)</text>'
    )

    # Legend inside top-right with safe margin
    leg_x, leg_y = width - 250, 140
    svg.append(
        f'  <rect x="{leg_x}" y="{leg_y}" width="185" height="86" fill="#f8fafc" stroke="#cbd5e1" stroke-width="1.5" rx="6"/>'
    )
    svg.append(f'  <circle cx="{leg_x + 20}" cy="{leg_y + 22}" r="5.5" fill="#2563eb" fill-opacity="0.82"/>')
    svg.append(f'  <text x="{leg_x + 36}" y="{leg_y + 26}" font-family="Helvetica, Arial, sans-serif" font-size="13" fill="#1e293b">Season 1 (n=160)</text>')
    svg.append(f'  <circle cx="{leg_x + 20}" cy="{leg_y + 46}" r="5.5" fill="#16a34a" fill-opacity="0.82"/>')
    svg.append(f'  <text x="{leg_x + 36}" y="{leg_y + 50}" font-family="Helvetica, Arial, sans-serif" font-size="13" fill="#1e293b">Season 2 (n=160)</text>')
    svg.append(f'  <circle cx="{leg_x + 20}" cy="{leg_y + 68}" r="3.0" fill="#cbd5e1"/>')
    svg.append(f'  <text x="{leg_x + 36}" y="{leg_y + 72}" font-family="Helvetica, Arial, sans-serif" font-size="12" fill="#64748b">Floor Projection</text>')

    svg.append("</svg>")
    _svg_to_png(svg, output_png)


def print_pca_summary(name: str, pca_res: dict, feature_cols: list[str], is_corr: bool):
    print(f"\n{'=' * 78}")
    print(f"PCA SUMMARY: {name}")
    print(f"{'=' * 78}")

    eigenvalues = pca_res["eigenvalues"]
    explained = pca_res["explained_pct"]
    cum_pct = pca_res["cumulative_pct"]
    eigenvectors = pca_res["eigenvectors"]
    loadings = pca_res["factor_loadings"]
    p = len(feature_cols)

    print("\n1. Eigenvalues Table:")
    print(f"{'PC':<6} {'Eigenvalue':>14} {'Variance (%)':>16} {'Cumulative (%)':>18}")
    print("-" * 58)
    for k in range(p):
        print(f"PC{k + 1:<4} {eigenvalues[k]:>14.4f} {explained[k]:>15.2f}% {cum_pct[k]:>17.2f}%")

    for thresh in (90.0, 95.0, 99.0):
        num_pcs = next(k + 1 for k, c in enumerate(cum_pct) if c >= thresh)
        print(f"-> PCs needed for >= {thresh:.0f}% variability: {num_pcs} (Cumulative = {cum_pct[num_pcs - 1]:.2f}%)")

    print(f"-> First 2 PCs explain: {cum_pct[1]:.2f}% of total variability")
    print(f"-> First 3 PCs explain: {cum_pct[2]:.2f}% of total variability")

    print("\n2. Eigenvectors (Component Weights / Linear Combination Coefficients):")
    header = f"{'Variable':<10}" + "".join(f"{'PC' + str(k + 1):>10}" for k in range(p))
    print(header)
    print("-" * len(header))
    for i, col in enumerate(feature_cols):
        row_str = f"{col:<10}" + "".join(f"{eigenvectors[k][i]:>10.4f}" for k in range(p))
        print(row_str)

    print("\n3. Factor Loadings (Correlations between Original Variables and PCs):")
    print(header)
    print("-" * len(header))
    for i, col in enumerate(feature_cols):
        row_str = f"{col:<10}" + "".join(f"{loadings[i][k]:>10.4f}" for k in range(p))
        print(row_str)

    print("\n4. Mathematical Equations of First Three Principal Components:")
    var_suffix = "*" if is_corr else "'"
    for k in range(3):
        terms = []
        for i, col in enumerate(feature_cols):
            coeff = eigenvectors[k][i]
            sign = "+" if coeff >= 0 else "-"
            if i == 0:
                terms.append(f"{coeff:.4f}·{col}{var_suffix}")
            else:
                terms.append(f" {sign} {abs(coeff):.4f}·{col}{var_suffix}")
        print(f"  PC{k + 1} = " + "".join(terms))


def main():
    rows, seasons, _ = load_data(DATA_PATH, FEATURE_COLS)
    means, stds, cov_mat, corr_mat, centered_rows, standardized_rows = compute_stats_and_matrices(
        rows, FEATURE_COLS
    )

    print("Feature Means and Standard Deviations:")
    for col, m, s in zip(FEATURE_COLS, means, stds):
        print(f"  {col:<5}: mean = {m:8.4f}, std = {s:8.4f}")

    # Run PCA on Covariance Matrix and Correlation Matrix
    pca_cov = perform_pca(cov_mat, centered_rows, stds, is_corr=False)
    pca_corr = perform_pca(corr_mat, standardized_rows, stds, is_corr=True)

    print_pca_summary("A) PCA USING COVARIANCE MATRIX", pca_cov, FEATURE_COLS, is_corr=False)
    print_pca_summary("B) PCA USING CORRELATION MATRIX", pca_corr, FEATURE_COLS, is_corr=True)

    # Save all required plots
    save_heatmap_image(
        cov_mat,
        FEATURE_COLS,
        title="Sample Covariance Matrix",
        subtitle="Soybean Cultivars Features (Target 'GY' Omitted, Log-Scaled Colors)",
        output_png=os.path.join(BASE_DIR, "covariance_matrix.png"),
        is_corr=False,
    )
    save_heatmap_image(
        corr_mat,
        FEATURE_COLS,
        title="Pearson Correlation Matrix",
        subtitle="Soybean Cultivars Features (Target 'GY' Omitted)",
        output_png=os.path.join(BASE_DIR, "correlation_matrix.png"),
        is_corr=True,
    )

    # Covariance PCA plots
    save_scree_plot(
        pca_cov,
        title="Scree Plot - PCA on Covariance Matrix",
        subtitle="Individual and Cumulative Explained Variance across Principal Components",
        output_png=os.path.join(BASE_DIR, "cov_scree_plot.png"),
    )
    save_2d_pca_plot(
        pca_cov,
        seasons,
        FEATURE_COLS,
        title="2D PCA Projection (PC1 vs PC2) - Covariance Matrix",
        subtitle=f"First Two PCs Explain {pca_cov['cumulative_pct'][1]:.2f}% of Total Variance",
        output_png=os.path.join(BASE_DIR, "cov_pca_2d.png"),
    )
    save_3d_pca_plot(
        pca_cov,
        seasons,
        title="3D PCA Projection (PC1, PC2, PC3) - Covariance Matrix",
        subtitle=f"First Three PCs Explain {pca_cov['cumulative_pct'][2]:.2f}% of Total Variance",
        output_png=os.path.join(BASE_DIR, "cov_pca_3d.png"),
    )

    # Correlation PCA plots
    save_scree_plot(
        pca_corr,
        title="Scree Plot - PCA on Correlation Matrix",
        subtitle="Individual and Cumulative Explained Variance across Principal Components",
        output_png=os.path.join(BASE_DIR, "corr_scree_plot.png"),
    )
    save_2d_pca_plot(
        pca_corr,
        seasons,
        FEATURE_COLS,
        title="2D PCA Projection (PC1 vs PC2) - Correlation Matrix",
        subtitle=f"First Two PCs Explain {pca_corr['cumulative_pct'][1]:.2f}% of Total Variance",
        output_png=os.path.join(BASE_DIR, "corr_pca_2d.png"),
    )
    save_3d_pca_plot(
        pca_corr,
        seasons,
        title="3D PCA Projection (PC1, PC2, PC3) - Correlation Matrix",
        subtitle=f"First Three PCs Explain {pca_corr['cumulative_pct'][2]:.2f}% of Total Variance",
        output_png=os.path.join(BASE_DIR, "corr_pca_3d.png"),
    )

    print("\nAll 8 plots generated successfully in workspace directory.")


if __name__ == "__main__":
    main()

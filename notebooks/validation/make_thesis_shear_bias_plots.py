#!/usr/bin/env python
"""Create thesis-ready BATSim-GalSim shear-bias comparison plots.

This script mirrors the relevant analysis in ``new_shear_accuracy.ipynb``:
it loads the matched FPFS auto-detection catalogs, bins galaxies by COSMOS
magnitude, single-Sérsic index, and half-light radius, then plots
BATSim-GalSim differences in multiplicative and additive bias.
"""

from __future__ import annotations

import argparse
import os
import pickle
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import galsim
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm


SCALE = 0.2
NN = 96
LENS_SHEAR = 0.02
N_ROT = 4
N_GAL_STAMP = 1600
TILES_PER_SIDE = 40
PSF_M20_NORM = -0.3192

MAG_EDGES = np.array([18, 20, 21, 22, 23, 24, 25, 25.3], dtype=float)
SERSIC_EDGES = np.array([0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0], dtype=float)
HLR_EDGES = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.2, 2.0], dtype=float)

IBM_BLUE = "#648FFF"
IBM_MAGENTA = "#DC267F"

POS_Y = "y"
POS_X = "x"


def apply_plot_style() -> None:
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["cmr10"]
    plt.rcParams["mathtext.fontset"] = "cm"
    plt.rcParams["figure.facecolor"] = "white"
    plt.rc("axes", unicode_minus=False)
    plt.rc("axes.formatter", use_mathtext=True)


def detection_to_galaxy_k(cat: np.ndarray, nn: int, n_rot: int, tiles_per_side: int) -> np.ndarray:
    tile_size = int(nn * np.sqrt(n_rot))
    y = cat[POS_Y].astype(int)
    x = cat[POS_X].astype(int)
    return (y // tile_size) * tiles_per_side + (x // tile_size)


def detection_to_rot_k(cat: np.ndarray, nn: int, rot_per_side: int) -> np.ndarray:
    y = cat[POS_Y].astype(int)
    x = cat[POS_X].astype(int)
    return (y // nn) * rot_per_side + (x // nn)


def keep_one_per_tile(cat: np.ndarray, k_array: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    _, first_idx = np.unique(k_array, return_index=True)
    first_idx = np.sort(first_idx)
    return cat[first_idx], k_array[first_idx]


def measure_shear(
    gal_array: np.ndarray,
    psf_array: np.ndarray,
    npix: int,
    pixel_scale: float,
    noise_variance: float,
    noise_array: np.ndarray | None,
    detection: np.ndarray | None,
) -> np.ndarray:
    import anacal

    fpfs_config = anacal.fpfs.FpfsConfig(
        sigma_shapelets=0.52,
        sigma_shapelets2=0.45,
        npix=npix,
    )
    return anacal.fpfs.process_image(
        fpfs_config=fpfs_config,
        mag_zero=30.0,
        gal_array=gal_array,
        psf_array=psf_array,
        pixel_scale=pixel_scale,
        noise_variance=max(noise_variance, 0.23),
        noise_array=noise_array,
        detection=detection,
    )


def jackknife_over_tiles(E_tiles: np.ndarray, R_tiles: np.ndarray, lens_shear: float) -> dict:
    E_tiles = np.asarray(E_tiles)
    R_tiles = np.asarray(R_tiles)
    E_tot = E_tiles.sum()
    R_tot = R_tiles.sum()
    g1 = E_tot / R_tot
    m = g1 / lens_shear - 1
    n = len(E_tiles)
    with np.errstate(divide="ignore", invalid="ignore"):
        g1_loo = (E_tot - E_tiles) / (R_tot - R_tiles)
    m_loo = g1_loo / lens_shear - 1
    g1_err = np.sqrt((n - 1) / n * np.sum((g1_loo - g1_loo.mean()) ** 2))
    m_err = np.sqrt((n - 1) / n * np.sum((m_loo - m_loo.mean()) ** 2))
    return {"g1": g1, "g1_err": g1_err, "m": m, "m_err": m_err}


def jackknife_additive_over_tiles(E2_tiles: np.ndarray, R2_tiles: np.ndarray) -> dict:
    E2_tiles = np.asarray(E2_tiles)
    R2_tiles = np.asarray(R2_tiles)
    E2_tot = E2_tiles.sum()
    R2_tot = R2_tiles.sum()
    c = E2_tot / R2_tot
    n = len(E2_tiles)
    with np.errstate(divide="ignore", invalid="ignore"):
        c_loo = (E2_tot - E2_tiles) / (R2_tot - R2_tiles)
    c_err = np.sqrt((n - 1) / n * np.sum((c_loo - c_loo.mean()) ** 2))
    return {"c": c, "c_err": c_err}


def load_cosmos_properties(cosmos_catalog_dir: Path | None) -> dict:
    np.random.seed(14)
    if cosmos_catalog_dir is None:
        default_catalog_dir = Path("/home/b7009348/share/COSMOS_25.2_training_sample")
        cosmos_catalog_dir = default_catalog_dir if default_catalog_dir.exists() else None

    if cosmos_catalog_dir is None:
        cosmos = galsim.COSMOSCatalog()
    else:
        cosmos = galsim.COSMOSCatalog(dir=str(cosmos_catalog_dir))
    all_inds = np.arange(len(cosmos))
    stamp_inds = np.array_split(all_inds, int(np.ceil(len(cosmos) / N_GAL_STAMP)))
    all_gal_inds = np.concatenate(stamp_inds)
    records = cosmos.getParametricRecord(index=all_gal_inds)
    stamp_offsets = np.concatenate([[0], np.cumsum([len(s) for s in stamp_inds])])

    return {
        "cosmos": cosmos,
        "stamp_inds": stamp_inds,
        "stamp_offsets": stamp_offsets,
        "all_gal_inds": all_gal_inds,
        "mag": records["mag_auto"],
        "hlr": records["hlr"][:, 0],
        "sersic_n": records["sersicfit"][:, 2],
        "is_bulgefit": records["use_bulgefit"].astype(bool),
    }


def load_matched_catalogs(
    cosmos_image_path: Path,
    stamp_inds: list[np.ndarray],
    force_remeasure: bool,
) -> tuple[list[np.ndarray], list[np.ndarray], list[np.ndarray]]:
    stamp_size = int(NN * np.sqrt(N_ROT * N_GAL_STAMP))
    rot_per_side = TILES_PER_SIDE * int(np.sqrt(N_ROT))
    psf = galsim.Moffat(beta=3.5, fwhm=0.7)
    psf_array = psf.drawImage(nx=NN, ny=NN, scale=SCALE).array

    gs_catalogs = []
    bt_catalogs = []
    cat_gal_k = []

    for i, inds in enumerate(tqdm(stamp_inds, desc="matched catalogs")):
        cache_file = cosmos_image_path / f"fpfs_autodet_cat_v6_{i}.pkl"
        n_gal = len(inds)

        if cache_file.exists() and not force_remeasure:
            with cache_file.open("rb") as f:
                gs_cat, bt_cat, k, *_ = pickle.load(f)
        else:
            gs_stamp = galsim.fits.read(
                file_name=cosmos_image_path / f"galsim_stamp_{i}_{stamp_size}x{stamp_size}.fits"
            )
            bt_stamp = galsim.fits.read(
                file_name=cosmos_image_path / f"batsim_stamp_{i}_{stamp_size}x{stamp_size}.fits"
            )
            gs_cat = measure_shear(
                gs_stamp.array,
                psf_array,
                NN,
                SCALE,
                noise_variance=0.23,
                noise_array=None,
                detection=None,
            )
            bt_cat = measure_shear(
                bt_stamp.array,
                psf_array,
                NN,
                SCALE,
                noise_variance=0.23,
                noise_array=None,
                detection=None,
            )

            k_rot_gs = detection_to_rot_k(gs_cat, NN, rot_per_side)
            k_rot_bt = detection_to_rot_k(bt_cat, NN, rot_per_side)
            gs_cat, k_rot_gs = keep_one_per_tile(gs_cat, k_rot_gs)
            bt_cat, k_rot_bt = keep_one_per_tile(bt_cat, k_rot_bt)

            k_gs = detection_to_galaxy_k(gs_cat, NN, N_ROT, TILES_PER_SIDE)
            k_bt = detection_to_galaxy_k(bt_cat, NN, N_ROT, TILES_PER_SIDE)

            valid_gs = k_gs < n_gal
            gs_cat, k_rot_gs, k_gs = gs_cat[valid_gs], k_rot_gs[valid_gs], k_gs[valid_gs]
            valid_bt = k_bt < n_gal
            bt_cat, k_rot_bt, k_bt = bt_cat[valid_bt], k_rot_bt[valid_bt], k_bt[valid_bt]

            common_rot = np.array(sorted(set(k_rot_gs) & set(k_rot_bt)))
            gs_mask = np.isin(k_rot_gs, common_rot)
            bt_mask = np.isin(k_rot_bt, common_rot)
            gs_ord = np.argsort(k_rot_gs[gs_mask])
            bt_ord = np.argsort(k_rot_bt[bt_mask])
            gs_cat = gs_cat[gs_mask][gs_ord]
            bt_cat = bt_cat[bt_mask][bt_ord]
            k = k_gs[gs_mask][gs_ord]

            with cache_file.open("wb") as f:
                pickle.dump((gs_cat, bt_cat, k, np.array([]), np.array([]), np.array([])), f)

        gs_catalogs.append(gs_cat)
        bt_catalogs.append(bt_cat)
        cat_gal_k.append(k)

    return gs_catalogs, bt_catalogs, cat_gal_k


def compute_trace_mask(
    bt_catalogs: list[np.ndarray],
    cat_gal_k: list[np.ndarray],
    stamp_offsets: np.ndarray,
    n_galaxies: int,
    threshold: float,
    psf_m20_norm: float,
) -> tuple[np.ndarray, np.ndarray]:
    mean_m20 = np.zeros(n_galaxies)
    mean_m00 = np.zeros(n_galaxies)
    n_obs = np.zeros(n_galaxies, dtype=int)

    for i, (cat, k_arr) in enumerate(zip(bt_catalogs, cat_gal_k)):
        g_idx = stamp_offsets[i] + k_arr
        np.add.at(mean_m20, g_idx, cat["fpfs_m20"])
        np.add.at(mean_m00, g_idx, cat["fpfs_m00"])
        np.add.at(n_obs, g_idx, 1)

    has_obs = n_obs > 0
    mean_m20[has_obs] /= n_obs[has_obs]
    mean_m00[has_obs] /= n_obs[has_obs]

    resolution = np.full(n_galaxies, np.nan)
    resolution[has_obs] = psf_m20_norm / (mean_m20[has_obs] / mean_m00[has_obs])

    trace_mask = np.zeros(n_galaxies, dtype=bool)
    trace_mask[has_obs] = resolution[has_obs] > threshold
    return trace_mask, resolution


def compute_bias_for_mask(
    galaxy_mask: np.ndarray,
    gs_catalogs: list[np.ndarray],
    bt_catalogs: list[np.ndarray],
    cat_gal_k: list[np.ndarray],
    stamp_offsets: np.ndarray,
) -> tuple[dict | None, dict | None]:
    E_gs, R_gs, E_bt, R_bt = [], [], [], []
    for i in range(len(gs_catalogs)):
        k_map = cat_gal_k[i]
        det_mask = galaxy_mask[stamp_offsets[i] + k_map]
        for E_list, R_list, cat in ((E_gs, R_gs, gs_catalogs[i]), (E_bt, R_bt, bt_catalogs[i])):
            E_list.append((cat["fpfs_w"] * cat["fpfs_e1"])[det_mask].sum())
            R_list.append(
                (
                    cat["fpfs_dw_dg1"] * cat["fpfs_e1"]
                    + cat["fpfs_w"] * cat["fpfs_de1_dg1"]
                )[det_mask].sum()
            )
    if np.sum(R_gs) == 0 or np.sum(R_bt) == 0:
        return None, None
    return jackknife_over_tiles(E_gs, R_gs, LENS_SHEAR), jackknife_over_tiles(E_bt, R_bt, LENS_SHEAR)


def compute_additive_bias_for_mask(
    galaxy_mask: np.ndarray,
    gs_catalogs: list[np.ndarray],
    bt_catalogs: list[np.ndarray],
    cat_gal_k: list[np.ndarray],
    stamp_offsets: np.ndarray,
) -> tuple[dict | None, dict | None]:
    E2_gs, R2_gs, E2_bt, R2_bt = [], [], [], []
    for i in range(len(gs_catalogs)):
        k_map = cat_gal_k[i]
        det_mask = galaxy_mask[stamp_offsets[i] + k_map]
        for E_list, R_list, cat in ((E2_gs, R2_gs, gs_catalogs[i]), (E2_bt, R2_bt, bt_catalogs[i])):
            E_list.append((cat["fpfs_w"] * cat["fpfs_e2"])[det_mask].sum())
            R_list.append(
                (
                    cat["fpfs_dw_dg2"] * cat["fpfs_e2"]
                    + cat["fpfs_w"] * cat["fpfs_de2_dg2"]
                )[det_mask].sum()
            )
    if np.sum(R2_gs) == 0 or np.sum(R2_bt) == 0:
        return None, None
    return jackknife_additive_over_tiles(E2_gs, R2_gs), jackknife_additive_over_tiles(E2_bt, R2_bt)


def run_binned_analysis(
    prop: np.ndarray,
    bin_edges: np.ndarray,
    base_mask: np.ndarray | None,
    min_bin_galaxies: int,
    kind: str,
    gs_catalogs: list[np.ndarray],
    bt_catalogs: list[np.ndarray],
    cat_gal_k: list[np.ndarray],
    stamp_offsets: np.ndarray,
    label: str,
) -> dict:
    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    n_bins = len(centers)
    out = {
        "centers": centers,
        "bin_edges": bin_edges,
        "n_gals": np.zeros(n_bins, dtype=int),
        f"{kind}_gs": np.full(n_bins, np.nan),
        f"{kind}_gs_err": np.full(n_bins, np.nan),
        f"{kind}_bt": np.full(n_bins, np.nan),
        f"{kind}_bt_err": np.full(n_bins, np.nan),
    }

    for b in tqdm(range(n_bins), desc=label):
        mask = (prop >= bin_edges[b]) & (prop < bin_edges[b + 1])
        if base_mask is not None:
            mask &= base_mask
        out["n_gals"][b] = mask.sum()
        if out["n_gals"][b] < min_bin_galaxies:
            continue

        if kind == "m":
            res_gs, res_bt = compute_bias_for_mask(mask, gs_catalogs, bt_catalogs, cat_gal_k, stamp_offsets)
        elif kind == "c":
            res_gs, res_bt = compute_additive_bias_for_mask(mask, gs_catalogs, bt_catalogs, cat_gal_k, stamp_offsets)
        else:
            raise ValueError("kind must be 'm' or 'c'")

        if res_gs is None:
            continue
        out[f"{kind}_gs"][b] = res_gs[kind]
        out[f"{kind}_gs_err"][b] = res_gs[f"{kind}_err"]
        out[f"{kind}_bt"][b] = res_bt[kind]
        out[f"{kind}_bt_err"][b] = res_bt[f"{kind}_err"]

    out[f"d{kind}"] = out[f"{kind}_bt"] - out[f"{kind}_gs"]
    out[f"d{kind}_err"] = np.sqrt(out[f"{kind}_gs_err"] ** 2 + out[f"{kind}_bt_err"] ** 2)
    return out


def get_valid(res: dict, ref_err: float, kind: str, max_err_factor: float) -> np.ndarray:
    return (
        ~np.isnan(res[f"d{kind}"])
        & (res[f"{kind}_gs_err"] < max_err_factor * ref_err)
        & (res[f"{kind}_bt_err"] < max_err_factor * ref_err)
    )


def build_sample_results(
    sample_mask: np.ndarray,
    properties: dict,
    gs_catalogs: list[np.ndarray],
    bt_catalogs: list[np.ndarray],
    cat_gal_k: list[np.ndarray],
    min_bin_galaxies: int,
) -> dict:
    stamp_offsets = properties["stamp_offsets"]
    single_sersic_mask = sample_mask & ~properties["is_bulgefit"]

    m_full_gs, m_full_bt = compute_bias_for_mask(
        sample_mask, gs_catalogs, bt_catalogs, cat_gal_k, stamp_offsets
    )
    c_full_gs, c_full_bt = compute_additive_bias_for_mask(
        sample_mask, gs_catalogs, bt_catalogs, cat_gal_k, stamp_offsets
    )

    return {
        "m_full_gs": m_full_gs,
        "m_full_bt": m_full_bt,
        "c_full_gs": c_full_gs,
        "c_full_bt": c_full_bt,
        "m_results": [
            run_binned_analysis(
                properties["mag"],
                MAG_EDGES,
                sample_mask,
                min_bin_galaxies,
                "m",
                gs_catalogs,
                bt_catalogs,
                cat_gal_k,
                stamp_offsets,
                "magnitude / m",
            ),
            run_binned_analysis(
                properties["sersic_n"],
                SERSIC_EDGES,
                single_sersic_mask,
                min_bin_galaxies,
                "m",
                gs_catalogs,
                bt_catalogs,
                cat_gal_k,
                stamp_offsets,
                "Sérsic n / m",
            ),
            run_binned_analysis(
                properties["hlr"],
                HLR_EDGES,
                sample_mask,
                min_bin_galaxies,
                "m",
                gs_catalogs,
                bt_catalogs,
                cat_gal_k,
                stamp_offsets,
                "half-light radius / m",
            ),
        ],
        "c_results": [
            run_binned_analysis(
                properties["mag"],
                MAG_EDGES,
                sample_mask,
                min_bin_galaxies,
                "c",
                gs_catalogs,
                bt_catalogs,
                cat_gal_k,
                stamp_offsets,
                "magnitude / c",
            ),
            run_binned_analysis(
                properties["sersic_n"],
                SERSIC_EDGES,
                single_sersic_mask,
                min_bin_galaxies,
                "c",
                gs_catalogs,
                bt_catalogs,
                cat_gal_k,
                stamp_offsets,
                "Sérsic n / c",
            ),
            run_binned_analysis(
                properties["hlr"],
                HLR_EDGES,
                sample_mask,
                min_bin_galaxies,
                "c",
                gs_catalogs,
                bt_catalogs,
                cat_gal_k,
                stamp_offsets,
                "half-light radius / c",
            ),
        ],
    }


def symmetric_limits(results: list[dict], kind: str, valid_masks: list[np.ndarray]) -> tuple[float, float]:
    vals = []
    for res, valid in zip(results, valid_masks):
        y = res[f"d{kind}"][valid] * 1e3
        yerr = res[f"d{kind}_err"][valid] * 1e3
        if y.size:
            vals.extend((y - yerr).tolist())
            vals.extend((y + yerr).tolist())
    if not vals:
        return -1, 1
    vmax = max(abs(np.nanmin(vals)), abs(np.nanmax(vals)))
    vmax = 1.1 * vmax if vmax > 0 else 1
    return -vmax, vmax


def plot_bias_difference_grid(
    sample: dict,
    title: str,
    output_base: Path,
    max_err_factor: float,
) -> None:
    titles = ["Magnitude", r"S$\acute{\mathrm{e}}$rsic n", "Half-light radius"]
    xlabels = [r"$m_{\rm auto}$", r"$n$", r"$r_e$ [arcsec]"]
    log_x = [False, False, True]

    m_valid = [
        get_valid(res, sample["m_full_gs"]["m_err"], "m", max_err_factor)
        for res in sample["m_results"]
    ]
    c_valid = [
        get_valid(res, sample["c_full_gs"]["c_err"], "c", max_err_factor)
        for res in sample["c_results"]
    ]

    m_ylim = symmetric_limits(sample["m_results"], "m", m_valid)
    c_ylim = symmetric_limits(sample["c_results"], "c", c_valid)

    fig, axes = plt.subplots(2, 3, figsize=(7.4, 4.8), sharey="row", constrained_layout=True)

    for col, (panel_title, xlabel, use_log) in enumerate(zip(titles, xlabels, log_x)):
        for row, kind, results, valid_masks, ylim in (
            (0, "m", sample["m_results"], m_valid, m_ylim),
            (1, "c", sample["c_results"], c_valid, c_ylim),
        ):
            ax = axes[row, col]
            res = results[col]
            valid = valid_masks[col]
            x = res["centers"][valid]
            y = res[f"d{kind}"][valid] * 1e3
            yerr = res[f"d{kind}_err"][valid] * 1e3
            colour = IBM_BLUE if kind == "m" else IBM_MAGENTA

            ax.axhline(0, color="0.2", lw=0.8, ls="--", zorder=1)
            ax.errorbar(
                x,
                y,
                yerr=yerr,
                fmt="o-",
                ms=3.8,
                lw=1.0,
                elinewidth=0.8,
                capsize=2.5,
                color=colour,
                mfc=colour,
                mec=colour,
                zorder=2,
            )
            ax.set_ylim(*ylim)
            ax.tick_params(direction="in", top=True, right=True, labelsize=10)
            if use_log:
                ax.set_xscale("log")
            if row == 0:
                ax.set_title(panel_title, fontsize=12)
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel(xlabel, fontsize=11)
            if col == 0 and kind == "m":
                ax.set_ylabel(r"$\Delta m$ [$10^{-3}$]", fontsize=11)
            if col == 0 and kind == "c":
                ax.set_ylabel(r"$\Delta c$ [$10^{-3}$]", fontsize=11)

    for suffix in (".pdf", ".png"):
        save_path = output_base.with_suffix(suffix)
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Saved {save_path}")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cosmos-image-path",
        type=Path,
        default=repo_root / "validation" / "COSMOS_lensed",
        help="Directory containing COSMOS_lensed FITS stamps and FPFS catalog caches.",
    )
    parser.add_argument(
        "--cosmos-catalog-dir",
        type=Path,
        default=None,
        help=(
            "Directory containing the GalSim COSMOS_25.2_training_sample catalog. "
            "Defaults to /home/b7009348/share/COSMOS_25.2_training_sample when present."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=script_dir,
        help="Directory for output figures.",
    )
    parser.add_argument(
        "--trace-threshold",
        type=float,
        default=0.45,
        help="Resolution/trace-radius threshold used for the trace cut.",
    )
    parser.add_argument(
        "--psf-m20-norm",
        type=float,
        default=PSF_M20_NORM,
        help="PSF fpfs_m20/fpfs_m00 normalisation used by the notebook.",
    )
    parser.add_argument(
        "--min-bin-galaxies",
        type=int,
        default=100,
        help="Minimum galaxy count required before estimating a binned bias.",
    )
    parser.add_argument(
        "--max-error-factor",
        type=float,
        default=5.0,
        help="Hide bins with errors larger than this factor times the full-sample error.",
    )
    parser.add_argument(
        "--force-remeasure",
        action="store_true",
        help="Ignore cached FPFS catalogs and remeasure from the large FITS stamps.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    apply_plot_style()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    properties = load_cosmos_properties(args.cosmos_catalog_dir)
    gs_catalogs, bt_catalogs, cat_gal_k = load_matched_catalogs(
        args.cosmos_image_path,
        properties["stamp_inds"],
        args.force_remeasure,
    )

    n_galaxies = len(properties["all_gal_inds"])
    full_mask = np.ones(n_galaxies, dtype=bool)
    trace_mask, _ = compute_trace_mask(
        bt_catalogs,
        cat_gal_k,
        properties["stamp_offsets"],
        n_galaxies,
        args.trace_threshold,
        args.psf_m20_norm,
    )

    print(f"Full sample: {full_mask.sum()} galaxies")
    print(f"Trace cut:   {trace_mask.sum()} galaxies pass R > {args.trace_threshold:.2f}")

    full_sample = build_sample_results(
        full_mask,
        properties,
        gs_catalogs,
        bt_catalogs,
        cat_gal_k,
        args.min_bin_galaxies,
    )
    trace_sample = build_sample_results(
        trace_mask,
        properties,
        gs_catalogs,
        bt_catalogs,
        cat_gal_k,
        args.min_bin_galaxies,
    )

    plot_bias_difference_grid(
        full_sample,
        "BATSim - GalSim shear bias differences: full sample",
        args.output_dir / "thesis_shear_bias_differences_full_sample",
        args.max_error_factor,
    )
    plot_bias_difference_grid(
        trace_sample,
        f"BATSim - GalSim shear bias differences: trace cut R > {args.trace_threshold:.2f}",
        args.output_dir / "thesis_shear_bias_differences_trace_cut",
        args.max_error_factor,
    )


if __name__ == "__main__":
    main()

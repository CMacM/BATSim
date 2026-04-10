"""Benchmarks to check computation time and memory usage for batsim."""
import batsim.stamp as batstamp
import batsim.transforms as batforms
import batsim
import contextlib
import galsim
import io
import numpy as np
import time

def time_shear_speed(nn=64, scale=0.2):
    
    # create galaxy
    gal = galsim.Sersic(n=1.5, half_light_radius=1.5, flux=40)
    
    # start timing for batsim
    start = time.time()
    
    # generate stamp, shear and sample
    gal_stamp = batstamp.Stamp(nn=nn, scale=scale, centering='galsim')
    lens = batforms.LensTransform(gamma1=0.2, gamma2=0, kappa=0,
                                  center=[-0.5*0.2, -0.5*0.2])
    gal_stamp.transform_grids(lens)
    bat_array = gal_stamp.sample_galaxy(gal)
    
    end = time.time()
    
    bat_time = end-start
    
    # start timing for galsim
    start = time.time()
    
    gal_shear = gal.shear(g1=0.2)
    gal_array = gal_shear.drawImage(nx=nn, ny=nn, scale=scale).array
    
    end = time.time()
    
    gal_time = end-start
    
    return {'batsim time': bat_time, 'galsim time' : gal_time}

def time_ia_speed(nn=128, scale=0.1):
    '''
        Times the time taken to apply an IA shear to a 128x128 pixel
        image and compares it with the time required to apply an
        affine shear.
    '''
    # initialise a galaxy object
    gal = galsim.Sersic(n=1.5, half_light_radius=1.5, flux=40)
    
    # start timing for ia transform
    ia_start = time.time()
    
    # crete stamp, apply transform, and sample gal obj
    ia_stamp = batstamp.Stamp(nn=nn, scale=scale)
    ia = batforms.IaTransform(scale=scale, hlr=1.5, phi=0.2)
    ia_stamp.transform_grids(ia)
    ia_gal = ia_stamp.sample_galaxy(gal)
    
    # stop timing
    ia_end = time.time()
    ia_time = ia_end - ia_start
    
    # get equivalent shear at hlr to pass to lens transform
    g1, g2 = ia.get_g1g2(1.5,0)
    
    # start timing for lens
    aff_start = time.time()
    
    # crete stamp, apply transform, and sample gal obj
    lens_stamp = batstamp.Stamp(nn=nn, scale=scale)
    lens = batforms.LensTransform(gamma1=g1, gamma2=g2, kappa=0)
    lens_stamp.transform_grids(lens)
    lens_gal = lens_stamp.sample_galaxy(gal)
    
    # stop timing
    aff_end = time.time()
    aff_time = aff_end - aff_start
    
    return {'IA time' : ia_time, 'Lens time' : aff_time}


def _parse_simulate_profile_logs(log_lines):
    stats = {}
    timings = {}
    for line in log_lines:
        msg = line.split("] ", 1)[-1]
        if msg.startswith("stats "):
            for token in msg[6:].split():
                if "=" not in token:
                    continue
                key, value = token.split("=", 1)
                try:
                    stats[key] = int(value)
                except ValueError:
                    try:
                        stats[key] = float(value)
                    except ValueError:
                        stats[key] = value
        elif "=" in msg:
            key, value = msg.split("=", 1)
            value = value[:-1] if value.endswith("s") else value
            try:
                timings[key] = float(value)
            except ValueError:
                timings[key] = value
    return {"timings": timings, "stats": stats}


def _extract_parametric_profile_info(cosmos_catalog, catalog_index, gal_obj):
    info = {
        "catalog_index": int(catalog_index),
        "gsobject_type": type(gal_obj).__name__,
    }
    for attr in ("flux", "nyquist_scale"):
        if hasattr(gal_obj, attr):
            try:
                info[attr] = float(getattr(gal_obj, attr))
            except Exception:
                pass
    param_cat = getattr(cosmos_catalog, "param_cat", None)
    if param_cat is None:
        return info
    keys = []
    if hasattr(param_cat, "colnames"):
        keys = ["mag_auto", "flux_radius", "zphot"] + [k for k in ("use_bulgefit", "viable_sersic") if k in param_cat.colnames]
    elif hasattr(param_cat, "dtype") and param_cat.dtype.names:
        keys = [k for k in ("mag_auto", "flux_radius", "zphot", "use_bulgefit", "viable_sersic") if k in param_cat.dtype.names]
    if not keys:
        return info
    row = param_cat[int(catalog_index)]
    for key in keys:
        try:
            value = row[key]
            if hasattr(value, "item"):
                value = value.item()
            info[key] = value
        except Exception:
            pass
    return info


def benchmark_parametric_cosmos_profiles(
    n_galaxies=5,
    ngrid=128,
    pix_scale=0.2,
    psf_obj=None,
    draw_method="auto",
    truncate_ratio=1.0,
    maximum_num_grids=4096,
    force_ngrid=False,
    seed=1234,
    cosmos_catalog=None,
):
    """Run a lightweight per-galaxy benchmark using parametric COSMOS profiles.

    Returns a list of dictionaries containing profile metadata, parsed
    `simulate_galaxy(profile=True)` logs, and end-to-end elapsed time.
    """
    cosmos_catalog = cosmos_catalog or galsim.COSMOSCatalog()
    rng = np.random.RandomState(seed)
    indices = rng.choice(len(cosmos_catalog), size=n_galaxies, replace=(n_galaxies > len(cosmos_catalog)))

    records = []
    for i, idx in enumerate(indices):
        gal = cosmos_catalog.makeGalaxy(index=int(idx), gal_type="parametric")
        profile_info = _extract_parametric_profile_info(cosmos_catalog, idx, gal)

        log_buf = io.StringIO()
        t0 = time.perf_counter()
        with contextlib.redirect_stdout(log_buf):
            image = batsim.simulate_galaxy(
                ngrid=ngrid,
                pix_scale=pix_scale,
                gal_obj=gal,
                psf_obj=psf_obj,
                truncate_ratio=truncate_ratio,
                maximum_num_grids=maximum_num_grids,
                draw_method=draw_method,
                force_ngrid=force_ngrid,
                profile=True,
            )
        elapsed_s = time.perf_counter() - t0

        profile_logs = [line for line in log_buf.getvalue().splitlines() if line.startswith("[simulate_galaxy]")]
        parsed_logs = _parse_simulate_profile_logs(profile_logs)
        record = {
            "galaxy_number": i,
            "profile": profile_info,
            "logger": parsed_logs,
            "elapsed_s": elapsed_s,
            "image_shape": tuple(image.shape),
            "image_sum": float(np.sum(image)),
        }
        records.append(record)

        print(
            f"[benchmark_parametric_cosmos_profiles] i={i} idx={int(idx)} "
            f"nn={parsed_logs['stats'].get('nn')} downsample_ratio={parsed_logs['stats'].get('downsample_ratio')} "
            f"elapsed_s={elapsed_s:.4e}"
        )
        print(f"[benchmark_parametric_cosmos_profiles] profile={profile_info}")
        for line in profile_logs:
            print(line)
    return records

if __name__ == "__main__":
    time_shear_speed()
    time_ia_speed()
    
    

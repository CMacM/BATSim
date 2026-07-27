"""Component-level BATSim benchmarks."""

import numpy as np

import batsim
from batsim.fft import (
    _apply_centering_phase,
    _apply_pixel_response,
    _compensate_block_integration,
    _extract_centered_coarse_image,
    _rfft_centered_image,
)
from batsim.sampling import _sample_galaxy_profile, _sample_psf_spectrum
from batsim.stamp import Stamp

try:
    from benchmarks import utils
except ImportError:
    import utils


class TimeGridSizing:
    """Benchmark supersampling and fine-grid decisions."""

    params = (
        ["typical", "compact", "elliptical"],
        [1, 2, 4],
    )
    param_names = ["profile", "integration_order"]

    def setup(self, profile, integration_order):
        if profile == "typical":
            self.galaxy = utils.sersic_galaxy()
        elif profile == "compact":
            self.galaxy = utils.compact_sersic_galaxy()
        elif profile == "elliptical":
            self.galaxy = utils.elliptical_sersic_galaxy()
        else:
            raise ValueError(f"Unknown profile: {profile}")

    def time_grid_decision(self, profile, integration_order):
        utils.grid_decision(self.galaxy, integration_order=integration_order)


class TimeStampConstruction:
    """Benchmark coordinate stamp construction on NumPy."""

    params = [256, 512, 1024]
    param_names = ["nn"]

    def time_stamp_construction(self, nn):
        Stamp(nn=nn, scale=utils.SCALE, backend=np, dtype=np.float64)


class TimeTransformApplication:
    """Benchmark public transforms over flattened coordinate arrays."""

    params = (
        ["lens", "ia", "flexion"],
        [256, 512],
    )
    param_names = ["transform", "nn"]

    def setup(self, transform, nn):
        self.coords = utils.coordinate_grid(nn, dtype=np.float64)
        self.transform = utils.transform_for_name(transform)

    def time_transform(self, transform, nn):
        self.transform.transform(self.coords)


class TimeCppSampling:
    """Benchmark GalSim profile sampling through BATSim's C++ bridge."""

    params = (
        ["gaussian", "sersic"],
        [128, 256],
    )
    param_names = ["profile", "nn"]

    def setup(self, profile, nn):
        self.galaxy = utils.gaussian_galaxy() if profile == "gaussian" else utils.sersic_galaxy()
        self.coords = utils.coordinate_grid(nn, dtype=np.float64)

    def time_sample_galaxy_profile(self, profile, nn):
        _sample_galaxy_profile(
            gal_obj=self.galaxy,
            coords=self.coords,
            fine_scale=utils.SCALE,
            real_dtype=np.float64,
        )


class TimeFFTGrid:
    """Benchmark FFT helper operations on square NumPy grids."""

    params = [512, 1024]
    param_names = ["n"]

    def setup(self, n):
        rng = np.random.default_rng(12345)
        self.image = rng.normal(size=(n, n)).astype(np.float32)
        self.spectrum = np.fft.rfft2(self.image).astype(np.complex64)

    def time_rfft_centered_image(self, n):
        _rfft_centered_image(
            xp=np,
            image=self.image,
            real_dtype=np.float32,
            complex_dtype=np.complex64,
            center_index=n // 2,
        )

    def time_apply_centering_phase(self, n):
        spectrum = self.spectrum.copy()
        _apply_centering_phase(np, spectrum, n, center_index=n // 2)

    def time_apply_pixel_response(self, n):
        spectrum = self.spectrum.copy()
        _apply_pixel_response(
            xp=np,
            spectrum=spectrum,
            n=n,
            fine_scale=utils.SCALE / 4,
            pix_scale=utils.SCALE,
        )

    def time_compensate_block_integration(self, n):
        spectrum = self.spectrum.copy()
        _compensate_block_integration(
            xp=np,
            spectrum=spectrum,
            n=n,
            fine_scale=utils.SCALE / 4,
            mode="quadrature",
            integration_order=2,
        )

    def time_extract_centered_coarse_image(self, n):
        _extract_centered_coarse_image(
            xp=np,
            image=self.image,
            downsample_ratio=4,
            center_index=n // 2,
        )


class TimePsfSpectrum:
    """Benchmark analytic PSF spectrum sampling through the C++ bridge."""

    params = (
        ["gaussian", "moffat"],
        [512, 1024],
    )
    param_names = ["psf", "n"]

    def setup(self, psf, n):
        self.psf = utils.gaussian_psf() if psf == "gaussian" else utils.moffat_psf()

    def time_sample_psf_spectrum(self, psf, n):
        _sample_psf_spectrum(
            psf_obj=self.psf,
            n=n,
            fine_scale=utils.SCALE / 4,
            complex_dtype=np.complex64,
        )


class TrackGridDecisions:
    """Track grid sizing decisions for representative profiles."""

    params = ["typical", "compact", "elliptical"]
    param_names = ["profile"]

    def setup(self, profile):
        if profile == "typical":
            self.galaxy = utils.sersic_galaxy()
        elif profile == "compact":
            self.galaxy = utils.compact_sersic_galaxy()
        elif profile == "elliptical":
            self.galaxy = utils.elliptical_sersic_galaxy()
        else:
            raise ValueError(f"Unknown profile: {profile}")

        self.decision = utils.grid_decision(self.galaxy, integration_order=2)

    def track_supersample(self, profile):
        return self.decision["supersample"]

    def track_fft_supersample(self, profile):
        return self.decision["fft_supersample"]

    def track_integration_order(self, profile):
        return self.decision["integration_order"]

    def track_fine_ngrid(self, profile):
        return self.decision["fine_ngrid"]

    def track_fine_compact(self, profile):
        return self.decision["fine_compact"]


class TrackRenderFlux:
    """Track flux residual for a small stable render."""

    def setup(self):
        self.kwargs = utils.render_kwargs(
            gal_obj=utils.gaussian_galaxy(),
            ngrid=64,
            psf_obj=utils.gaussian_psf(),
            backend="np",
        )

    def track_flux_residual(self):
        image = batsim.simulate_galaxy(**self.kwargs)
        return utils.flux_residual(image)

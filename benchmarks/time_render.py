"""End-to-end ``simulate_galaxy`` benchmarks."""

import batsim

try:
    from benchmarks import utils
except ImportError:
    import utils


class TimeRenderScenarios:
    """Benchmark representative public render paths."""

    params = (
        ["minimal", "pixel", "psf_kvalue", "psf_real", "lens", "ia", "flexion"],
        [64, 128, 256],
    )
    param_names = ["scenario", "ngrid"]
    timeout = 180

    def setup(self, scenario, ngrid):
        self.galaxy = utils.sersic_galaxy()
        self.psf = None
        self.transform = None
        self.draw_method = "auto"
        self.psf_mode = "kvalue"

        if scenario == "minimal":
            self.draw_method = "no_pixel"
        elif scenario == "pixel":
            pass
        elif scenario == "psf_kvalue":
            self.psf = utils.gaussian_psf()
        elif scenario == "psf_real":
            self.psf = utils.gaussian_psf()
            self.psf_mode = "real"
        elif scenario in ("lens", "ia", "flexion"):
            self.psf = utils.gaussian_psf()
            self.transform = utils.transform_for_name(scenario)
        else:
            raise ValueError(f"Unknown render scenario: {scenario}")

        self.kwargs = utils.render_kwargs(
            gal_obj=self.galaxy,
            ngrid=ngrid,
            transform_obj=self.transform,
            psf_obj=self.psf,
            draw_method=self.draw_method,
            psf_mode=self.psf_mode,
            backend="np",
        )

    def teardown(self, scenario, ngrid):
        utils.cleanup_backend("np")

    def time_render(self, scenario, ngrid):
        batsim.simulate_galaxy(**self.kwargs)


class TimeIntegrationScaling:
    """Benchmark block-integration order in a controlled medium render."""

    params = [1, 2, 4]
    param_names = ["integration_order"]
    timeout = 180

    def setup(self, integration_order):
        self.kwargs = utils.render_kwargs(
            gal_obj=utils.compact_sersic_galaxy(),
            ngrid=128,
            psf_obj=utils.gaussian_psf(),
            integration_order=integration_order,
            backend="np",
        )

    def teardown(self, integration_order):
        utils.cleanup_backend("np")

    def time_render_integration_order(self, integration_order):
        batsim.simulate_galaxy(**self.kwargs)


class TimeFluxOption:
    """Benchmark the cost of final input-flux normalization."""

    params = [True, False]
    param_names = ["force_input_flux"]

    def setup(self, force_input_flux):
        self.kwargs = utils.render_kwargs(
            ngrid=128,
            psf_obj=utils.gaussian_psf(),
            force_input_flux=force_input_flux,
            backend="np",
        )

    def time_render_force_input_flux(self, force_input_flux):
        batsim.simulate_galaxy(**self.kwargs)


class TimeCenteringOption:
    """Benchmark true-center versus integer-center grid alignment."""

    params = [True, False]
    param_names = ["use_true_center"]

    def setup(self, use_true_center):
        self.kwargs = utils.render_kwargs(
            ngrid=128,
            psf_obj=utils.gaussian_psf(),
            use_true_center=use_true_center,
            backend="np",
        )

    def time_render_use_true_center(self, use_true_center):
        batsim.simulate_galaxy(**self.kwargs)


class TimeIntegrationCompensation:
    """Benchmark integration-compensation modes on a compact profile."""

    params = [None, "exact_sinc", "quadrature"]
    param_names = ["compensate_integration"]

    def setup(self, compensate_integration):
        self.kwargs = utils.render_kwargs(
            gal_obj=utils.compact_sersic_galaxy(),
            ngrid=128,
            psf_obj=utils.gaussian_psf(),
            integration_order=2,
            compensate_integration=compensate_integration,
            backend="np",
        )

    def time_render_compensate_integration(self, compensate_integration):
        batsim.simulate_galaxy(**self.kwargs)


class TimePrecisionOption:
    """Benchmark single versus double precision in a medium PSF render."""

    params = ["single", "double"]
    param_names = ["precision"]

    def setup(self, precision):
        self.kwargs = utils.render_kwargs(
            ngrid=128,
            psf_obj=utils.moffat_psf(),
            precision=precision,
            backend="np",
        )

    def time_render_precision(self, precision):
        batsim.simulate_galaxy(**self.kwargs)


class TimeCuPyRender:
    """Optional GPU render benchmark that skips when CuPy is unavailable."""

    params = [128]
    param_names = ["ngrid"]

    def setup(self, ngrid):
        utils.get_cupy()
        self.kwargs = utils.render_kwargs(
            ngrid=ngrid,
            psf_obj=utils.gaussian_psf(),
            backend="cp",
        )

    def teardown(self, ngrid):
        utils.cleanup_backend("cp")

    def time_render_cupy(self, ngrid):
        batsim.simulate_galaxy(**self.kwargs)

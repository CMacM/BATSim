"""Peak-memory benchmarks for representative CPU renders."""

import batsim

try:
    from benchmarks import utils
except ImportError:
    import utils


class PeakMemRenderSize:
    """Track memory growth for medium and large renders."""

    params = [128, 256]
    param_names = ["ngrid"]
    timeout = 180

    def setup(self, ngrid):
        self.kwargs = utils.render_kwargs(
            gal_obj=utils.sersic_galaxy(),
            ngrid=ngrid,
            psf_obj=utils.gaussian_psf(),
            backend="np",
        )

    def teardown(self, ngrid):
        utils.cleanup_backend("np")

    def peakmem_render_size(self, ngrid):
        batsim.simulate_galaxy(**self.kwargs)


class PeakMemPsfPath:
    """Compare memory for PSF-convolved and no-PSF render paths."""

    params = [False, True]
    param_names = ["with_psf"]
    timeout = 180

    def setup(self, with_psf):
        self.kwargs = utils.render_kwargs(
            gal_obj=utils.sersic_galaxy(),
            ngrid=128,
            psf_obj=utils.gaussian_psf() if with_psf else None,
            draw_method="auto",
            backend="np",
        )

    def teardown(self, with_psf):
        utils.cleanup_backend("np")

    def peakmem_render_psf_path(self, with_psf):
        batsim.simulate_galaxy(**self.kwargs)


class PeakMemIntegrationOrder:
    """Compare memory for low and high block-integration order."""

    params = [1, 4]
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

    def peakmem_render_integration_order(self, integration_order):
        batsim.simulate_galaxy(**self.kwargs)

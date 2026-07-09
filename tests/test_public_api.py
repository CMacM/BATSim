import batsim


def test_top_level_public_api_exports():
    """Check that cleaned top-level exports preserve the intended public API."""
    assert batsim.simulate_galaxy.__name__ == "simulate_galaxy"
    assert batsim.clear_backend_memory.__name__ == "clear_backend_memory"
    assert batsim.Stamp.__name__ == "Stamp"

    assert batsim.LensTransform is batsim.AffineLensingTransform
    assert batsim.IaTransform is batsim.IATransform
    assert batsim.FlexionTransform.__name__ == "FlexionTransform"

    assert hasattr(batsim, "_gsinterface")

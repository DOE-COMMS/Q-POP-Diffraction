"""
CuPy path check: the GPU result matches NumPy on a small random field.

Skipped when CuPy or a GPU is not available.

    pytest tests/test_gpu.py
"""

import numpy as np
import pytest

cp = pytest.importorskip("cupy")
if cp.cuda.runtime.getDeviceCount() == 0:
    pytest.skip("no GPU", allow_module_level=True)


def test_fftn_matches_numpy():
    rng = np.random.default_rng(0)
    a = rng.standard_normal((24, 16, 10)) + 1j * rng.standard_normal((24, 16, 10))
    np.testing.assert_allclose(cp.fft.fftn(cp.asarray(a)).get(), np.fft.fftn(a), rtol=1e-12, atol=1e-10)

# --- Compatibility patch for deprecated numpy aliases ---
import numpy as np

if not hasattr(np, 'float'):
    np.float = float
if not hasattr(np, 'int'):
    np.int = int
if not hasattr(np, 'bool'):
    np.bool = bool
if not hasattr(np, 'object'):
    np.object = object
# --------------------------------------------------------

from . import chunks, components, fitters, lib, models, template, tellurics, bad_pixels, plot_lib, timeseries

# __all__ = ['chunks', 'fitters', 'models', 'components']  # FIXME: Not working?

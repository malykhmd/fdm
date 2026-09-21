##########################
# Polynomialization and quadratization for FDM
##########################

_fdm_pd_poly_loaded = False

for _fdm_pd_poly_file in [
    'sage/polynomialization.sage',
    'polynomialization.sage',
]:
    if not _fdm_pd_poly_loaded:
        try:
            load(_fdm_pd_poly_file)
            _fdm_pd_poly_loaded = True
        except Exception:
            pass

if not _fdm_pd_poly_loaded:
    raise IOError("Cannot find polynomialization.sage. Run Sage from the FDM root folder or from the sage folder.")

_fdm_pd_quad_loaded = False

for _fdm_pd_quad_file in [
    'sage/quadratization.sage',
    'quadratization.sage',
]:
    if not _fdm_pd_quad_loaded:
        try:
            load(_fdm_pd_quad_file)
            _fdm_pd_quad_loaded = True
        except Exception:
            pass

if not _fdm_pd_quad_loaded:
    raise IOError("Cannot find quadratization.sage. Run Sage from the FDM root folder or from the sage folder.")

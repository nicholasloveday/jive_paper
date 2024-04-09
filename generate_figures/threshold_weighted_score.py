"""
Implementation of some basic threshold weighted scoring functions.
See Taggart (2021) https://doi.org/10.1002/qj.4206

In the future, this will be implemented in https://github.com/nci/scores
"""

import functools
from typing import Callable, Optional, Sequence, Tuple, Union

import numpy as np
import xarray as xr

SCORING_FUNCS = ["squared_error"]
AuxFuncType = Callable[[xr.DataArray], xr.DataArray]
EndpointType = Union[int, float, xr.DataArray]


def threshold_weighted_squared_error(
    fcst: xr.DataArray,
    obs: xr.DataArray,
    interval_where_one: Tuple[EndpointType, EndpointType],
    dims: Optional[Sequence[str]] = None,
) -> xr.DataArray:
    """
    Returns the scores computed using a threshold weighted scoring
    function derived from the squared error function. Scores are averaged
    over the `dims` dimension.

    This is a convenience function for `threshold_weighted_score`
    with the option `scoring_func="squared_error"`. The weight is
    1 on the specified interval and 0 elsewhere. For more flexible weighting schemes,
    see `threshold_weighted_score` and `consistent_expectile_score`.

    Args:
        fcst: array of forecast values.
        obs: array of corresponding observation values.
        interval_where_one: endpoints of the interval where the weights are 1.
            Must be increasing. Infinite endpoints are permissible. By supplying a tuple of
            arrays, endpoints can vary with dimension.
        dims: dimensions to preserve in the output.
            All other dimensions are collapsed by taking the mean.

    Returns:
        xarray data array of threshold weighted scores, with the dimensions specified by `dims`.
        If `dims` is `None`, the returned DataArray will have only one entry, the overall mean score.

    Raises:
        ValueError: if `interval_where_one` is not length 2.
        ValueError: if `interval_where_one` is not increasing.

    References:
        Taggart, R. (2021). "Evaluation of point forecasts for extreme events
        using consistent scoring functions", Q. J. Royal Meteorol. Soc.,
        https://doi.org/10.1002/qj.4206
    """
    _check_tws_args("squared_error", interval_where_one, None)

    g, phi, phi_prime = _auxiliary_funcs(fcst, obs, interval_where_one, None)

    result = consistent_expectile_score(fcst, obs, 0.5, phi, phi_prime, dims=dims)

    return result


def consistent_expectile_score(
    fcst: xr.DataArray,
    obs: xr.DataArray,
    alpha: float,
    phi: Callable[[xr.DataArray], xr.DataArray],
    phi_prime: Callable[[xr.DataArray], xr.DataArray],
    dims: Optional[Sequence[str]] = None,
) -> xr.DataArray:
    """
    Returns the score using a scoring function that is consistent for the
    alpha-expectile functional, based on a supplied convex function phi.
    See Geniting (2011), or Equation (10) from Taggart (2021).

    Args:
        fcst: array of forecast values.
        obs: array of corresponding observation values.
        alpha: expectile level. Must be strictly between 0 and 1.
        phi: a convex function on the real numbers, accepting a single array like argument.
        phi_prime: a subderivative of `phi`, accepting a single array like argument.
        dims: dimensions to preserve in the output.
            All other dimensions are collapsed by taking the mean.

    Returns:
        array of (mean) scores that is consistent for alpha-expectile functional,
        with the dimensions specified by `dims`. If `dims` is `None`, the returned DataArray will have
        only one entry, the overall mean score.

    Raises:
        InvalidInputError: if `alpha` is not strictly between 0 and 1.

    Note:
        .. math::

            S(x, y) =
            \\begin{cases}
            (1 - \\alpha)(\\phi(y) - \\phi(x) - \\phi'(x)(y-x)), & y < x \\\\
            \\alpha(\\phi(y) - \\phi(x) - \\phi'(x)(y-x)), & x \\leq y
            \\end{cases}

        where

            - :math:`x` is the forecast
            - :math:`y` is the observation
            - :math:`\\alpha` is the expectile level
            - :math:`\\phi` is a convex function of a single variable
            - :math:`\\phi'` is the subderivative of :math:`\\phi`
            - :math:`S(x,y)` is the score.

        Note that if :math:`\\phi` is differentiable then `\\phi'` is its derivative.

    References:
        -   Gneiting, T. (2011). "Making and evaluating point forecasts",
            J. Amer. Statist. Assoc.,
            https://doi.org/10.1198/jasa.2011.r10138
        -   Taggart, R. (2021). "Evaluation of point forecasts for extreme events
            using consistent scoring functions", Q. J. Royal Meteorol. Soc.,
            https://doi.org/10.1002/qj.4206

    """
    check_alpha(alpha)

    score_overfcst = (1 - alpha) * (
        phi(obs) - phi(fcst) - phi_prime(fcst) * (obs - fcst)
    )
    score_underfcst = alpha * (phi(obs) - phi(fcst) - phi_prime(fcst) * (obs - fcst))
    result = score_overfcst.where(obs < fcst, score_underfcst)

    if dims is None or (set(dims) < set(fcst.dims)):
        dims_to_collapse = dims_complement(result, dims=dims)
        result = result.mean(dim=dims_to_collapse)

    return result


def _check_tws_args(scoring_func, interval_where_one, interval_where_positive):
    """
    Some argument checks for `threshold_weighted_score`.
    Checks for valid interval endpoints are done in `_auxiliary_funcs`.
    """
    if len(interval_where_one) != 2:
        raise ValueError("`interval_where_one` must have length 2")

    if interval_where_positive is not None and len(interval_where_positive) != 2:
        raise ValueError("`interval_where_positive` must be length 2 when not `None`")

    if scoring_func not in SCORING_FUNCS:
        raise ValueError("`scoring_func` must be one of:" + ", ".join(SCORING_FUNCS))


def _auxiliary_funcs(
    fcst: xr.DataArray,
    obs: xr.DataArray,
    interval_where_one: Tuple[EndpointType, EndpointType],
    interval_where_positive: Optional[Tuple[EndpointType, EndpointType]],
) -> Tuple[AuxFuncType, AuxFuncType, AuxFuncType]:
    """
    Returns the three auxiliary functions g, phi and phi_prime
    which are used to construct quantile, expectile or Huber mean scoring functions.
    See Equations (8), (10) and (11) from Taggart (2022) for the role of g, phi and phi_prime.
    """

    if interval_where_positive is None:  # rectangular weight
        a, b = interval_where_one

        if isinstance(a, (float, int)):
            a = xr.DataArray(a)
            b = xr.DataArray(b)

        if (a >= b).any():
            raise ValueError(
                "left endpoint of `interval_where_one` must be strictly less than right endpoint"
            )

        # safest to work with finite a and b
        a = a.where(a > -np.inf, float(min(fcst.min(), obs.min(), b.min())) - 1)
        b = b.where(b < np.inf, float(max(fcst.max(), obs.max(), a.max())) + 1)

        g = functools.partial(_g_j_rect, a, b)
        phi = functools.partial(_phi_j_rect, a, b)
        phi_prime = functools.partial(_phi_j_prime_rect, a, b)

    else:  # trapezoidal weight
        a, d = interval_where_positive
        b, c = interval_where_one

        if isinstance(a, (float, int)):
            a = xr.DataArray(a)
            d = xr.DataArray(d)

        if isinstance(b, (float, int)):
            b = xr.DataArray(b)
            c = xr.DataArray(c)

        if (b >= c).any():
            raise ValueError(
                "left endpoint of `interval_where_one` must be strictly less than right endpoint"
            )

        if (np.isinf(a) & (a != b)).any() or (np.isinf(d) & (c != d)).any():
            raise ValueError(
                "`interval_where_positive` endpoint can only be infinite when "
                "corresponding `interval_where_one` endpoint is infinite."
            )

        if not ((a < b) | ((a == b) & np.isinf(a))).all():
            raise ValueError(
                "left endpoint of `interval_where_positive` must be less than "
                "left endpoint of `interval_where_one`, unless both are `-numpy.inf`."
            )

        if not ((c < d) | ((c == d) & np.isinf(c))).all():
            raise ValueError(
                "right endpoint of `interval_where_positive` must be greater than "
                "right endpoint of `interval_where_one`, unless both are `numpy.inf`."
            )

        # safest to work with finite intervals
        b = b.where(b > -np.inf, min(fcst.min(), obs.min(), c.min()) - 1)
        a = a.where(a > -np.inf, b.min() - 1)
        c = c.where(c < np.inf, max(fcst.max(), obs.max(), b.max()) + 1)
        d = d.where(d < np.inf, c.max() + 1)

        g = functools.partial(_g_j_trap, a, b, c, d)
        phi = functools.partial(_phi_j_trap, a, b, c, d)
        phi_prime = functools.partial(_phi_j_prime_trap, a, b, c, d)

    return g, phi, phi_prime


def _g_j_rect(a: EndpointType, b: EndpointType, x: xr.DataArray) -> xr.DataArray:
    """
    Returns values of a nondecreasing function g_j, where g_j is obtained by integrating
    a rectangular weight function. The weight is 1 on the interval [a, b) and 0 otherwise.
    The formula is based on the first row of Table B1 from Taggart (2022).

    Args:
        a: left endpoint of interval where weight = 1. Can be `-np.inf`.
        b: right endpoint of the interval where weight = 1. Can be `np.inf`.
        x: points where g_j is to be evaluated.

    Returns:
        array of function values for the function g_j.

    Note:
        Requires a < b. This is tested in `_auxiliary_funcs`.
    """
    # results correspond to each case in the first row of Table B1, Taggart (2022).
    result1 = 0
    result2 = x - a
    result3 = b - a

    result = result2.where(x < b, result3).where(x >= a, result1)
    result = result.where(~np.isnan(x), np.nan)

    return result


def _phi_j_rect(a: EndpointType, b: EndpointType, x: xr.DataArray) -> xr.DataArray:
    """
    Returns values of a convex function phi_j, where phi_j is obtained by integrating
    a rectangular weight function. The weight is 1 on the interval [a, b) and 0 otherwise.
    The formula is based on the second row of Table B1 from Taggart (2022).

    Args:
        a: left endpoint of interval where weight = 1. Can be `-np.inf`.
        b: right endpoint of the interval where weight = 1. Can be `np.inf`.
        x: points where phi_j is to be evaluated.

    Returns:
        array of function values for the function phi_j.

    Note:
        Requires a < b. This is tested in `_auxiliary_funcs`.
    """
    # results correspond to each case in the second row of Table B1, Taggart (2022).
    result1 = 0
    result2 = 2 * (x - a) ** 2
    result3 = 4 * (b - a) * x + 2 * (a**2 - b**2)

    result = result2.where(x < b, result3).where(x >= a, result1)
    result = result.where(~np.isnan(x), np.nan)

    return result


def _phi_j_prime_rect(
    a: EndpointType, b: EndpointType, x: xr.DataArray
) -> xr.DataArray:
    """
    The subderivative of `_phi_j_rect(a, b, x)` w.r.t. x.
    """
    return 4 * _g_j_rect(a, b, x)


def _g_j_trap(
    a: EndpointType, b: EndpointType, c: EndpointType, d: EndpointType, x: xr.DataArray
) -> xr.DataArray:
    """
    Returns values of a nondecreasing function g_j, where g_j is obtained by integrating
    a trapezoidal weight function. The weight is 1 on the interval (b, c) and 0 outside
    the interval (a,d). The formula is based on the third row of Table B1 from Taggart (2022).

    Args:
        a: left endpoint of interval where weight > 0.
        b: left endpoint of the interval where weight = 1.
        c: right endpoint of the interval where weight = 1.
        d: right endpoint of interval where weight > 0.
        x: points where g_j is to be evaluated.

    Returns:
        array of function values for the function g_j.

    Note:
        Requires a < b < c < d. This is tested in `_auxiliary_funcs`.
    """
    # results correspond to each case in the third row of Table B1, Taggart (2022).
    result0 = 0
    result1 = (x - a) ** 2 / (2 * (b - a))
    result2 = x - (b + a) / 2
    result3 = -((d - x) ** 2) / (2 * (d - c)) + (d + c - a - b) / 2
    result4 = (d + c - a - b) / 2
    result = (
        result1.where(x >= a, result0)
        .where(x < b, result2)
        .where(x < c, result3)
        .where(x < d, result4)
        .where(~np.isnan(x), np.nan)
    )
    return result


def _phi_j_trap(
    a: EndpointType, b: EndpointType, c: EndpointType, d: EndpointType, x: xr.DataArray
) -> xr.DataArray:
    """
    Returns values of a convex function phi_j, where phi_j is obtained by integrating
    a trapezoidal weight function. The weight is 1 on the interval (b, c) and 0 outside
    the interval (a,d). The formula is based on the fourth row of Table B1 from Taggart (2022).

    Args:
        a: left endpoint of interval where weight > 0.
        b: left endpoint of the interval where weight = 1.
        c: right endpoint of the interval where weight = 1.
        d: right endpoint of interval where weight > 0.
        x: points where phi_j is to be evaluated.

    Returns:
        array of function values for the function phi_j.

    Note:
        Requires a < b < c < d. This is tested in `_auxiliary_funcs`.
    """
    # results correspond to each case in the fourth row of Table B1, Taggart (2022).
    result0 = 0
    result1 = 2 * (x - a) ** 3 / (3 * (b - a))
    result2 = 2 * x**2 - 2 * (a + b) * x + 2 * (b - a) ** 2 / 3 + 2 * a * b
    result3 = (
        2 * (d - x) ** 3 / (3 * (d - c))
        + 2 * (d + c - a - b) * x
        + 2 * ((b - a) ** 2 + 3 * a * b - (d - c) ** 2 - 3 * c * d) / 3
    )
    result4 = (
        2 * (d + c - a - b) * x
        + 2 * ((b - a) ** 2 + 3 * a * b - (d - c) ** 2 - 3 * c * d) / 3
    )
    result = (
        result1.where(x >= a, result0)
        .where(x < b, result2)
        .where(x < c, result3)
        .where(x < d, result4)
        .where(~np.isnan(x), np.nan)
    )
    return result


def _phi_j_prime_trap(
    a: EndpointType, b: EndpointType, c: EndpointType, d: EndpointType, x: xr.DataArray
) -> xr.DataArray:
    """
    The subderivative of `_phi_j_trap(a, b, c, d, x)` w.r.t. x.
    """
    return 4 * _g_j_trap(a, b, c, d, x)


def check_alpha(alpha: float) -> None:
    """Raises if quantile or expectile level `alpha` not in the open interval (0,1)."""
    if alpha <= 0 or alpha >= 1:
        raise ValueError("`alpha` must be strictly between 0 and 1")


def dims_complement(data, dims=None):
    """
    Returns the complement of data.dims and dims

    Args:
        data: an xarray DataArray or Dataset
        dims: an Iterable of strings corresponding to dimension names

    Returns:
        A sorted list of dimension names, the complement of data.dims and dims
    """
    if dims is None:
        dims = []

    complement = set(data.dims) - set(dims)
    return sorted(list(complement))

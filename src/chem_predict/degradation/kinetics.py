"""First-order thermal exposure with user-supplied kinetic parameters."""

from collections.abc import Sequence
from dataclasses import dataclass
from math import exp, expm1, isfinite


@dataclass(frozen=True, slots=True)
class ThermalStage:
    temperature_c: float
    duration_h: float


@dataclass(frozen=True, slots=True)
class FirstOrderExposure:
    integrated_rate: float
    remaining_fraction: float
    converted_fraction: float
    stage_rate_constants_per_h: tuple[float, ...]


def first_order_exposure(stages: Sequence[ThermalStage], *, k_reference_per_h: float,
                         reference_temperature_c: float,
                         activation_energy_kj_mol: float) -> FirstOrderExposure:
    """Integrate piecewise constant-temperature Arrhenius first-order loss.

    Requires measured/fitted parameters for the actual process. Conversion is
    parent disappearance, not yield of any particular degradation product.
    """
    for name, value in (("k_reference_per_h", k_reference_per_h),
                        ("activation_energy_kj_mol", activation_energy_kj_mol)):
        if not isfinite(value) or value < 0:
            raise ValueError(f"{name} must be finite and non-negative")
    if not isfinite(reference_temperature_c) or reference_temperature_c <= -273.15:
        raise ValueError("Reference temperature must be above absolute zero")
    rates, integrated = [], 0.0
    for stage in stages:
        if not isfinite(stage.temperature_c) or stage.temperature_c <= -273.15:
            raise ValueError("Stage temperature must be above absolute zero")
        if not isfinite(stage.duration_h) or stage.duration_h < 0:
            raise ValueError("Stage duration_h must be finite and non-negative")
        exponent = activation_energy_kj_mol * 1000 / 8.314462618 * (
            1 / (reference_temperature_c + 273.15) - 1 / (stage.temperature_c + 273.15))
        try:
            rate = k_reference_per_h * exp(exponent) if k_reference_per_h else 0.0
        except OverflowError as exc:
            raise ValueError("Arrhenius extrapolation exceeds numerical range") from exc
        if not isfinite(rate):
            raise ValueError("Arrhenius rate exceeds numerical range")
        rates.append(rate)
        integrated += rate * stage.duration_h
    return FirstOrderExposure(integrated, exp(-integrated), -expm1(-integrated), tuple(rates))

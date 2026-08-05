"""
Collisionless plume-model dispatch and the shared local field state.

PyRPOD ships two collisionless analytical plume models, both from Cai & Wang
2012 and both already verified in this repository:

* :class:`~pyrpod.plume.RarefiedPlumeGasKinetics.SimplifiedGasKinetics` --
  the far-field simplification (Eq. 13's ``Q'``), and
* :class:`~pyrpod.plume.RarefiedPlumeGasKinetics.CollisionlessGasKinetics` --
  the full model, which integrates the exact special factor Q over the finite
  exit disk and therefore stays valid in the near field.

The second is a subclass of the first with an identical constructor, so
selecting between them is a class lookup, not a plugin framework. This module
is that lookup (:data:`PLUME_MODELS`, :func:`create_model`) plus the two
pieces of model-INDEPENDENT physics the strike pipeline needs:

``local_field_state(model)``
    Reduces whichever model was selected to one common
    :class:`LocalFieldState` -- number density, mass density, axial/radial
    velocity, velocity magnitude, temperature and local speed ratio at the
    evaluated point. Every model-specific getter is called through the
    instance, so a model that overrides the field solutions is honoured
    without this function knowing which one it holds.

``maxwellian_surface_loads(state, ...)``
    Applies the Maxwellian (Shen) gas-surface interaction to that common
    state. It is written once here, so pressure, shear and heat-transfer
    logic is never duplicated per model.

Scope
-----
Both models are COLLISIONLESS. Nothing here consumes a Knudsen number,
applies a collisional correction, or reads DSMC data; the study layer records
Kn purely as derived metadata (see :mod:`pyrpod.mdao.study_config`).

Naming
------
Two vocabularies meet here and are deliberately kept distinct:

* the *model name* is the Python class name, which is what a study
  configuration writes (``plume_model.name: CollisionlessGasKinetics``);
* the *kinetics key* is the short token a case's ``config.ini`` carries in
  ``[pm] kinetics`` (``Simplified``, ``Collisionless``, or ``None`` to
  disable surface loads entirely).

:func:`model_name_for_kinetics` and :func:`kinetics_key_for` convert between
them. ``Simplified`` keeps meaning exactly what it always meant, so every
existing case and its results are unchanged.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np

from pyrpod.plume.RarefiedPlumeGasKinetics import (
    AVOGADROS_NUMBER,
    CollisionlessGasKinetics,
    Scalar,
    SimplifiedGasKinetics,
    get_maxwellian_heat_transfer,
    get_maxwellian_pressure,
    get_maxwellian_shear_pressure,
)

__all__ = [
    "DEFAULT_PLUME_MODEL",
    "KINETICS_DISABLED",
    "KINETICS_KEYS",
    "LocalFieldState",
    "PLUME_MODELS",
    "PlumeModelError",
    "create_model",
    "kinetics_key_for",
    "local_field_state",
    "maxwellian_surface_loads",
    "model_name_for_kinetics",
    "resolve_model_name",
]


class PlumeModelError(ValueError):
    """Raised when an unknown plume model or kinetics key is requested."""


#: The collisionless plume models the strike pipeline can dispatch to, keyed
#: by class name (what a study configuration names).
PLUME_MODELS: dict[str, type[SimplifiedGasKinetics]] = {
    "SimplifiedGasKinetics": SimplifiedGasKinetics,
    "CollisionlessGasKinetics": CollisionlessGasKinetics,
}

#: Model used when nothing selects one -- the historical behavior of every
#: existing case and of the whole strike pipeline before model dispatch.
DEFAULT_PLUME_MODEL = "SimplifiedGasKinetics"

#: ``[pm] kinetics`` value that disables the gas-kinetics surface loads.
KINETICS_DISABLED = "None"

#: ``config.ini`` ``[pm] kinetics`` token -> model class name.
KINETICS_KEYS: dict[str, str] = {
    "Simplified": "SimplifiedGasKinetics",
    "Collisionless": "CollisionlessGasKinetics",
}


def resolve_model_name(name: str | None) -> str:
    """Validate a plume-model class name, defaulting when none is given.

    Raises
    ------
    PlumeModelError
        If ``name`` is not one of :data:`PLUME_MODELS`. Unknown models are
        never silently replaced by the default.
    """
    if name is None:
        return DEFAULT_PLUME_MODEL
    model_name = str(name)
    if model_name not in PLUME_MODELS:
        raise PlumeModelError(
            f"unknown plume model {model_name!r}; supported models are "
            f"{sorted(PLUME_MODELS)}")
    return model_name


def kinetics_key_for(model_name: str) -> str:
    """``[pm] kinetics`` token that selects ``model_name``."""
    model_name = resolve_model_name(model_name)
    for key, name in KINETICS_KEYS.items():
        if name == model_name:
            return key
    raise PlumeModelError(  # pragma: no cover - PLUME_MODELS/KINETICS_KEYS agree
        f"no [pm] kinetics key is mapped to plume model {model_name!r}")


def model_name_for_kinetics(kinetics: str) -> str:
    """Plume-model class name selected by a ``[pm] kinetics`` token.

    ``'None'`` has no model (surface loads are disabled) and is rejected
    here; callers test for it before asking.
    """
    key = str(kinetics)
    if key == KINETICS_DISABLED:
        raise PlumeModelError(
            "[pm] kinetics = None disables the gas-kinetics surface loads; "
            "there is no plume model to select")
    if key not in KINETICS_KEYS:
        raise PlumeModelError(
            f"unknown [pm] kinetics value {key!r}; supported values are "
            f"{sorted(KINETICS_KEYS)} (or {KINETICS_DISABLED!r} to disable "
            "gas-kinetics surface loads)")
    return KINETICS_KEYS[key]


def create_model(model_name: str, distance: Scalar, theta: Scalar,
                 thruster_characteristics: Mapping[str, Any], T_w: float,
                 sigma: float) -> SimplifiedGasKinetics:
    """Instantiate the named plume model at one evaluation point.

    Every supported model shares the constructor signature
    ``(distance, theta, thruster_characteristics, T_w, sigma)``, so the
    factory is a class lookup and a call.

    Parameters
    ----------
    model_name : str
        A key of :data:`PLUME_MODELS`.
    distance : float
        Distance from the nozzle exit center to the evaluated point (m).
    theta : float
        Plume-centerline off-axis angle of the evaluated point (rad).
    thruster_characteristics : mapping
        The case's thruster-definition entry (``d``, ``ve``, ``R``,
        ``gamma``, ``Te``, ``n``).
    T_w, sigma : float
        Surface temperature (K) and diffuse-reflection fraction.
    """
    return PLUME_MODELS[resolve_model_name(model_name)](
        distance, theta, thruster_characteristics, T_w, sigma)


@dataclass(frozen=True)
class LocalFieldState:
    """Plume flow state at one point, common to every collisionless model.

    This is the interface between a plume model and the gas-surface
    interaction: once a model has been reduced to these numbers, the surface
    loads no longer depend on which model produced them.

    Attributes
    ----------
    number_density : float
        Local number density n (particles / m^3).
    mass_density : float
        Local mass density rho = n * M / N_A (kg / m^3).
    axial_velocity : float
        Macroscopic velocity component along the plume axis, U (m/s).
    radial_velocity : float
        Macroscopic velocity component transverse to the plume axis, W
        (m/s). Zero on the centerline.
    velocity_magnitude : float
        |(U, W)| (m/s); the speed the Maxwellian wall formulas use.
    temperature : float
        Local translational temperature T (K).
    speed_ratio : float
        Local molecular speed ratio S = |U| * beta(T), beta = 1/sqrt(2 R T).
    on_centerline : bool
        Whether the point was evaluated with the exact centerline closed
        forms (theta == 0) rather than the off-axis field solutions.
    """

    number_density: float
    mass_density: float
    axial_velocity: float
    radial_velocity: float
    velocity_magnitude: float
    temperature: float
    speed_ratio: float
    on_centerline: bool = False

    @property
    def velocity(self) -> tuple[float, float]:
        """(axial, radial) velocity components in the plume frame (m/s)."""
        return (self.axial_velocity, self.radial_velocity)

    def to_dict(self) -> dict[str, float | bool]:
        """Plain-data form, for recording a sampled field state."""
        return {
            "number_density": self.number_density,
            "mass_density": self.mass_density,
            "axial_velocity": self.axial_velocity,
            "radial_velocity": self.radial_velocity,
            "velocity_magnitude": self.velocity_magnitude,
            "temperature": self.temperature,
            "speed_ratio": self.speed_ratio,
            "on_centerline": self.on_centerline,
        }


def local_field_state(model: SimplifiedGasKinetics) -> LocalFieldState:
    """Reduce any collisionless plume model to the common local field state.

    Off the centerline the model's own field solutions are used; ON the
    centerline (``theta == 0``) the exact closed forms are, exactly as the
    models' own ``get_pressure`` / ``get_shear_pressure`` / ``get_heat_flux``
    do. Both branches call through the instance, so an overriding model (the
    full :class:`CollisionlessGasKinetics`) supplies its own field values
    while the reduction itself stays model-independent.

    Note that the two field getters are normalized differently by the models
    (``get_U_normalized`` returns U * sqrt(beta_0), the centerline form
    returns U * beta_0); the division by ``beta_0`` reproduces the
    established pipeline behavior in both branches unchanged.
    """
    if model.theta != 0:                       # off the plume centerline
        number_density = model.n_0 * model.get_num_density_ratio()
        temperature = model.T_0 * model.get_temp_ratio()
        axial = model.get_U_normalized() / model.beta_0
        radial = model.get_W_normalized() / model.beta_0
        speed = float(np.sqrt(axial ** 2 + radial ** 2))
        on_centerline = False
    else:                                      # exact centerline closed forms
        number_density = model.n_0 * model.get_num_density_centerline()
        temperature = model.T_0 * model.get_temp_centerline()
        speed = model.get_velocity_centerline() / model.beta_0
        axial, radial = speed, 0.0
        on_centerline = True

    mass_density = number_density * model.molar_mass / AVOGADROS_NUMBER
    return LocalFieldState(
        number_density=float(number_density),
        mass_density=float(mass_density),
        axial_velocity=float(axial),
        radial_velocity=float(radial),
        velocity_magnitude=float(speed),
        temperature=float(temperature),
        speed_ratio=float(speed * model.get_beta(temperature)),
        on_centerline=on_centerline,
    )


def maxwellian_surface_loads(state: LocalFieldState, *, sigma: float,
                             T_w: float, R: float, gamma: float,
                             incidence: Scalar) -> tuple[float, float, float]:
    """Maxwellian (Shen) surface loads for one local field state.

    The single place pressure, shear and heat transfer are computed from a
    plume field, so no model adapter reimplements them.

    Parameters
    ----------
    state : LocalFieldState
        Plume flow state at the face, from :func:`local_field_state`.
    sigma : float
        Fraction of diffuse molecular reflections, in [0, 1].
    T_w : float
        Surface (wall) temperature (K).
    R : float
        Specific gas constant (J / kg / K).
    gamma : float
        Ratio of specific heats.
    incidence : float
        Angle between the LOCAL FLOW DIRECTION (radial from the nozzle exit,
        because the flow is collisionless) and the face unit normal (rad).
        This is the true incidence angle, not the positional off-axis angle
        that locates the face in the plume field.

    Returns
    -------
    (float, float, float)
        Pressure (Pa), shear stress (Pa, signed as the Shen formula returns
        it -- callers take the magnitude) and heat flux (W/m^2).
    """
    pressure = get_maxwellian_pressure(state.mass_density,
                                       state.velocity_magnitude,
                                       state.speed_ratio, sigma, incidence,
                                       state.temperature, T_w)
    shear = get_maxwellian_shear_pressure(state.mass_density,
                                          state.velocity_magnitude,
                                          state.speed_ratio, sigma, incidence)
    heat_flux = get_maxwellian_heat_transfer(state.mass_density,
                                             state.speed_ratio, sigma,
                                             incidence, state.temperature,
                                             T_w, R, gamma)
    return pressure, shear, heat_flux

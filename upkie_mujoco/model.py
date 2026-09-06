"""Build a MuJoCo model from the official Upkie robot description."""

from __future__ import annotations

import mujoco
import upkie_description

TIMESTEP = 0.002
INITIAL_BASE_HEIGHT = 0.6
LEG_KP = 80.0
LEG_KD = 2.0
LEG_TORQUE_LIMIT = 16.0
WHEEL_TORQUE_LIMIT = 1.7

LEG_JOINTS = (
    "left_hip",
    "left_knee",
    "right_hip",
    "right_knee",
)
WHEEL_JOINTS = ("left_wheel", "right_wheel")


def _gain_parameters(value: float) -> list[float]:
    """Return a MuJoCo gain/bias parameter vector."""
    return [value, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]


def build_model() -> mujoco.MjModel:
    """Compile an Upkie floating-base model with a floor and six actuators.

    The official URDF is loaded through ``upkie_description``. Fixed URDF links
    are fused by MuJoCo at compile time, preserving their combined mass and
    inertia while leaving the six actuated joints in the model.
    """
    spec = mujoco.MjSpec.from_file(upkie_description.URDF_PATH)
    spec.modelname = "upkie_mujoco"
    spec.option.timestep = TIMESTEP
    spec.option.integrator = mujoco.mjtIntegrator.mjINT_RK4
    spec.option.gravity = [0.0, 0.0, -9.81]

    base = spec.body("base")
    base.pos = [0.0, 0.0, INITIAL_BASE_HEIGHT]
    # The URDF base is a massless virtual link. MuJoCo requires a small,
    # positive inertia once a free joint is attached to it.
    base.mass = 0.001
    base.inertia = [1e-6, 1e-6, 1e-6]
    base.add_freejoint(name="base_free_joint")

    spec.worldbody.add_geom(
        name="floor",
        type=mujoco.mjtGeom.mjGEOM_PLANE,
        size=[0.0, 0.0, 0.05],
        friction=[1.0, 0.005, 0.0001],
        solref=[0.004, 1.0],
        rgba=[0.18, 0.18, 0.18, 1.0],
    )

    for joint_name in LEG_JOINTS:
        joint = spec.joint(joint_name)
        joint.armature = 0.01
        joint.damping = [0.05, 0.0, 0.0]
        position_limit = 1.26 if joint_name.endswith("hip") else 2.51
        spec.add_actuator(
            name=f"{joint_name}_position",
            trntype=mujoco.mjtTrn.mjTRN_JOINT,
            target=joint_name,
            gaintype=mujoco.mjtGain.mjGAIN_FIXED,
            gainprm=_gain_parameters(LEG_KP),
            biastype=mujoco.mjtBias.mjBIAS_AFFINE,
            biasprm=[0.0, -LEG_KP, -LEG_KD, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            ctrllimited=True,
            ctrlrange=[-position_limit, position_limit],
            forcelimited=True,
            forcerange=[-LEG_TORQUE_LIMIT, LEG_TORQUE_LIMIT],
        )

    for joint_name in WHEEL_JOINTS:
        joint = spec.joint(joint_name)
        joint.armature = 0.001
        joint.damping = [0.001, 0.0, 0.0]
        spec.add_actuator(
            name=f"{joint_name}_torque",
            trntype=mujoco.mjtTrn.mjTRN_JOINT,
            target=joint_name,
            gear=[1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            ctrllimited=True,
            ctrlrange=[-WHEEL_TORQUE_LIMIT, WHEEL_TORQUE_LIMIT],
        )

    return spec.compile()


def set_neutral_leg_targets(model: mujoco.MjModel, data: mujoco.MjData) -> None:
    """Command the four leg joints to Upkie's zero-angle configuration."""
    for joint_name in LEG_JOINTS:
        data.ctrl[model.actuator(f"{joint_name}_position").id] = 0.0

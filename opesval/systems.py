"""Validation systems, taken from the OPES and OPES Explore papers.

OPES paper:    M. Invernizzi and M. Parrinello, J. Phys. Chem. Lett. 11, 2731 (2020)
               https://github.com/invemichele/opes
Explore paper: M. Invernizzi and M. Parrinello, J. Chem. Theory Comput. 18, 3988 (2022)
               https://github.com/invemichele/OPES-explore
"""

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
SYSTEMS_DIR = ROOT / "systems"
KB = 0.0083144626181532  # kJ/(mol K), same as GROMACS


@dataclass(frozen=True)
class Method:
    name: str  # used in paths and reports
    plumed: str  # PLUMED input file in the system directory
    paper: str  # where the setup comes from
    optional: bool = False  # run only when explicitly requested


@dataclass(frozen=True)
class Tier:
    replicas: int
    nsteps: int


@dataclass(frozen=True)
class TwoBasins:
    """DeltaF = F(B) - F(A) between two regions of the first CV, as in the papers.

    Transitions are counted between the cores of the two basins, with hysteresis.
    """

    cv: str  # COLVAR field of the CV defining the basins
    in_B: callable  # cv -> bool, region B for DeltaF (A is the complement)
    core_A: callable  # cv -> bool
    core_B: callable  # cv -> bool
    grid: tuple  # FES grid for each CV of the bias: (min, max, nbins, periodic)
    ref_deltaF: callable  # () -> float, reference DeltaF in kBT units


@dataclass(frozen=True)
class Ala4Basins:
    """The 8 metastable states of alanine tetrapeptide, from the phi dihedrals.

    Same definition as get_ala4_basin.sh in the Explore paper repository.
    """

    phis: tuple = ("phi1", "phi2", "phi3")
    to_0: tuple = (-2.2, -1.8)  # enter state 0 when phi is in this interval
    to_1: tuple = (0.8, 1.2)  # enter state 1 when phi is in this interval


@dataclass(frozen=True)
class System:
    name: str
    engine: str  # "ves" (plumed ves_md_linearexpansion) or "gromacs"
    description: str
    kbt: float  # in the energy units of the simulation
    time_unit: str
    methods: tuple
    tiers: dict = field(hash=False)
    analysis: object = None
    paper_data: dict = field(default_factory=dict, hash=False)  # method -> {kind: file}

    @property
    def path(self):
        return SYSTEMS_DIR / self.name

    def method(self, name):
        for m in self.methods:
            if m.name == name:
                return m
        raise KeyError(f"unknown method '{name}' for system '{self.name}'")


def _ref_from_1d_fes(filename, in_B, comments="#"):
    def ref():
        data = np.loadtxt(SYSTEMS_DIR / filename, comments=comments)
        x, fes = data[:, 0], data[:, 1]
        return np.logaddexp.reduce(-fes[~in_B(x)]) - np.logaddexp.reduce(-fes[in_B(x)])

    return ref


def _mueller_ref():
    data = np.loadtxt(SYSTEMS_DIR / "mueller/ref_fes.grid", comments="#!")
    x, fes = data[:, 0], data[:, 1]
    inside = (x >= -1.3) & (x <= 1.0)  # walls
    x, fes = x[inside], fes[inside]
    return np.logaddexp.reduce(-fes[x < 0]) - np.logaddexp.reduce(-fes[x > 0])


def _ala2_ref():
    # reference for the masterclass-22-03 setup (Verlet scheme), from a long unbiased run
    data = np.loadtxt(SYSTEMS_DIR / "ala2/ref_fes2D.dat", comments="#!")
    phi, fes = data[:, 0], data[:, 2] / (KB * 300)
    finite = np.isfinite(fes)
    phi, fes = phi[finite], fes[finite]
    in_B = (phi > 0) & (phi < 2.3)
    return np.logaddexp.reduce(-fes[~in_B]) - np.logaddexp.reduce(-fes[in_B])


ALA2_IN_B = lambda phi: (phi > 0) & (phi < 2.3)  # noqa: E731, as in the OPES paper

SYSTEMS = {
    s.name: s
    for s in [
        System(
            name="wolfe-quapp",
            engine="ves",
            description="modified Wolfe-Quapp 2D potential, biasing x (OPES paper, Fig. 3, S2)",
            kbt=1.0,
            time_unit="",
            methods=(Method("opes", "plumed-opes.dat", "OPES"),),
            tiers={"full": Tier(100, 20_000_000), "quick": Tier(100, 4_000_000), "smoke": Tier(4, 400_000)},
            analysis=TwoBasins(
                cv="p.x",
                in_B=lambda x: x > 0,
                core_A=lambda x: x < -1.2,
                core_B=lambda x: x > 1.2,
                grid=((-3.0, 3.0, 100, False),),
                ref_deltaF=_ref_from_1d_fes("wolfe-quapp/ref_fes.dat", lambda x: x > 0),
            ),
            paper_data={"opes": {"deltaF_fes": "deltaF_fes-opes.dat"}},
        ),
        System(
            name="mueller",
            engine="ves",
            description="Mueller-Brown 2D potential, biasing x (Explore paper, Fig. 2, 3)",
            kbt=1.0,
            time_unit="",
            methods=(
                Method("opes", "plumed-opes.dat", "EXPLORE"),
                Method("explore", "plumed-explore.dat", "EXPLORE"),
            ),
            tiers={"full": Tier(25, 200_000_000), "quick": Tier(25, 20_000_000), "smoke": Tier(4, 1_000_000)},
            analysis=TwoBasins(
                cv="p.x",
                in_B=lambda x: x > 0,
                core_A=lambda x: x < -0.35,
                core_B=lambda x: x > 0.35,
                grid=((-1.3, 1.0, 101, False),),
                ref_deltaF=_mueller_ref,
            ),
            paper_data={
                m: {"deltaF_fes": f"deltaF_fes-{m}.dat", "deltaF_rew": f"deltaF_rew-{m}.dat"}
                for m in ("opes", "explore")
            },
        ),
        System(
            name="ala2",
            engine="gromacs",
            description="alanine dipeptide in vacuum, biasing phi,psi (OPES paper Fig. 1, S4; Explore paper Fig. 4)",
            kbt=KB * 300,
            time_unit="ps",
            methods=(
                Method("opes_jpcl", "plumed-opes_jpcl.dat", "OPES"),
                Method("opes", "plumed-opes.dat", "EXPLORE"),
                Method("explore", "plumed-explore.dat", "EXPLORE"),
            ),
            tiers={"full": Tier(10, 5_000_000), "quick": Tier(10, 5_000_000), "smoke": Tier(4, 250_000)},
            analysis=TwoBasins(
                cv="phi",
                in_B=ALA2_IN_B,
                core_A=lambda phi: (phi < -1.0) | (phi > 2.6),
                core_B=lambda phi: (phi > 0.7) & (phi < 1.6),
                grid=((-np.pi, np.pi, 100, True), (-np.pi, np.pi, 100, True)),
                ref_deltaF=_ala2_ref,
            ),
            paper_data={"opes_jpcl": {"deltaF_fes": "deltaF_fes-opes_jpcl.dat"}},
        ),
        System(
            name="ala4",
            engine="gromacs",
            description="alanine tetrapeptide in vacuum, biasing the psi dihedrals (Explore paper, Fig. 6)",
            kbt=KB * 300,
            time_unit="ps",
            methods=(
                Method("opes", "plumed-opes.dat", "EXPLORE"),
                Method("explore", "plumed-explore.dat", "EXPLORE"),
                Method("opes_jpcl", "plumed-opes_jpcl.dat", "OPES", optional=True),
            ),
            tiers={"full": Tier(10, 10_000_000), "quick": Tier(10, 10_000_000), "smoke": Tier(4, 250_000)},
            analysis=Ala4Basins(),
        ),
    ]
}


def select(systems=None, methods=None, include_optional=False):
    """Yield (system, method) pairs matching the selection."""
    for sname in systems or SYSTEMS:
        if sname not in SYSTEMS:
            raise KeyError(f"unknown system '{sname}', available: {', '.join(SYSTEMS)}")
        system = SYSTEMS[sname]
        for m in system.methods:
            if methods and m.name not in methods:
                continue
            if m.optional and not (include_optional or methods):
                continue
            yield system, m

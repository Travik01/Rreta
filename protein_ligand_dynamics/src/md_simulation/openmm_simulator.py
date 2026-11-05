"""Молекулярная динамика на OpenMM."""

import logging

try:
    import openmm as mm
    from openmm import app, unit
except ImportError:
    mm = app = unit = None

logger = logging.getLogger(__name__)


class OpenMMSimulator:
    def __init__(
        self,
        complex_pdb: str,
        forcefield: str = "amber14-all.xml",
        water_model: str = "amber14/tip3p.xml",
        use_gpu: bool = True,
    ):
        if mm is None:
            raise ImportError("Для OpenMM установите: pip install openmm")

        self.complex_pdb = complex_pdb
        self.forcefield_files = [forcefield, water_model]
        self.use_gpu = use_gpu

        pdb = app.PDBFile(complex_pdb)
        self.topology = pdb.topology
        self.positions = pdb.positions
        self.system = None
        self.simulation = None

    def prepare(self, box_padding: float = 1.0, ionic_strength: float = 0.15):
        logger.info("Подготовка системы OpenMM")
        modeller = app.Modeller(self.topology, self.positions)
        forcefield = app.ForceField(*self.forcefield_files)
        modeller.addSolvent(
            forcefield,
            padding=box_padding * unit.nanometer,
            ionicStrength=ionic_strength * unit.molar,
        )

        self.topology = modeller.topology
        self.positions = modeller.positions
        logger.info("Атомов после сольватации: %d", self.topology.getNumAtoms())

        self.system = forcefield.createSystem(
            self.topology,
            nonbondedMethod=app.PME,
            nonbondedCutoff=1.0 * unit.nanometer,
            constraints=app.HBonds,
            rigidWater=True,
        )
        integrator = mm.LangevinMiddleIntegrator(
            300 * unit.kelvin,
            1.0 / unit.picosecond,
            2.0 * unit.femtosecond,
        )

        if self.use_gpu:
            try:
                platform = mm.Platform.getPlatformByName("CUDA")
                properties = {"CudaPrecision": "mixed"}
                logger.info("Платформа: CUDA")
            except mm.OpenMMException:
                platform = mm.Platform.getPlatformByName("CPU")
                properties = {}
                logger.info("Платформа: CPU")
        else:
            platform = mm.Platform.getPlatformByName("CPU")
            properties = {}
            logger.info("Платформа: CPU")

        self.simulation = app.Simulation(
            self.topology,
            self.system,
            integrator,
            platform,
            properties,
        )
        self.simulation.context.setPositions(self.positions)
        self.simulation.minimizeEnergy(maxIterations=1000)
        logger.info("Система подготовлена")

    def equilibrate(
        self,
        npt_steps: int = 10000,
        nvt_steps: int = 10000,
        temperature: float = 300.0,
        pressure: float = 1.0,
    ):
        logger.info("Уравновешивание: NVT %d, NPT %d шагов", nvt_steps, npt_steps)
        self.simulation.context.setVelocitiesToTemperature(temperature * unit.kelvin)
        self.simulation.step(nvt_steps)

        self.system.addForce(
            mm.MonteCarloBarostat(
                pressure * unit.bar,
                temperature * unit.kelvin,
                25,
            )
        )
        self.simulation.context.reinitialize(preserveState=True)
        self.simulation.step(npt_steps)

    def run_production(
        self,
        num_steps: int = 500000,
        report_interval: int = 1000,
        trajectory_file: str = "trajectory.dcd",
        log_file: str = "simulation.log",
    ) -> dict:
        logger.info("Production MD: %d шагов", num_steps)
        self.simulation.reporters.append(app.DCDReporter(trajectory_file, report_interval))
        self.simulation.reporters.append(
            app.StateDataReporter(
                log_file,
                report_interval,
                step=True,
                time=True,
                potentialEnergy=True,
                kineticEnergy=True,
                temperature=True,
                volume=True,
                speed=True,
            )
        )
        self.simulation.step(num_steps)
        logger.info("Траектория сохранена: %s", trajectory_file)
        return {"trajectory_file": trajectory_file, "num_steps": num_steps}

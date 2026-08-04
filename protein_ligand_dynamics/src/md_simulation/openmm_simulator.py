"""
Классическая MD с OpenMM (полная версия с отсутствовавшими деталями).
"""

import numpy as np
from pathlib import Path
from typing import Optional, Dict

# OpenMM загружается лениво, чтобы модуль можно было импортировать
# без установленного OpenMM (ошибка возникнет только при использовании).
try:
    import openmm as mm
    from openmm import app, unit
    _OPENMM_AVAILABLE = True
except ImportError:
    mm = app = unit = None
    _OPENMM_AVAILABLE = False


class OpenMMSimulator:
    """
    MD-симулятор на базе OpenMM с AMBER силовым полем.
    """
    
    def __init__(
        self,
        complex_pdb: str,
        forcefield: str = "amber14-all.xml",
        water_model: str = "amber14/tip3p.xml",
        use_gpu: bool = True
    ):
        if not _OPENMM_AVAILABLE:
            raise ImportError(
                "Для OpenMM-симуляции установите: pip install openmm"
            )

        self.complex_pdb = complex_pdb
        self.forcefield_files = [forcefield, water_model]
        self.use_gpu = use_gpu
        
        # Загрузка структуры
        self.pdb = app.PDBFile(complex_pdb)
        self.topology = self.pdb.topology
        self.positions = self.pdb.positions
        
        self.system = None
        self.simulation = None
    
    def prepare(
        self,
        box_padding: float = 1.0,
        ionic_strength: float = 0.15,
        pH: float = 7.0
    ):
        """Подготовка системы (см. полную версию выше)."""
        print(f"[OpenMM] Подготовка системы...")
        
        # Добавление растворителя
        modeller = app.Modeller(self.topology, self.positions)
        
        forcefield = app.ForceField(*self.forcefield_files)
        
        modeller.addSolvent(
            forcefield,
            padding=box_padding * unit.nanometer,
            ionicStrength=ionic_strength * unit.molar
        )
        
        self.topology = modeller.topology
        self.positions = modeller.positions
        
        print(f"          Атомов: {self.topology.getNumAtoms():,}")
        
        # Создание системы
        self.system = forcefield.createSystem(
            self.topology,
            nonbondedMethod=app.PME,
            nonbondedCutoff=1.0 * unit.nanometer,
            constraints=app.HBonds,
            rigidWater=True
        )
        
        # Интегратор
        integrator = mm.LangevinMiddleIntegrator(
            300 * unit.kelvin,
            1.0 / unit.picosecond,
            2.0 * unit.femtosecond
        )
        
        # Платформа
        if self.use_gpu:
            try:
                platform = mm.Platform.getPlatformByName('CUDA')
                properties = {'CudaPrecision': 'mixed'}
                print(f"          Платформа: CUDA")
            except:
                platform = mm.Platform.getPlatformByName('CPU')
                properties = {}
                print(f"          Платформа: CPU (медленно!)")
        else:
            platform = mm.Platform.getPlatformByName('CPU')
            properties = {}
        
        # Симуляция
        self.simulation = app.Simulation(
            self.topology,
            self.system,
            integrator,
            platform,
            properties
        )
        
        self.simulation.context.setPositions(self.positions)
        
        # Минимизация
        print(f"          Минимизация энергии...")
        self.simulation.minimizeEnergy(maxIterations=1000)
        
        print(f"[OpenMM] ✓ Система готова")
    
    def equilibrate(
        self,
        npt_steps: int = 10000,
        nvt_steps: int = 10000,
        temperature: float = 300.0,
        pressure: float = 1.0
    ):
        """Уравновешивание (NVT + NPT)."""
        print(f"[OpenMM] Уравновешивание...")
        
        # NVT
        self.simulation.context.setVelocitiesToTemperature(temperature * unit.kelvin)
        self.simulation.step(nvt_steps)
        
        # NPT
        self.system.addForce(mm.MonteCarloBarostat(
            pressure * unit.bar,
            temperature * unit.kelvin,
            25
        ))
        
        self.simulation.context.reinitialize(preserveState=True)
        self.simulation.step(npt_steps)
        
        print(f"[OpenMM] ✓ Уравновешивание завершено")
    
    def run_production(
        self,
        num_steps: int = 500000,
        report_interval: int = 1000,
        trajectory_file: str = "trajectory.dcd",
        log_file: str = "simulation.log"
    ) -> Dict:
        """Production MD."""
        print(f"[OpenMM] Production MD ({num_steps:,} шагов)...")
        
        # Репортеры
        self.simulation.reporters.append(
            app.DCDReporter(trajectory_file, report_interval)
        )
        
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
                speed=True
            )
        )
        
        # Запуск
        self.simulation.step(num_steps)
        
        print(f"[OpenMM] ✓ Симуляция завершена")
        print(f"          Траектория: {trajectory_file}")
        
        return {
            "trajectory_file": trajectory_file,
            "num_steps": num_steps
        }

from openmm import *
from openmm.app import *
from openmm.unit import *
import sys
from sys import stdout

print("Preparing input files...")

platform = Platform.getPlatform('CUDA')
output=sys.argv[2]
protein = PDBFile(str(sys.argv[1]))
prot_top = protein.topology
ff = ForceField('amber19-all.xml', 'amber19/tip3pfb.xml')
modeller = Modeller(protein.topology, protein.positions)

print("Adding solvent...")

modeller.addSolvent(ff, ionicStrength = 0.15 * molar, positiveIon = 'Na+', negativeIon = 'Cl-', padding = 10.0 * angstrom, boxShape = 'dodecahedron', model = 'tip3p', neutralize=True)

print("Creating system...")

system = ff.createSystem(modeller.topology, nonbondedMethod = PME, nonbondedCutoff = 1.0 * nanometer, constraints = HBonds, rigidWater=True)

print("Implementing an integrator...")

integrator = LangevinMiddleIntegrator(310 * kelvin, 1/picoseconds, 0.004 * picoseconds)

print("Setting positions...")

simulation = Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy(tolerance=10.0 * kilojoules_per_mole / nanometer)
simulation.reporters.append(PDBReporter(str(output), 2000, enforcePeriodicBox=False))
simulation.reporters.append(StateDataReporter(stdout, 2000, step=True, potentialEnergy=True, temperature=True))

print("Minimizing...")

simulation.step(2000)
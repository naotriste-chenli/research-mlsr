# << Energy Minimization >> Based from OpenMM documentation: https://docs.openmm.org/latest/userguide/

from openmm import *
from openmm.unit import *
from openmm.app import *
import sys
from sys import stdout

platform = Platform.getPlatform('CUDA')

protein = PDBFile(str(sys.argv[1]))
ff = ForceField('amber14-all.xml', 'implicit/gbn2.xml')
modeller = Modeller(protein.topology, protein.positions)

print("Minimizing: " + str(sys.argv[1]))

system = ff.createSystem(modeller.topology, nonbondedMethod = NoCutoff, constraints = HBonds)
integrator = LangevinMiddleIntegrator(311.15 * kelvin, 1/picosecond, 0.001 * picosecond)

# << Harmonic restraint >> Adapted from OpenMM documents for custom forces: https://opensimconfluence.atlassian.net/wiki/spaces/OpenSim/overview

hookes_restraint = CustomExternalForce('0.5 * k * ((x - xi)^2 + (y - yi)^2 + (z - zi)^2)')
hookes_restraint.addGlobalParameter('k', 1000)
hookes_restraint.addPerParticleParameter('xi')
hookes_restraint.addPerParticleParameter('yi')
hookes_restraint.addPerParticleParameter('zi')

for atoms in modeller.topology.atoms():
	if atoms.element.symbol != 'H':
		coor = modeller.positions[atoms.index]
		hookes_restraint.addParticle(atoms.index, [coor.x, coor.y, coor.z])

system.addForce(hookes_restraint)

simulation = Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(modeller.positions)

simulation.minimizeEnergy(tolerance = 1.0 * kilojoules_per_mole / nanometer)
simulation.reporters.append(PDBReporter(str(sys.argv[2]), 2500))

simulation.step(2500)

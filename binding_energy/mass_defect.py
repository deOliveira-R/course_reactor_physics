#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 31 05:50:43 2021

@author: rodrigo

Source of constants: 2018 CODATA recommended values
Source of mass: https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
                https://www.nndc.bnl.gov/nudat3/
"""

import numpy as np

kg_per_Da = 1.66053906660E-27
J_per_Da = 1.49241808560E-10
MeV_per_Da = 931.49410242
J_per_MeV = 1.602176634E-13


def m_to_s(m):
    return m*60


def h_to_s(h):
    return h*m_to_s(60)


def d_to_s(d):
    return d*h_to_s(24)


def y_to_s(y):
    return y*d_to_s(365.25)


def mean_t_to_lambda(meant_t):
    return 1/meant_t


def mean_t_to_half_life(mean_t):
    return -np.log(1/2)*mean_t


def half_life_to_lambda(half_life):
    return -np.log(1/2)/half_life


class Particle:
    def __init__(self, name, mass, stable=True, mean_lifetime=None):
        self.name = name
        self.mass = mass
        self.stable = stable
        self.mean_lifetime = mean_lifetime


# Eventually we add Hadron and Fermion class with properties for inheritance
# but it's not relevant at the moment
class Baryon(Particle):
    def __init__(self, A, Z, name, mass, stable=True, mean_lifetime=None):
        Particle.__init__(self, name, mass, stable, mean_lifetime)
        self.Z = Z
        self.A = A
        self.N = self.A - self.Z

    def mass_excess(self):
        return self.mass - 1


class Lepton(Particle):
    def __init__(self, name, mass, stable=True, mean_lifetime=None):
        Particle.__init__(self, name, mass, stable, mean_lifetime)
        # add then some other stuff specific to Leptons


class Nuclide:
    def __init__(self, A, Z, E, name, mass, half_life=None, daughter=None, fy=None):
        self.A = A
        self.Z = Z
        self.E = E
        self.N = self.A - self.Z
        self.name = name
        self.mass = mass

        self.half_life = half_life

        # both attributes should be replaced by decay modes with branch
        self.daughter = daughter
        self.fy = fy

    def predicted_mass(self):
        return self.Z * H1.mass + self.N * neutron.mass

    def mass_defect(self):
        return self.predicted_mass() - self.mass

    def mass_defect_nucleon(self):
        return self.mass_defect() / self.A

    def binding_energy(self):
        return self.mass_defect() * MeV_per_Da

    def binding_energy_nucleon(self):
        return self.mass_defect_nucleon() * MeV_per_Da

    def mass_excess(self):
        return self.mass - self.A

    def report_data(self):
        print(f"{self.name.capitalize()} data:")
        print("Predicted mass:", self.predicted_mass())
        print("Mass defect:", self.mass_defect())
        print("MDN:", self.mass_defect_nucleon())
        print("BEN:", self.binding_energy_nucleon())
        print("Mass excess:", self.mass_excess(), '\n')

    def decay_mass(self):
        # getattr(self, 'name')
        return self.mass - globals()[self.daughter].mass

    def Q_beta(self):
        return self.decay_mass() * MeV_per_Da


proton = Baryon(1, 1, 'proton', 1.007276466621)
neutron = Baryon(1, 0, 'neutron', 1.00866491588, stable=False, mean_lifetime=879.4)

electron = Lepton('electron', 0.000548579909070)

H1 = Nuclide(1, 1, 'H', 'hydrogen', 1.00782503223)

U238 = Nuclide(238, 92, 'U', 'uranium', 238.0507883)
U236 = Nuclide(236, 92, 'U', 'uranium', 236.045568)
U235 = Nuclide(235, 92, 'U', 'uranium', 235.0439299)

Rn222 = Nuclide(222, 86, 'Rn', 'radon', 222.0175763)

Ce140 = Nuclide(140, 58, 'Ce', 'cerium', 139.9054387)

La140 = Nuclide(140, 57, 'La', 'lanthanum', 139.9094776, half_life=d_to_s(1.6781), daughter='Ce140')

Ba141 = Nuclide(141, 56, 'Ba', 'barium', 140.914411, half_life=m_to_s(18.27), daughter='La141')
Ba140 = Nuclide(140, 56, 'Ba', 'barium', 139.910605, half_life=d_to_s(12.752), daughter='La140')

Cs135 = Nuclide(135, 55, 'Cs', 'cesium', 134.9059770)
Cs133 = Nuclide(133, 55, 'Cs', 'cesium', 132.905451933)

Xe135 = Nuclide(135, 54, 'Xe', 'Xenon', 134.907227, half_life=h_to_s(9.14), daughter='Cs135')
Xe131 = Nuclide(131, 54, 'Xe', 'Xenon', 130.9050824)

I135 = Nuclide(135, 53, 'I', 'iodine', 134.910048, half_life=h_to_s(6.57), daughter='Xe135')
I131 = Nuclide(131, 53, 'I', 'iodine', 130.9061246, half_life=d_to_s(8.02070), daughter='Xe131')

Tc99m = Nuclide(99, 43, 'Tc', 'technetium', 98.9062547, half_life=h_to_s(6.01), daughter='Tc99')

Mo99 = Nuclide(99, 42, 'Mo', 'molybdenum', 98.9077119, half_life=h_to_s(65.94), daughter='Tc99m')

Zr96 = Nuclide(96, 40, 'Zr', 'zirconium', 95.9082734)
Zr92 = Nuclide(92, 40, 'Zr', 'zirconium', 91.9050408)

Y92 = Nuclide(92, 39, 'Y', 'yttrium', 91.908949, half_life=h_to_s(3.54), daughter='Zr92')

Sr92 = Nuclide(92, 38, 'Sr', 'strontium', 91.911038, half_life=h_to_s(2.66), daughter='Y92')

Rb92 = Nuclide(92, 37, 'Rb', 'rubidium', 91.919729, half_life=4.492, daughter='Sr92')

Kr92 = Nuclide(92, 36, 'Kr', 'kryptonium', 91.926156, half_life=1.840, daughter='Rb92')

Fe56 = Nuclide(56, 26, 'Fe', 'iron', 55.9349375)


def capture_mass(capturing, captured):
    return (capturing.mass_excess() + neutron.mass_excess()) - captured.mass_excess()


def capture_mass_old(capturing, captured):
    return (capturing.mass + neutron.mass) - captured.mass

def Q_capture(capturing, captured):
    return capture_mass(capturing, captured)*MeV_per_Da


def fission_multiplicity(fissile, product1, product2):
    return (fissile.A + 1) - (product1.A + product2.A)


def fission_mass(fissile, product1, product2):
    """
    Excess mass can be used because a nuclear reaction conserves nucleons, so
    M = Au + deltaMass, but Au is the same in both sides, so only deltas are left

    Using mass excess, we deal with smaller values
    """
    n_released = fission_multiplicity(fissile, product1, product2)

    mass_released = (fissile.mass_excess() + neutron.mass_excess()) \
                   - (product1.mass_excess() + product2.mass_excess() + n_released*neutron.mass_excess())

    return mass_released


def fission_mass_old(fissile, product1, product2):
    n_released = fission_multiplicity(fissile, product1, product2)

    mass_released = (fissile.mass + neutron.mass) \
                   - (product1.mass + product2.mass + n_released*neutron.mass)

    return mass_released


def Q_fission(fissile, product1, product2):
    return fission_mass(fissile, product1, product2)*MeV_per_Da


# These are to demonstrate the impact that ignoring an electron has
# by using the proton, instead of hydrogen, on the mass defect calculation
def predicted_mass_noe(nuc):
    return nuc.Z*proton.mass + nuc.N*neutron.mass


def mass_defect_noe(nuc):
    return predicted_mass_noe(nuc) - nuc.mass


def mass_defect_nucleon_noe(nuc):
    return mass_defect_noe(nuc)/nuc.A


def binding_energy_nucleon_noe(nuc):
    return mass_defect_nucleon_noe(nuc)*MeV_per_Da


def report_data_wrong(nuc):
    print(f"{nuc.name.capitalize()} data (WRONG!):")
    print("Predicted mass:", predicted_mass_noe(nuc))
    print("Mass defect:", mass_defect_noe(nuc))
    print("MDN:", mass_defect_nucleon_noe(nuc))
    print("BEN:", binding_energy_nucleon_noe(nuc), '\n')


report_data_wrong(U235)
U235.report_data()

print("The correct value, given by NuDat3, is 7.5909155 MeV\n")

Fe56.report_data()

print('U235 (n,f) Mo99 + I135')
print('Neutrons released: ', fission_multiplicity(U235, Mo99, I135))
print('Mass converted: ', fission_mass(U235, Mo99, I135))
print('Mass converted (old): ', fission_mass_old(U235, Mo99, I135))
print('fission Q: ', Q_fission(U235, Mo99, I135), '\n')

print('U235 (n,f) Rb92 + Ba141')
print('Neutrons released: ', fission_multiplicity(U235, Rb92, Ba141))
print('Mass converted: ', fission_mass(U235, Rb92, Ba141))
print('fission Q: ', Q_fission(U235, Rb92, Ba141), '\n')

print('U235 (n,gamma) U236')
print('Mass converted: ', capture_mass(U235, U236))
print('Mass converted (old): ', capture_mass_old(U235, U236))
print('capture Q: ', Q_capture(U235, U236), '\n')

print('Kr92 beta decay')
print('Mass converted: ', Kr92.decay_mass())
print('Q beta: ', Kr92.Q_beta(), '\n')

print('Rb92 decay')
print('Mass converted: ', Rb92.decay_mass())
print('Q beta: ', Rb92.Q_beta(), '\n')

print('Sr92 decay')
print('Mass converted: ', Sr92.decay_mass())
print('Energy released: ', Sr92.Q_beta(), '\n')

print('Y92 decay')
print('Mass converted: ', Y92.decay_mass())
print('Energy released: ', Y92.Q_beta(), '\n')

# A 3 GW reactor example
print('A 3 GW reactor')
power = 3E9  # W
average_fission_energy = 200  # MeV
fission_rate = power / (average_fission_energy * J_per_MeV)
print('Produces: ', fission_rate, ' fissions/second')

Da_per_second = power / J_per_Da
kg_per_second = Da_per_second * kg_per_Da
print('Converts: ', Da_per_second, ' Da/second')
print('Converts: ', kg_per_second, ' kg/second')
seconds_per_year = d_to_s(365)
kg_per_year = kg_per_second * seconds_per_year
print('Converts: ', kg_per_year, ' kg/year')

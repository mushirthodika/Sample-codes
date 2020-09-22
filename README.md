# MolSSI
Sample codes for MolSSI application

He_HF.f90 --> Hartree-Fock code for the energy calculation of Helium ground state.

dvr.f90 --> Discrete Variable Representation based on Colbert and Miller's paper for the calculation of eigenfunctions for a given potential. (J. Chem. Phys. 96 (3), 1992). 

findirc.f90 --> This code retrieves the xyz coordinates of geometries along the IRC from a Gaussian09 IRC calculation.

abs_spec-new.f90 --> Generates the absorption spectra from a set of energies and oscillator strengths. This code is a part of our in-house lab code called SArCASM (Simulating Absorption Spectra for Complex And Solvated Molecules).

stab_analysis_555.m --> Code for performing (5,5,5) Generalized Pade approximants analytic continuation. Used in the determination of complex energies for electronic resonances from stabilization calculations.

cap-int.py --> Subroutine for calculating complex absorbing potential matrix elements in the atomic basis. (Reference - J. Chem. Phys. 115, 6853, 2001)

asymmdw_kh.f90 --> Calculates a time-averaged Kramers-Henneberger potential for a given molecular potential. (Reference - Phys. Rev. Lett. 21, 838, 1968)

# Equilibration of coupled system of Phonons and Electrons

This repository contains codes pertaining to the problem of real-time equilibration of Phonons and electrons. We are investigating the dynamics of a coupled electron-
phonon system within Keldysh formalism, where both electrons and phonons evolve in time self-consistently, as opposed to the standard paradigm of one set of 
constituents forming a static bath for the other. 

Experimentally, this is relevant for various physical models. For example, in devices such as bolometers (measures power of incident radiation by noting the 
concomitant temperature rise in the incident material), the general idea is that charged particles impinge and excite electrons in lattice, which subsequently 
thermalise by redistributing received energy. This would involve a chain-reaction b/w electrons and phonons, each exchanging energy with the other while evolving 
towards steady state. Various results in such systems depend crucially on the specificities of this non-equilibrium dynamics. More interestingly, at low 
temperatures, this dynamics might be interesting as other subsidiary phenomena such as superconductivity (e.g. s-wave) involve electron-phonon interaction to enable 
an entirely new phase of matter.  


We hope to study this phenomena in detail by doing an honest calculation sans popular approximations. Codes for the E-Ph thermalisation problem are both in python and Julia.

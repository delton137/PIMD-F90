!-------------------------------------------------------------------------------------
! Copyright (c) 2016 Daniel C. Elton
!
! This software is licensed under The MIT License (MIT)
! Permission is hereby granted, free of charge, to any person obtaining a copy of this
! software and associated documentation files (the "Software"), to deal in the Software
! without restriction, including without limitation the rights to use, copy, modify, merge,
! publish, distribute, sublicense, and/or sell copies of the Software, and to permit
! persons to whom the Software is furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
! BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
! DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!-------------------------------------------------------------------------------------
module MD
 use dans_timer
 use NormalModes, only: EvolveRing
 use pot_mod
 use dans_timer
 use system_mod
 use consts
 use InputOutput
 use force_calc
 use main_stuff

 contains

!-----------------------------------------------------------------
!--------   Standard full PIMD Velocity-Verlet integration ------
!-----------------------------------------------------------------
subroutine PIMD_VV
 use mpi
 implicit none

 if (pid .eq. 0) then
	!Propagate NH chains
	if (BEADTHERMOSTAT) call bead_thermostat

	if (THERMOSTAT)  then
		call Nose_Hoover(s, uk, global_chain_length, vxi_global, tau, delt2, 3*Natoms*Nbeads, temp*Nbeads)
		PPt = PPt*s
	endif

	!update momenta a half step w/ old forces
	PPt = PPt - MASSCON*dRRt*delt2

	!Normal modes stuff
	if (Nbeads .gt. 1) then
		do i = 1, Nwaters
			!Evolve ring a full step
			Call EvolveRing(RRt(:,3*i-2,:), PPt(:,3*i-2,:), Nbeads, massO)
			Call EvolveRing(RRt(:,3*i-1,:), PPt(:,3*i-1,:), Nbeads, massH)
			Call EvolveRing(RRt(:,3*i-0,:), PPt(:,3*i-0,:), Nbeads, massH)
		enddo
	else
		do i = 1,Nwaters
			do k = 1,Nbeads
				RRt(:,3*i-2,k) = RRt(:,3*i-2,k) + imassO*PPt(:,3*i-2,k)*delt
				RRt(:,3*i-1,k) = RRt(:,3*i-1,k) + imassH*PPt(:,3*i-1,k)*delt
				RRt(:,3*i-0,k) = RRt(:,3*i-0,k) + imassH*PPt(:,3*i-0,k)*delt
			enddo
		enddo
	endif

	!calculate centroid positions
	RRc = sum(RRt,3)/Nbeads

	!calculate centroid momenta
	PPc = sum(PPt,3)/Nbeads

	!check PBCs
	call PBCs(RRt, RRc)

 endif!if (pid .eq. 0)

 !call force routine
 call full_bead_forces

 if (pid .eq. 0) then

	!update kinetic energy
	call calc_uk

	!write stuff out if necessary
	call start_timer("WritingOut")
	call write_out
	call stop_timer("WritingOut")

	!update momenta a half step w/ new forces
	PPt = PPt - MASSCON*dRRt*delt2

	!Propagate NH chains
	if (BEADTHERMOSTAT) call bead_thermostat

	call calc_uk

	if (THERMOSTAT) then
		call Nose_Hoover(s, uk, global_chain_length, vxi_global, tau, delt2, 3*Natoms*Nbeads, temp*Nbeads)
		PPt = PPt*s
	endif

	if (BAROSTAT) call Pcouple

 endif!(pid .eq. 0)

end subroutine PIMD_VV


!-----------------------------------------------------------------
!--------  " Monomer PIMD" --------------------------------------
! THIS IS DEPRECATED --- EVERYTHING IS DONE WITH CONTRACTED ROUTINE BELOW
!-----------------------------------------------------------------
subroutine monomer_PIMD
 use system_mod
 use consts
 use main_stuff
 use InputOutput
 use NormalModes
 use mpi
 use force_calc
 implicit none
 double precision :: e1, virialmon, virialcmon, virial_mon_c
 double precision, dimension(3,Natoms,Nbeads)  ::  dRRmon
 double precision, dimension(3,Natoms) :: dRRc
 double precision, dimension(3,3)      :: dr1, r1
 double precision, dimension(3)        :: roh1, roh2

 if (t .eq. 1) dRRmon = 0
 Upott  = 0
 virialmon = 0
 virialcmon = 0

 if (pid .eq. 0) then
	!Propagate NH chains
	if (BEADTHERMOSTAT) call bead_thermostat

	if (THERMOSTAT)  then
		call Nose_Hoover(s, uk, global_chain_length, vxi_global, tau, delt2, 3*Natoms*Nbeads, temp*Nbeads)
		PPt = PPt*s
	endif

	!update momenta a half step w/ old forces
	PPt = PPt - MASSCON*dRRt*delt2

	!Normal modes stuff
	if (Nbeads .gt. 1) then
		do i = 1, Nwaters
			!Evolve ring a full step
			call EvolveRing(RRt(:,3*i-2,:), PPt(:,3*i-2,:), Nbeads, massO)
			call EvolveRing(RRt(:,3*i-1,:), PPt(:,3*i-1,:), Nbeads, massH)
			call EvolveRing(RRt(:,3*i-0,:), PPt(:,3*i-0,:), Nbeads, massH)
		enddo
	else
		do i = 1,Nwaters
			do k = 1,Nbeads
				RRt(:,3*i-2,k) = RRt(:,3*i-2,k) + imassO*PPt(:,3*i-2,k)*delt
				RRt(:,3*i-1,k) = RRt(:,3*i-1,k) + imassH*PPt(:,3*i-1,k)*delt
				RRt(:,3*i-0,k) = RRt(:,3*i-0,k) + imassH*PPt(:,3*i-0,k)*delt
			enddo
		enddo
	endif

	!calculate centroid positions
	RRc = sum(RRt,3)/Nbeads

	!calculate centroid momenta
	PPc = sum(PPt,3)/Nbeads

	!check PBCs
	call PBCs(RRt, RRc)

	!intermolecular force calculation
	call potential(RRc, RRc, Upot, dRRc, virt, virialc, dip_momI, dip_momE, chg, t, BAROSTAT, sys_label)

	!remove the intramolecular (PS) part from the centroid force - it is already applied on the beads
	call remove_centroid_monomer_term(RRc, dRRc, Upot, virial_mon_c)

    Upott(:) = Upot  !potential energy for the ENTIRE system (all images)

    !intramolecular force calculation
	do j = 1, Nbeads
		do iw = 1, Nwaters
			iO=3*iw-2; iH1 = 3*iw-1; iH2=3*iw-0


			!PS potential energy surface calculation
			r1(1:3, 1:3) = RRt(1:3, (/iO, iH1, iH2/), j)

			if (SIESTA_MON_CALC) then
				call siesta_monomer(r1, dr1, e1)
			else
				call pot_nasa(r1, dr1, e1, box, boxi)
			endif

			dRRmon(1:3, (/iO, iH1, iH2/), j) = dr1

			Upott(j) = Upott(j) + e1

			!monomer centroid virial
			roh1 = RRc(1:3, iH1) - RRc(1:3, iO)
			roh1 = roh1 - box*anint(roh1*boxi) !PBC
			roh2 = RRc(1:3, iH2) - RRc(1:3, iO)
			roh2 = roh2 - box*anint(roh2*boxi) !PBC
			virialcmon = virialcmon + dot_product(roh1, dr1(:,2))
			virialcmon = virialcmon + dot_product(roh2, dr1(:,3))

			!monomer normal virial
			roh1 = RRt(1:3, iH1, j) - RRt(1:3, iO, j)
			roh1 = roh1 - box*anint(roh1*boxi) !PBC
			roh2 = RRt(1:3, iH2, j) - RRt(1:3, iO, j)
			roh2 = roh2 - box*anint(roh2*boxi) !PBC
			virialmon = virialmon + dot_product(roh1, dr1(:,2))
			virialmon = virialmon + dot_product(roh2, dr1(:,3))

 		enddo
	enddo

	!update dRRt
	do j = 1, Nbeads
		dRRt(:,:,j) = dRRc + dRRmon(:,:,j)
	enddo

	!update Upot
    Upot    = sum(Upott)    	!potential energy for the ENTIRE system
	virial  = virialmon + virt(1,1) + virt(2,2) + virt(3,3) - virial_mon_c !virial for the ENTIRE system (all images)
	virialc = virialcmon/Nbeads + virialc - virial_mon_c

    call calc_monomer_dip_moments(dip_momIt, RRt)

	!calculate electronic polarization dipoles using TTM method
	if (pot_model .eq. 6) call dip_ttm(RRc, dip_momI, dip_momE, chg, t)

	!add electronic polarization dipoles to monomer dipoles
	do j = 1, Nbeads
		dip_momIt(:,:,j) = dip_momIt(:,:,j) + dip_momE
		dip_momEt(:,:,j) = dip_momE
	enddo

	!update kinetic energy
	call calc_uk

	!write stuff out if necessary
	call start_timer("WritingOut")
	call write_out
	call stop_timer("WritingOut")

	!update momenta a half step w/ new forces
	PPt = PPt - MASSCON*dRRt*delt2

	!Propagate NH chains
	if (BEADTHERMOSTAT) call bead_thermostat

	call calc_uk

	if (THERMOSTAT) then
		call Nose_Hoover(s, uk, global_chain_length, vxi_global, tau, delt2, 3*Natoms*Nbeads, temp*Nbeads)
		PPt = PPt*s
	endif

	if (BAROSTAT) call Pcouple

 endif!(pid .eq. 0)

end subroutine monomer_PIMD


!--------------------------------------------------------------------------------------------
!------------- Contracted MD with intermolecular forces on monomer with multiple time step
!-- This subroutine performs evaluates only the intramolecular (fast) forces on the beads and
!-- evaluates the intermolecular (slow) forces on the centroid. It also uses a multiple timestep method
!-- (see Tuckerman, pg 118). The momentum and positions will be updated with the intramolecular forces
!-- every intra_timesteps times. For instance if the 'outer timestep' (normal timestep) is .5 ps and intra_timesteps = 5
!-- then the inner timestep is .1 ps
!--------------------------------------------------------------------------------------------
subroutine contracted_MD
 use main_stuff
 use NormalModes
 use pot_mod
 use dans_timer
 use system_mod
 use consts
 use InputOutput
 use force_calc
 Implicit None
 double precision :: e1, virialmon, virialcmon, virial_mon_c
 double precision, dimension(3,Natoms,Nbeads)  ::  dRRfast
 double precision, dimension(3,Natoms) :: dRRc
 double precision, dimension(3,3)      :: dr1, r1
 double precision, dimension(3)        :: roh1, roh2
 integer :: tintra, iM

 if (t .eq. 1) dRRfast =0
 Upott  = 0
 virialmon = 0
 virialcmon = 0

 if (pid .eq. 0) then
	!Propagate NH chains
	if (BEADTHERMOSTAT) call bead_thermostat

	if (THERMOSTAT)  then
		call Nose_Hoover(s, uk, global_chain_length, vxi_global, tau, delt2, 3*Natoms*Nbeads, temp*Nbeads)
		PPt = PPt*s
	endif

	!update momenta a half step w/ old forces
	PPt = PPt - MASSCON*dRRt*delt2

	!calculate centroid positions
	RRc = sum(RRt,3)/Nbeads
		!calculate centroid momenta
	PPc = sum(PPt,3)/Nbeads

	!check PBCs
	call PBCs(RRt, RRc)

	call start_timer("MonomerPIMD")

	!---  intramolecular (fast) forces -------------------------------------------------
    do tintra = 1, intra_timesteps

		!update momenta with fast forces (delt2fast = deltfast/2)
		PPt = PPt - MASSCON*dRRfast*delt2fast

		!update positions with fast forces
		if (Nbeads .gt. 1) then
			do i = 1, Nwaters
				Call EvolveRing(RRt(:,3*i-2,:), PPt(:,3*i-2,:), Nbeads, massO)
				Call EvolveRing(RRt(:,3*i-1,:), PPt(:,3*i-1,:), Nbeads, massH)
				Call EvolveRing(RRt(:,3*i-0,:), PPt(:,3*i-0,:), Nbeads, massH)
			enddo
		else
			do i = 1,Nwaters
				do k = 1,Nbeads
					RRt(:,3*i-2,k) = RRt(:,3*i-2,k) + imassO*PPt(:,3*i-2,k)*deltfast
					RRt(:,3*i-1,k) = RRt(:,3*i-1,k) + imassH*PPt(:,3*i-1,k)*deltfast
					RRt(:,3*i-0,k) = RRt(:,3*i-0,k) + imassH*PPt(:,3*i-0,k)*deltfast
				enddo
			enddo
		endif

		!update fast forces (intramolecular forces)
		!masternode calcuates the intramolecular forces, puts them in dRRfast
		virialmon = 0
		virialcmon = 0

		do j = 1, Nbeads
			do iw = 1, Nwaters
				iO=3*iw-2; iH1 = 3*iw-1; iH2=3*iw-0

				r1(1:3, 1:3) = RRt(1:3, (/iO, iH1, iH2/), j)

				call pot_nasa(r1, dr1, e1, box, boxi)

				dRRfast(1:3, (/iO, iH1, iH2/), j) = dr1

				!if last timestep in loop update monomer energy
				!and calculate dipole moments using dip. mom. surface
				if (tintra .eq. intra_timesteps) then
					Upott(j) = Upott(j) + e1
					!monomer centroid virial first
					roh1 = RRc(1:3, iH1) - RRc(1:3, iO)
					roh1 = roh1 - box*anint(roh1*boxi) !PBC
					roh2 = RRc(1:3, iH2) - RRc(1:3, iO)
					roh2 = roh2 - box*anint(roh2*boxi) !PBC
					virialcmon = virialcmon + dot_product(roh1, dr1(:,2))
					virialcmon = virialcmon + dot_product(roh2, dr1(:,3))
					!monomer normal virial
					roh1 = RRt(1:3, iH1, j) - RRt(1:3, iO, j)
					roh1 = roh1 - box*anint(roh1*boxi) !PBC
					roh2 = RRt(1:3, iH2, j) - RRt(1:3, iO, j)
					roh2 = roh2 - box*anint(roh2*boxi) !PBC
					virialmon = virialmon + dot_product(roh1, dr1(:,2))
					virialmon = virialmon + dot_product(roh2, dr1(:,3))
				endif
			enddo
		enddo

		!update momenta with fast forces
		PPt = PPt - MASSCON*dRRfast*delt2fast

	enddo!tintra  = 1..
	!---  end intramolecular (fast) forces ----------------------------------------------

	!calculate centroid positions
	RRc = sum(RRt,3)/Nbeads

	!calculate centroid momenta
	PPc = sum(PPt,3)/Nbeads

	!check PBCs
	call PBCs(RRt, RRc)

	call stop_timer("MonomerPIMD")

	!intermolecular force calculation
	call potential(RRc, RRc, Upot, dRRc, virt, virialc, dip_momI, dip_momE, chg, t, BAROSTAT, sys_label)

	!remove the intramolecular (PS) part from the centroid force - it is already applied on the beads
	call remove_centroid_monomer_term(RRc, dRRc, Upot, virial_mon_c)

    Upott(:) = Upott(:) + Upot  !potential energy for the ENTIRE system (all images)

	!update dRRt
	do j = 1, Nbeads
		dRRt(:,:,j) = dRRc
	enddo

	call calc_monomer_dip_moments(dip_momIt, RRt)

!	!calculate electronic polarization dipoles using TTM method
	if (pot_model .eq. 6) call dip_ttm(RRc, dip_momE, chg, t)

    !write(*,*) dip_momE

  !add electronic polarization dipoles to monomer dipoles
	do j = 1, Nbeads
		dip_momIt(:,:,j) = dip_momIt(:,:,j) + dip_momE
		dip_momEt(:,:,j) = dip_momE
	enddo

	!update Upot, virial and virialc
    Upot    = sum(Upott)    	!potential energy for the ENTIRE system
	virial  = virialmon + virt(1,1) + virt(2,2) + virt(3,3) - virial_mon_c !virial for the ENTIRE system (all images)
	virialc = virialcmon/Nbeads + virialc - virial_mon_c
	!update kinetic energy
	call calc_uk

	!write stuff out if necessary
	call start_timer("WritingOut")
	call write_out
	call stop_timer("WritingOut")

	!update momenta a half step w/ slow forces
	PPt = PPt - MASSCON*dRRt*delt2

	!Propagate NH chains
	if (BEADTHERMOSTAT) call bead_thermostat

	call calc_uk

	if (THERMOSTAT) then
		call Nose_Hoover(s, uk, global_chain_length, vxi_global, tau, delt2, 3*Natoms*Nbeads, temp*Nbeads)
		PPt = PPt*s
	endif

	if (BAROSTAT) call Pcouple

endif!(pid .eq. 0)

end subroutine contracted_MD


!--------------------------------------------------------------------------------------------
!-- Remove the intramolecular (PS) contribution from the slow (centroid) force
!--
!-- In the contracted / monomer-PIMD schemes the intramolecular Partridge-Schwenke (PS)
!-- potential is evaluated on every bead (the fast force). The slow force evaluated at the
!-- centroid must therefore contain ONLY the intermolecular part of the potential. If the
!-- intramolecular part is left in, the intramolecular force constant is counted twice and
!-- every intramolecular frequency is shifted up by sqrt(2): 1611 -> 2294 cm^-1 (bend) and
!-- 3922 -> 5550 cm^-1 (stretch) for a single water. These were the unexplained extra peaks
!-- in Figs. 11-12 of Elton, Fritz & Fernandez-Serra, PCCP 21, 409 (2019).
!--
!-- The TTM potentials take care of this internally (pot_ttm skips pot_nasa when CONTRACTION
!-- is set). SIESTA returns the complete potential it was asked for - plain DFT, or DFT with
!-- the monomer correction  E_DFT + E_PS - E_fit  of the modified SIESTA - so the PS term is
!-- removed here:
!--                     V_slow(c) = V_SIESTA(c) - V_PS(c)
!-- With the monomer-corrected SIESTA this gives V_slow = E_DFT - E_fit, i.e. exactly the
!-- intermolecular part of the corrected potential, and with Nbeads = 1 the scheme reduces
!-- to classical MD on the SIESTA potential, as it should.
!-- NOTE: SIESTA must return the full potential (including its PS monomer term). Do not use
!-- a SIESTA build whose monomer correction returns E_DFT - E_fit only; that would remove the
!-- PS term twice.
!--------------------------------------------------------------------------------------------
subroutine remove_centroid_monomer_term(RRc, dRRc, Upot, virial_mon_c)
 Implicit None
 double precision, dimension(3,Natoms), intent(in)    :: RRc
 double precision, dimension(3,Natoms), intent(inout) :: dRRc
 double precision, intent(inout) :: Upot
 double precision, intent(out)   :: virial_mon_c
 double precision, dimension(3,3) :: r1, dr1
 double precision, dimension(3)   :: roh1, roh2
 double precision :: e1
 integer :: iw, jO, jH1, jH2

 virial_mon_c = 0d0

 if (pot_model .ne. 6) return  !TTM potentials already exclude the monomer term when CONTRACTION is set

 do iw = 1, Nwaters
	jO = 3*iw-2; jH1 = 3*iw-1; jH2 = 3*iw

	r1(1:3, 1:3) = RRc(1:3, (/jO, jH1, jH2/))

	call pot_nasa(r1, dr1, e1, box, boxi)

	dRRc(1:3, (/jO, jH1, jH2/)) = dRRc(1:3, (/jO, jH1, jH2/)) - dr1
	Upot = Upot - e1

	!corresponding (centroid) virial contribution
	roh1 = r1(:,2) - r1(:,1)
	roh1 = roh1 - box*anint(roh1*boxi) !PBC
	roh2 = r1(:,3) - r1(:,1)
	roh2 = roh2 - box*anint(roh2*boxi) !PBC
	virial_mon_c = virial_mon_c + dot_product(roh1, dr1(:,2)) + dot_product(roh2, dr1(:,3))
 enddo

end subroutine remove_centroid_monomer_term








end module MD

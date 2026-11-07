!---------------------------------------------------------------------------
! MOCK of the fsiesta pipe interface, for testing PIMD-F90's SIESTA code path
! without SIESTA.  It emulates the monomer-corrected SIESTA used in the 2019
! PCCP paper on a set of NON-interacting water molecules:
!     E_SIESTA = sum_mol [ E_DFT + (E_PS - E_fit) ]  with  E_DFT ~= E_fit
!              = sum_mol E_PS(mol)
! so it returns the Partridge-Schwenke monomer energy/forces (in eV, eV/Ang),
! exactly what the real monomer-corrected SIESTA returns for a lone molecule.
! Like the real SIESTA it receives the RAW coordinates; it counts how many
! molecules arrive split across the box edge (|x_H - x_O| > L/2 in any
! direction), which is what tripped the 2.5 eV gate in the real code.
!---------------------------------------------------------------------------
module fsiesta
  implicit none
  integer, parameter :: dp = kind(1.d0)
  double precision, parameter :: EV2KCAL = 23.0600d0
  integer, save :: ncalls = 0, nmol = 0, nsplit = 0
  public
contains
  subroutine siesta_launch(siesta_cmd, label, nnodes, mpi_comm, mpi_launcher)
    character(len=*), intent(in) :: siesta_cmd, label
    integer, optional, intent(in) :: nnodes, mpi_comm
    character(len=*), intent(in), optional :: mpi_launcher
    write(*,*) 'MOCK SIESTA: launched, label = ', trim(label)
  end subroutine siesta_launch

  subroutine siesta_units(length, energy)
    character(len=*), intent(in) :: length, energy
  end subroutine siesta_units

  subroutine siesta_forces(label, na, xa, cell, energy, fa, stress)
    character(len=*),   intent(in) :: label
    integer,            intent(in) :: na
    real(dp),           intent(in) :: xa(3,na)
    real(dp), optional, intent(in) :: cell(3,3)
    real(dp), optional, intent(out):: energy
    real(dp), optional, intent(out):: fa(3,na)
    real(dp), optional, intent(out):: stress(3,3)
    double precision, dimension(3,3) :: r1, dr1
    double precision, dimension(3) :: box, boxi, d
    double precision :: e1
    integer :: im, ih
    box = (/ cell(1,1), cell(2,2), cell(3,3) /)
    boxi = 1d0/box
    energy = 0d0
    fa = 0d0
    ncalls = ncalls + 1
    do im = 1, na/3
       r1(:,:) = xa(:, 3*im-2:3*im)
       nmol = nmol + 1
       do ih = 2, 3
          d = r1(:,ih) - r1(:,1)
          if (any(abs(d) .gt. 0.5d0*box)) then
             nsplit = nsplit + 1
             exit
          endif
       enddo
       call pot_nasa(r1, dr1, e1, box, boxi)   ! pot_nasa applies the minimum image itself
       energy = energy + e1/EV2KCAL
       fa(:, 3*im-2:3*im) = -dr1/EV2KCAL
    enddo
    if (present(stress)) stress = 0d0
  end subroutine siesta_forces

  subroutine siesta_quit(label)
    character(len=*), intent(in) :: label
    write(*,'(a,i10,a,i10,a,i10)') ' MOCK SIESTA: calls =', ncalls, '  molecules seen =', nmol, &
         '  molecules arriving SPLIT across the box edge =', nsplit
    open(unit=99, file='mock_siesta_split_count.dat', status='replace')
    write(99,*) ncalls, nmol, nsplit
    close(99)
  end subroutine siesta_quit
end module fsiesta

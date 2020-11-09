module single_patt_mod
implicit none
private
type mol_patt
    private
    integer(4),public :: n = 0
    integer(4),public,allocatable :: patt(:)
  !contains
end type

interface mol_patt
    module procedure :: init_mol_patt
end interface

public :: mol_patt

contains
!----------------------------------------------------
type(mol_patt) function init_mol_patt(patt) result(res)
    integer(4),intent(in) :: patt(:)
    integer(4) :: ps
    ps = size(patt)
    res%n = ps
    allocate(res%patt(ps))
    res%patt = patt
end function
end module
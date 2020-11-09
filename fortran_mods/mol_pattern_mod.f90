module mol_pattern_mod
use single_patt_mod
implicit none

private

type mol_patts
    private
    integer(4),public :: n = 0
    type(mol_patt),public,allocatable :: patts(:)
  contains
    procedure,public :: find_patts
    procedure,public :: add_patt
end type

public :: mol_patts

contains
!----------------------------------------------------
subroutine add_patt(this, patt)
    class(mol_patts) :: this
    integer(4) :: patt(:), i
    type(mol_patt),allocatable :: patt_tmp(:)
    
    if(allocated(this%patts)) then
        do i=1,this%n
            if(size(patt) == this%patts(i)%n) then
                if( all(patt(:) == this%patts(i)%patt(:)) ) return
            endif
        enddo
        allocate(patt_tmp, source=this%patts)
        deallocate(this%patts)
        i = size(this%patts) + 1
        allocate(this%patts(i))
        this%patts(:i-1) = patt_tmp
        deallocate(patt_tmp)
    else
        i = 1
        allocate(this%patts(i))
    endif

    this%patts(i) = mol_patt(patt)
    write(*,'(*(g0,1x))') "this pattern is added:",patt
    this%n = this%n + 1
end subroutine
!----------------------------------------------------
subroutine find_patts(this, tps, patt, opt)
    class(mol_patts) :: this
    integer(4),intent(in) :: patt(:)
    character,intent(in) :: opt(:)
    integer(4),intent(in) :: tps(:)
    integer(4) :: ntp, np
    integer(4) :: it, ip, it_p
    logical :: one

    ntp = size(tps)
    np = size(patt)

    if(np /= size(opt)) then
        write(*,*) "np /= size(opt)"
        return
    endif

    it = 1
    do while(it <= ntp)
        ip = 1
        it_p = 0
        one = .false.
        do while(ip <= np)
            select case(opt(ip))
            case(".")
                if(tps(it+it_p) /= patt(ip)) exit
                ip = ip + 1
                it_p = it_p + 1
            case("+")
                if(tps(it+it_p) == patt(ip)) then
                    if(.not. one) one = .true.
                    it_p = it_p + 1
                    cycle
                else
                    ip = ip + 1
                    if(one) cycle
                    exit
                endif
            case default
                stop "incorect opt in find_patts subroutine"
            end select
        enddo
        if(ip > np) then
            call this%add_patt(tps(it:it+it_p-1))
            it = it + max(it_p,1)
        else
            it = it + 1
        endif
    enddo

        
end subroutine
end module
module config_mod
implicit none

private

type config
    private
    integer(4),public :: n = 0
    integer(4),public :: n_mols = 0
    integer(4),public :: n_bonds = 0
    integer(4),public :: n_angles = 0
    integer(4),public,allocatable :: tp(:),mol(:)
    integer(4),public,allocatable :: tp_bonds(:), bonds(:,:)
    integer(4),public,allocatable :: tp_angles(:), angles(:,:)
    integer(4) :: iter_id, iter_mol, mol_last

    real(8),public,allocatable :: r(:,:), q(:), rad(:)
    real(8),public :: Lmin(3), Lmax(3), Lbox(3)

    logical,public,allocatable :: is_any(:)

  contains
    procedure,public :: get_MolSize

    procedure,public :: set_box_size
    procedure,public :: set_n
    
    procedure,public :: add_particle
    procedure,public :: add_particle_unsafe
    procedure,public :: add_radius_id_unsafe
    procedure,public :: add_charge
    procedure,public :: add_nmol
    procedure,public :: add_bond
    procedure,public :: add_angle

    procedure,public :: iter_reset
    procedure,public :: iter_NextMol
    procedure,public :: iter_hasNextMol

    procedure,public :: print_config
    procedure,public :: print_xyz

    procedure,public :: this_is_any
end type

public :: config

contains
!----------------------------------------------------
subroutine print_config(this, out_file)
    class(config) :: this
    integer(4) :: i, io
    character(len=*),optional,intent(in) :: out_file
    character(len=:),allocatable :: ufmt, fname
    
    ufmt = "(*(g0,1x))"

    if(present(out_file)) then
        fname = out_file
    else
        fname = "cg_config.data"
    endif

    write(*,*) "print to: ",fname

    open(newunit=io,file=fname)
    write(io,*) ""
    write(io,ufmt) this%n,"atoms"
    write(io,ufmt) this%n_bonds,"bonds"
    write(io,ufmt) this%n_angles,"angles"
    write(io,*) ""
    write(io,ufmt) maxval(this%tp),"atom types"
    if(allocated(this%tp_bonds)) write(io,ufmt) maxval(this%tp_bonds),"bond types"
    if(allocated(this%tp_angles)) write(io,ufmt) maxval(this%tp_angles),"angle types"
    write(io,*) ""
    write(io,ufmt) this%Lmin(1), this%Lmax(1),"xlo xhi"
    write(io,ufmt) this%Lmin(2), this%Lmax(2),"ylo yhi"
    write(io,ufmt) this%Lmin(3), this%Lmax(3),"zlo zhi"
    write(io,*) ""
    write(io,ufmt) "Atoms"
    write(io,*) ""
    do i = 1, this%n
        write(io,ufmt) i,this%mol(i),this%tp(i),this%q(i),this%r(:,i)
    enddo
    write(io,*) ""
    if(this%n_bonds /= 0) write(io,ufmt) "Bonds"
    write(io,*) ""
    do i=1,this%n_bonds
        write(io,ufmt) i,this%tp_bonds(i),this%bonds(:,i)
    enddo
    write(io,*) ""
    if(this%n_angles /= 0) write(io,ufmt) "Angles"
    write(io,*) ""
    do i=1,this%n_angles
        write(io,ufmt) i,this%tp_angles(i),this%angles(:,i)
    enddo
    close(io)
end subroutine
!----------------------------------------------------
subroutine print_xyz(this, out_file)
    class(config) :: this
    integer(4) :: i, io
    character(len=*),optional,intent(in) :: out_file
    character(len=:),allocatable :: ufmt, fname
    
    ufmt = "(*(g0,1x))"

    if(present(out_file)) then
        fname = out_file
    else
        fname = "config.xyz"
    endif

    write(*,*) "print to: ",fname

    open(newunit=io,file=fname)
    write(io,ufmt) this%n
    write(io,ufmt) "class(config). print_xyz output"
    do i = 1, this%n
        write(io,ufmt) this%tp(i),this%r(:,i)
    enddo
    close(io)
end subroutine
!----------------------------------------------------
subroutine add_charge(this,q,tp)
    class(config) :: this
    real(8),intent(in) :: q
    integer(4),intent(in) :: tp
    integer(4) :: i
    do concurrent(i = 1:this%n)
        if(this%tp(i) == tp) then
            this%q(i) = q
        endif
    enddo
end subroutine
!----------------------------------------------------
subroutine add_particle(this,rin, tpin, id)
    class(config) :: this
    real(8),intent(in) :: rin(3)
    integer(4),intent(in) :: tpin
    integer(4),optional,intent(in) :: id
    real(8),allocatable :: r_tmp(:,:),q_tmp(:), rad_tmp(:)
    integer(4),allocatable :: tp_tmp(:),mol_tmp(:)
    logical,allocatable :: is_any_tmp(:)
    integer(4) :: n_old, this_id

    n_old = this%n
    if(present(id)) then
        if(id > this%n) this%n = id
        this_id = id
    else
        this%n = this%n + 1
        this_id = this%n
    endif

    if(allocated(this%r)) then
        allocate(r_tmp,source=this%r)
        deallocate(this%r)
    endif
    if(allocated(this%tp)) then
        allocate(tp_tmp,source=this%tp)
        deallocate(this%tp)
    endif
    if(allocated(this%mol)) then
        allocate(mol_tmp,source=this%mol)
        deallocate(this%mol)
    endif
    if(allocated(this%rad)) then
        allocate(rad_tmp,source=this%rad)
        deallocate(this%rad)
    endif
    if(allocated(this%is_any)) then
        allocate(is_any_tmp,source=this%is_any)
        deallocate(this%is_any)
    endif
    if(allocated(this%q)) then
        allocate(q_tmp,source=this%q)
        deallocate(this%q)
    endif

    if(.not. allocated(this%r)) allocate(this%r(3,this%n))
    if(.not. allocated(this%q)) allocate(this%q(this%n))
    if(.not. allocated(this%tp)) allocate(this%tp(this%n))
    if(.not. allocated(this%mol)) allocate(this%mol(this%n))
    if(.not. allocated(this%rad)) allocate(this%rad(this%n))
    if(.not. allocated(this%is_any)) allocate(this%is_any(this%n))

    if(allocated(r_tmp)) this%r(:,:n_old) = r_tmp
    if(allocated(q_tmp)) this%q(:n_old) = q_tmp
    if(allocated(tp_tmp)) this%tp(:n_old) = tp_tmp
    if(allocated(mol_tmp)) this%mol(:n_old) = mol_tmp
    if(allocated(rad_tmp)) this%rad(:n_old) = rad_tmp
    if(allocated(is_any_tmp)) this%is_any(:n_old) = is_any_tmp

    this%r(:,this_id) = rin
    this%tp(this_id) = tpin
    this%mol(this_id) = this_id
    this%rad(this_id) = this_id
    this%is_any(this_id) = .false.
    this%q(this_id) = 0.d0

    if(allocated(r_tmp)) deallocate(r_tmp)
    if(allocated(q_tmp)) deallocate(q_tmp)
    if(allocated(tp_tmp)) deallocate(tp_tmp)
    if(allocated(mol_tmp)) deallocate(mol_tmp)
    if(allocated(rad_tmp)) deallocate(rad_tmp)
    if(allocated(is_any_tmp)) deallocate(is_any_tmp)

end subroutine
!----------------------------------------------------
!----------------------------------------------------
subroutine set_n(this,newn)
    class(config) :: this
    integer(4),intent(in) :: newn

    this%n = newn

    if(allocated(this%r)) deallocate(this%r)
    if(allocated(this%tp)) deallocate(this%tp)
    if(allocated(this%mol)) deallocate(this%mol)
    if(allocated(this%rad)) deallocate(this%rad)
    if(allocated(this%is_any)) deallocate(this%is_any)
    if(allocated(this%q)) deallocate(this%q)

    allocate(this%r(3,this%n))
    allocate(this%q(this%n))
    allocate(this%tp(this%n))
    allocate(this%mol(this%n))
    allocate(this%rad(this%n))
    allocate(this%is_any(this%n))

end subroutine
!----------------------------------------------------
subroutine add_particle_unsafe(this,rin, tpin, id)
    class(config) :: this
    real(8),intent(in) :: rin(3)
    integer(4),intent(in) :: tpin
    integer(4),intent(in) :: id

    if(id > this%n) stop

    this%r(:,id) = rin
    this%tp(id) = tpin
    this%mol(id) = id
    this%is_any(id) = .false.
    this%q(id) = 0.d0
end subroutine
!----------------------------------------------------
subroutine add_radius_id_unsafe(this,radin, id)
    class(config) :: this
    real(8),intent(in) :: radin
    integer(4),intent(in) :: id

    if(id > this%n) stop

    this%rad(id) = radin
end subroutine
!----------------------------------------------------
subroutine add_bond(this, i_in, j_in, tpb, id)
    class(config) :: this
    integer(4),intent(in) :: i_in,j_in,tpb
    integer(4),optional,intent(in) :: id

    integer(4),allocatable :: bonds_tmp(:,:)
    integer(4),allocatable :: tpb_tmp(:)
    integer(4) :: n_old, this_id

    n_old = this%n_bonds
    if(present(id)) then
        if(id > this%n) this%n_bonds = id
        this_id = id
    else
        this%n_bonds = this%n_bonds + 1
        this_id = this%n_bonds
    endif

    if(allocated(this%bonds)) then
        allocate(bonds_tmp,source=this%bonds)
        deallocate(this%bonds)
    endif

    if(allocated(this%tp_bonds)) then
        allocate(tpb_tmp,source=this%tp_bonds)
        deallocate(this%tp_bonds)
    endif

    if(.not. allocated(this%bonds)) allocate(this%bonds(2,this%n_bonds))
    if(.not. allocated(this%tp_bonds)) allocate(this%tp_bonds(this%n_bonds))

    if(allocated(bonds_tmp)) this%bonds(:,:n_old) = bonds_tmp
    if(allocated(tpb_tmp)) this%tp_bonds(:n_old) = tpb_tmp

    this%bonds(:,this_id) = [i_in, j_in]
    this%tp_bonds(this_id) = tpb

    if(allocated(bonds_tmp)) deallocate(bonds_tmp)
    if(allocated(tpb_tmp)) deallocate(tpb_tmp)

end subroutine
!----------------------------------------------------
subroutine add_angle(this, i_in, j_in, k_in, tpa, id)
    class(config) :: this
    integer(4),intent(in) :: i_in,j_in,k_in,tpa
    integer(4),optional,intent(in) :: id

    integer(4),allocatable :: angles_tmp(:,:), tpa_tmp(:)
    integer(4) :: n_old, this_id

    n_old = this%n_angles
    if(present(id)) then
        if(id > this%n) this%n_angles = id
        this_id = id
    else
        this%n_angles = this%n_angles + 1
        this_id = this%n_angles
    endif

    if(allocated(this%angles)) then
        allocate(angles_tmp,source=this%angles)
        deallocate(this%angles)
    endif

    if(allocated(this%tp_angles)) then
        allocate(tpa_tmp,source=this%tp_angles)
        deallocate(this%tp_angles)
    endif

    if(.not. allocated(this%angles)) allocate(this%angles(3,this%n_angles))
    if(.not. allocated(this%tp_angles)) allocate(this%tp_angles(this%n_angles))

    if(allocated(angles_tmp)) this%angles(:,:n_old) = angles_tmp
    if(allocated(tpa_tmp)) this%tp_angles(:n_old) = tpa_tmp

    this%angles(:,this_id) = [i_in, j_in, k_in]
    this%tp_angles(this_id) = tpa

    if(allocated(angles_tmp)) deallocate(angles_tmp)
    if(allocated(tpa_tmp)) deallocate(tpa_tmp)

end subroutine
!----------------------------------------------------
!subroutine add_nmol(this)
!    class(config) :: this
!    this%n_mols = this%n_mols + 1
!end subroutine
!----------------------------------------------------
integer(4) function add_nmol(this)
    class(config) :: this
    this%n_mols = this%n_mols + 1
    add_nmol = this%n_mols
end function
!----------------------------------------------------
subroutine iter_reset(this)
    class(config) :: this
    this%iter_mol = 0
    this%iter_id = 0
end subroutine
!----------------------------------------------------
subroutine this_is_any(this,id)
    class(config) :: this
    integer(4) :: id(:)
    if(any(this%is_any(id(:)))) stop "atom is already CGed"
    this%is_any(id(:)) = .true.
end subroutine
!----------------------------------------------------
integer(4) function get_MolSize(this) result(res)
    class(config) :: this
    integer(4) :: i
    res=1   
    do i = this%iter_id+1, this%n
        if(this%mol(i) /= this%mol(this%iter_id)) exit
        res = res + 1
    enddo
end function
!----------------------------------------------------
integer(4) function iter_NextMol(this) result(res)
    class(config) :: this
    integer(4) :: i
    if(this%iter_id == 0) then
        this%iter_id = 1
        this%iter_mol = this%mol(1)
        res = 1
        return
    endif
    do i = this%iter_id+1, this%n
        if(this%mol(i) /= this%mol(this%iter_id)) then
            this%iter_id = i
            this%iter_mol = this%mol(i)
            res = i
            exit
        endif
    enddo
    if(i == this%n+1) then
        call this%iter_reset()
        res = 0
    endif
end function
!----------------------------------------------------
pure logical function iter_hasNextMol(this)
    class(config),intent(in) :: this
    integer(4) :: i
    iter_hasNextMol = .false.
    do i = this%iter_id+1, this%n
        if(this%mol(i) /= this%mol(this%iter_id)) then
            iter_hasNextMol = .true.
            exit
        endif
    enddo
end function
!----------------------------------------------------
!subroutine config_set_values(one_line,id)
!    integer(4), intent(in) :: id
!    character(len=*),intent(in) :: one_line
!    integer(4) :: i_tmp
!    real(8) :: r_tmp
!    read(one_line,*) i_tmp,mol(id),tp(id),r_tmp,r(:,id)
!end subroutine
!----------------------------------------------------
subroutine set_box_size(this,d,lmn,lmx)
    class(config) :: this
    character, intent(in) :: d
    real(8),intent(in) :: lmn, lmx
    integer(4) :: i
    character :: opt(3) = ["x","y","z"]
    do i = 1, 3
        if(opt(i) == d) then
            this%Lmin(i) = lmn
            this%Lmax(i) = lmx
            this%Lbox(i) = lmx - lmn
            !write(*,'(2a,3(g0,1x))') opt(i)," set to",this%Lmin(i),this%Lmax(i),this%Lbox(i)
            exit
        endif
    enddo
end subroutine
end module
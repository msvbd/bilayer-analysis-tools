module lmp_traject_reader_mod
use config_mod
implicit none

private

integer(4) :: io
integer(4) :: at_num, item_num
integer(4),public,protected :: pre_step = 0
integer(4),public,protected :: cur_step = 0
integer(4),public,protected :: step_num = 0

integer(4) :: id_pos = 0
integer(4) :: mol_pos = 0
integer(4) :: rad_pos = 0
integer(4) :: type_pos = 0
integer(4) :: x_pos = 0
integer(4) :: y_pos = 0
integer(4) :: z_pos = 0

logical,public,protected :: is_end = .false.

public :: get_snap, open_trajectory, close_trajectory, reopen_trajectory

contains
!----------------------------------------------------------
type(config) function get_snap(typeFilter) result(conf)
    integer(4),optional :: typeFilter(:)
    step_num = step_num + 1
    
    call read_head(io, conf)
    if(present(typeFilter)) then
        call read_positions_filter(io, conf, typeFilter)
    else
        call read_positions(io, conf)
    endif
end function
!----------------------------------------------------------
subroutine open_trajectory(file_name)
    character(len=*),intent(in) :: file_name
    type(config) :: conf
    is_end = .false.
    pre_step = 0
    cur_step = 0
    step_num = 0
    open(newunit=io,file=file_name,action='read',status='old')
    call read_head(io, conf)

    write(*,*) "traject format: ",pack(["id   ","mol  ","type ","rad  ",&
                &"x    ","y    ","z    ","xu   ","yu   ","zu   "],&
                    [id_pos, mol_pos, type_pos, rad_pos,&
                    & x_pos, y_pos, z_pos, x_pos, y_pos, z_pos] /= 0)

    close(io)
    open(newunit=io,file=file_name,action='read',status='old')
end subroutine
!----------------------------------------------------------
subroutine reopen_trajectory()
    character(len=100) :: file_name
    inquire(unit=io,name=file_name)
    write(*,*) "Reopen: ",trim(file_name)
    call close_trajectory()
    call open_trajectory(trim(file_name))
end subroutine
!----------------------------------------------------------
subroutine close_trajectory()
    is_end = .true.
    close(io)
end subroutine
!----------------------------------------------------------
subroutine read_head(io, conf)
    integer(4),intent(in) :: io
    type(config),intent(inout) :: conf
    integer(4) :: istat, i, i_item
    real(8) :: Lb(3),Ub(3)
    character(len=100) :: line
    character(len=7) :: item
    
    do
        read(io,'(a)',iostat=istat) line
        !write(*,*) ">>>",trim(line)
        if(istat /= 0) then
            is_end = .true.
            return
        endif

        if(index(line,"ITEM: TIMESTEP") /= 0) then
            pre_step = cur_step
            read(io,*) cur_step
            !write(*,*) "ITEM: TIMESTEP = ",time_step
        elseif(index(line,"ITEM: NUMBER OF ATOMS") /= 0) then
            read(io,*) at_num
            if(at_num /= conf%n) call conf%set_n(at_num)
            !write(*,*) "ITEM: NUMBER OF ATOMS = ",at_num
        elseif(index(line,"ITEM: BOX BOUNDS") /= 0) then
            !write(*,*) "ITEM: BOX BOUNDS"
            do i=1,3
                read(io,*) lb(i), ub(i)
                !write(*,*) lb(i), ub(i)
            enddo
            call conf%set_box_size("x",lb(1),ub(1))
            call conf%set_box_size("y",lb(2),ub(2))
            call conf%set_box_size("z",lb(3),ub(3))
        elseif(index(line,"ITEM: ATOMS") /= 0) then
            !write(*,*) "ITEM: ATOMS"
            i = 13
            i_item = 0
            do
                read(line(i:),*) item
                i_item = i_item + 1
                if(trim(item) == "id") id_pos = i_item
                if(trim(item) == "mol") mol_pos = i_item
                if(trim(item) == "radius") rad_pos = i_item
                if(trim(item) == "type") type_pos = i_item
                if(trim(item) == "x") x_pos = i_item
                if(trim(item) == "y") y_pos = i_item
                if(trim(item) == "z") z_pos = i_item
                if(trim(item) == "xu") x_pos = i_item
                if(trim(item) == "yu") y_pos = i_item
                if(trim(item) == "zu") z_pos = i_item
                i = i + len_trim(item) + 1
                if(i > len_trim(line)) exit
            enddo
            item_num = max(id_pos, mol_pos, rad_pos, type_pos, &
                            & x_pos, y_pos, z_pos)
            if(type_pos == 0) type_pos = id_pos
            exit
        else
            write(*,*) "unknow ITEM -> ",trim(line)
            is_end = .true.
            exit
        endif
    enddo

end subroutine
!----------------------------------------------------------
subroutine read_positions(io, conf)
    integer(4),intent(in) :: io
    type(config),intent(inout) :: conf
    integer(4) :: i, istat
    real(8) :: items(0:20)

    items(0) = 0.d0
    !write(*,*) at_num
    do i = 1, at_num
        !if(mod(i,1000)==0) write(*,*) i,at_num
        !write(*,*) i,at_num
        read(io,*,iostat=istat) items(1:item_num)
        if(istat /= 0) then
            is_end = .true.
            exit
        endif

        call conf%add_particle_unsafe(&
            &[items(x_pos), items(y_pos), items(z_pos)],&
            nint(items(type_pos)),&
            nint(items(id_pos)))
        !if(rad_pos /= 0) &
        !    call conf%add_radius_id_unsafe(items(rad_pos), nint(items(id_pos)))
    enddo
end subroutine
!----------------------------------------------------------
subroutine read_positions_filter(io, conf, tf)
    integer(4),optional :: tf(:) ! Type filter
    integer(4),intent(in) :: io
    type(config),intent(inout) :: conf
    type(config) :: conf_tmp
    integer(4) :: i, istat, real_N
    real(8) :: items(0:20)

    call conf_tmp%set_n(conf%n)
    call read_positions(io, conf_tmp)

    items(0) = 0.d0
    real_N = 0
    do i = 1, conf_tmp%n
        if(any(conf_tmp%tp(i)==tf)) then    
            real_N = real_N + 1
            call conf%add_particle_unsafe(&
                 conf_tmp%r(:,i),&
                 conf_tmp%tp(i),&
                 real_N)
            !if(rad_pos /= 0) &
            !    call conf%add_radius_id_unsafe(conf_tmp%rad(i), real_N)
        endif
    enddo
    conf%n = real_N
end subroutine

end module

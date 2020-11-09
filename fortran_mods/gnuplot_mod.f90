module gnuplot_mod
    implicit none
    
    private
    
    type gnuplot
        private
        character(len=100),allocatable :: fce(:)
        character(len=100),allocatable :: data(:)
        character(len=100),allocatable :: sets(:)
        character(len=100),allocatable :: lines(:)
        character(len=200) :: terminal
        character(len=200) :: dir = "."
    
      contains
        procedure,public :: init

        procedure,public :: add
        procedure,public :: add_fce
        procedure,public :: add_line

        procedure,public :: rm_data
        procedure,public :: rm_sets
        procedure,public :: rm_lines

        procedure,public :: plot

        procedure,public :: set
        procedure,public :: set_terminal
        procedure,public :: set_dir
    end type
    
    public :: gnuplot
    
    contains
    !----------------------------------------------------
    subroutine init(this)
        class(gnuplot) :: this
        this%terminal = ""

        call this%set_terminal("pdfcairo",&
            'enhanced size 10cm,7cm font "Helvetica, 12" linewidth 2')
        call this%set("key","width 1. height 1. box")
        call this%set("style data","lines")

        call this%set("encoding","iso_8859_1")
        !call this%set("border","lw 1")

        call this%set("linetype 1","lc 'black'")
        call this%set("linetype 2","lc 'red'")
        call this%set("linetype 3","lc 'green'")
        call this%set("linetype 4","lc 'blue'")
        call this%set("linetype 5","lc 'purple'")
    end subroutine
    !----------------------------------------------------
    subroutine set_dir(this, dir)
        class(gnuplot) :: this
        character(len=*) :: dir
        this%dir = trim(dir)
    end subroutine
    !----------------------------------------------------
    subroutine rm_data(this)
        class(gnuplot) :: this
        if(allocated(this%data)) deallocate(this%data)
        if(allocated(this%fce)) deallocate(this%fce)
    end subroutine
    !----------------------------------------------------
    subroutine rm_lines(this)
        class(gnuplot) :: this
        if(allocated(this%data)) deallocate(this%lines)
    end subroutine
    !----------------------------------------------------
    subroutine rm_sets(this)
        class(gnuplot) :: this
        if(allocated(this%sets)) deallocate(this%sets)
    end subroutine
    !----------------------------------------------------
    subroutine add(this, in_file, u, t, w, a, s, note)
        class(gnuplot) :: this
        character(len=*),intent(in) :: in_file
        character(len=*),intent(in),optional :: u
        character(len=*),intent(in),optional :: t
        character(len=*),intent(in),optional :: w
        character(len=*),intent(in),optional :: a
        character(len=*),intent(in),optional :: s
        character(len=*),intent(in),optional :: note
        character(len=:),allocatable :: data_tmp(:)
        integer(4) sd

        if(allocated(this%data)) then
            allocate(data_tmp, source=this%data)
            deallocate(this%data)
            sd = size(data_tmp)
            allocate(character(len=100) :: this%data(sd+1))
            this%data(1:sd) = data_tmp
            deallocate(data_tmp)
        else
            sd = 0
            allocate(character(len=100) :: this%data(1))
        endif

        write(this%data(sd+1),'(*(g0))') trim(in_file),"'"
        if(present(u)) &
            write(this%data(sd+1),'(*(g0))') trim(this%data(sd+1))," u ",trim(u)
        if(present(w)) &
            write(this%data(sd+1),'(*(g0))') trim(this%data(sd+1))," w ",trim(w)
        if(present(s)) &
            write(this%data(sd+1),'(*(g0))') trim(this%data(sd+1))," smooth ",trim(s)
        if(present(a)) &
            write(this%data(sd+1),'(*(g0))') trim(this%data(sd+1))," axes ",trim(a)
        if(present(note)) &
            write(this%data(sd+1),'(*(g0))') trim(this%data(sd+1))," ",trim(note)
        if(present(t)) &
            write(this%data(sd+1),'(*(g0))') trim(this%data(sd+1))," t '",trim(t),"'"

    end subroutine
    !----------------------------------------------------
    subroutine add_line(this, note)
        class(gnuplot) :: this
        character(len=*),intent(in) :: note
        character(len=:),allocatable :: lines_tmp(:)
        integer(4) sd

        if(allocated(this%lines)) then
            sd = size(this%lines)
            allocate(lines_tmp, source=this%lines)
            deallocate(this%lines)
            
            allocate(character(len=100) :: this%lines(sd+1))
            this%lines(1:sd) = lines_tmp
            deallocate(lines_tmp)
        else
            sd = 0
            allocate(character(len=100) :: this%lines(1))
        endif

        write(this%lines(sd+1),'(*(g0))') trim(note)

    end subroutine
    !----------------------------------------------------
    subroutine add_fce(this, fce, t, w, a, note)
        class(gnuplot) :: this
        character(len=*),intent(in) :: fce
        character(len=*),intent(in),optional :: t
        character(len=*),intent(in),optional :: w
        character(len=*),intent(in),optional :: a
        character(len=*),intent(in),optional :: note
        character(len=:),allocatable :: fce_tmp(:)
        integer(4) sd

        if(allocated(this%fce)) then
            allocate(fce_tmp, source=this%fce)
            deallocate(this%fce)
            sd = size(fce_tmp)
            allocate(character(len=100) :: this%fce(sd+1))
            this%fce(1:sd) = fce_tmp
            deallocate(fce_tmp)
        else
            sd = 0
            allocate(character(len=100) :: this%fce(1))
        endif

        write(this%fce(sd+1),'(*(g0))') trim(fce)
        if(present(w)) &
            write(this%fce(sd+1),'(*(g0))') trim(this%fce(sd+1))," w ",trim(w)
        if(present(a)) &
            write(this%fce(sd+1),'(*(g0))') trim(this%fce(sd+1))," axes ",trim(a)
        if(present(note)) &
            write(this%fce(sd+1),'(*(g0))') trim(this%fce(sd+1))," ",trim(note)
        if(present(t)) &
            write(this%fce(sd+1),'(*(g0))') trim(this%fce(sd+1))," t '",trim(t),"'"

    end subroutine
    !----------------------------------------------------
    subroutine set(this, what, how)
        class(gnuplot) :: this
        character(len=*),intent(in) :: what
        character(len=*),intent(in),optional :: how
        character(len=:),allocatable :: sets_tmp(:)
        integer(4) sd, i
        logical :: there_is

        there_is = .false.

        if(allocated(this%sets)) then
            sd = size(this%sets)
            do i=1,sd
                if(this%sets(i)(5:5+len_trim(what)) == trim(what) .or. &
                    this%sets(i)(7:7+len_trim(what)) == trim(what)) then
                    there_is = .true.
                    exit
                endif
            enddo
            allocate(sets_tmp, source=this%sets)
            deallocate(this%sets)
            if(there_is) then
                allocate(character(len=100) :: this%sets(sd))
                if(i == 1) then
                    this%sets(1:sd-1) = sets_tmp(2:sd)
                elseif(i == sd) then
                    this%sets = sets_tmp
                else
                    this%sets(1:i-1) = sets_tmp(1:i-1)
                    this%sets(i:sd-1) = sets_tmp(i+1:sd)
                endif
                sd = sd - 1
            else
                allocate(character(len=100) :: this%sets(sd+1))
                this%sets(1:sd) = sets_tmp
            endif
            deallocate(sets_tmp)
        else
            sd = 0
            allocate(character(len=100) :: this%sets(1))
        endif

        if(present(how)) then
            write(this%sets(sd+1),'(*(g0))') "set ",trim(what)," ",trim(how)
        else
            write(this%sets(sd+1),'(*(g0))') "unset ",trim(what)
        endif

    end subroutine
    !----------------------------------------------------
    subroutine plot(this, out_file,xr, yr , xl, yl, rm)
        class(gnuplot) :: this
        character(len=*),optional :: xr
        character(len=*),optional :: yr
        character(len=*),optional :: xl
        character(len=*),optional :: yl
        character(len=*) :: out_file
        character(len=100) :: gpl_file
        integer(4) io, i, dot, gi(5)
        real gr(5)
        logical,optional :: rm
        logical :: rmp = .true.
        
        call random_seed()
        if(present(rm)) rmp = rm

        if(.not. allocated(this%data)) return

        dot = index(out_file, ".", .true.)
        call random_number(gr)
        gi=int(gr*10.)
        if(rmp) then
            write(gpl_file,'(*(g0))') trim(this%dir),"/","gnuplot_mod_",&
                &trim(out_file(:dot-1)),".gpl"
        else
            write(gpl_file,'(*(g0))') trim(this%dir),"/","gnuplot_mod_",&
                &trim(out_file(:dot-1)),"_",gi,".gpl"
        endif

        open(newunit=io, file=gpl_file)
        write(io,'(g0)') this%terminal

        if(allocated(this%sets)) then
            do i=1, size(this%sets)
                write(io,'(*(g0))') trim(this%sets(i))
            enddo
        endif
        if(allocated(this%lines)) then
            do i=1, size(this%lines)
                write(io,'(*(g0))') trim(this%lines(i))
            enddo
        endif

        if(present(xl)) write(io,*) "set xlabel '"//trim(xl)//"'"
        if(present(yl)) write(io,*) "set ylabel '"//trim(yl)//"'"
        if(present(xr)) write(io,*) "set xrange ["//trim(xr)//"]"
        if(present(yr)) write(io,*) "set yrange ["//trim(yr)//"]"

        write(io,'(*(g0))') "set output '",trim(this%dir),"/",trim(out_file),"'"

        write(io,'(g0)',advance='no') "plot "
        if(allocated(this%data)) then
            if(size(this%data)>0) &
                & write(io,'(*(g0))',advance='no') "'",trim(this%dir),"/",trim(this%data(1))
            do i=2, size(this%data)
                write(io,'(*(g0))',advance='no') ",'",trim(this%dir),"/",trim(this%data(i))
            enddo
        endif

        if(allocated(this%fce)) then
            do i=1, size(this%fce)
                write(io,'(*(g0))',advance='no') ",",trim(this%fce(i))
            enddo
        endif
        write(io,*) ""
        !write(io,'(*(g0))') trim(this%fce(i))

        write(io,'(g0)') "unset output"
        close(io)

        call execute_command_line('sed -i -E "s/@dir/'//trim(this%dir)//'/g" '//gpl_file)
        if(rmp) then
            call execute_command_line('gnuplot '//gpl_file//' && rm '//gpl_file)
        else
           call execute_command_line('gnuplot '//gpl_file)
        endif
        !call execute_command_line('gnuplot '//gpl_file)
        !call execute_command_line("rm gnuplot_mod_plot.gpl")

    end subroutine
    !----------------------------------------------------
    subroutine set_terminal(this, tp, note)
        class(gnuplot) :: this
        character(len=*) :: tp
        character(len=*),optional :: note


        write(this%terminal,'(*(g0))') "set terminal ",tp
        if(present(note)) write(this%terminal,'(*(g0))') trim(this%terminal)//" "//note
    end subroutine
end module
!========================================================
!program test
!use gnuplot_mod
!implicit none

!integer(4) :: io,i,n
!real(8) :: dx
!type(gnuplot) :: gplo

!n=20
!dx = 2.*3.14/dble(n)

!open(newunit=io, file='sinus.data')
!    do i=0,n
!        write(io,*) i*dx, sin(i*dx)
!    enddo
!close(io)

!open(newunit=io, file='cos.data')
!    do i=0,n
!        write(io,*) i*dx, cos(i*dx)
!    enddo
!close(io)

!call gplo%init()

!call gplo%set("title","'goniometricke fce'")
!call gplo%set("xlabel","'x'")
!call gplo%set("ylabel","'y'")

!call gplo%add('sinus.data',t="sinus",w="p")
!call gplo%add('cos.data',t="cosinus",w="p")
!call gplo%add_fce('sin(x)',t="sin(x)")
!call gplo%add_fce('cos(x)',t="cos(x)")

!call gplo%plot("goniom.png","[0:6]")

!end program
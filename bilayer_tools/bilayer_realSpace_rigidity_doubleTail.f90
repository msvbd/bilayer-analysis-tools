!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! The program computes rigidity of bilayer which contains a polymer with two tails.
!
! inputs:
!         - name of output directory - the first argument in command line
!         - name of trajectory file  - the second argument in command line
!         - setting variables:
!             - tFilter - the array of particle types of bilayer
!
! outputs:
!         - neighbots.data            - the number of neighbours polymer of type "i" around the type "j"
!         - PSi.data                  - the distribution of "S_i"
!         - RS_rigidity_fit_gpl.data  - It contains the number of each type of polymers, the area 
!                                       of the bilayer (Lx*Ly) and the results of the fit
!                                       of "P(S_i)" (by Gnuplot).
!         - PSi.png                   - plot of "P(S_i)" (by Gnuplot)
!         - PSi_fit.png               - plot of "P(S_i)" with fitting functions (by Gnuplot)
!
! Created by: Martin Svoboda (svobodam@icpf.cas.cz )
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program realSpaceRigidityDT
use config_mod
use gnuplot_mod
use lmp_traject_reader_mod
use mol_pattern_mod
implicit none

! === array dimension variables =====================================
integer(4),parameter :: nb = 200   ! half the number of bins in distributions
real(8),parameter    :: rcut = 2.0 ! cutoff for neighbour counting

! === setting =======================================================
integer(4),parameter :: tFilter(7) = [3,4,5,6,7,8,9]  ! diCTAC - trajectory

! === initiation of other variables ================================= 
character(len=500) blabla
character(len=100) charus, file_name

integer(4) num_of_args, id1, id2
integer(4) tp1, tp2, tpmin, tpmax, i, io
integer(4) np1, np2

real(8) :: rij(3), h, n1(3), n2(3), si(3)
real(8) :: neighborts(2,2), area
real(8) :: PSi(2,2,-nb:nb), dpsi

logical existuje

type(config) :: cnfg
type(mol_patts) :: patts
type(gnuplot) :: gnupl

! === read command line input =======================================
num_of_args = command_argument_count()
if(num_of_args < 2) stop "num_of_args < 2"

call get_command_argument(1, charus)
call get_command_argument(2, file_name)

write(*,*) "charus: ",trim(charus)
write(*,*) "file name: ",trim(file_name)

! === find a pattern in the trajectory file =========================
call open_trajectory(trim(file_name))
cnfg = get_snap(tFilter)
call patts%find_patts(cnfg%tp, [3,4,5], [".","+","."] ) ! FOL 
call patts%find_patts(cnfg%tp, &
                        [6,7,8,9,7,8,9], &
                        [".",".","+",".",".","+","."] ) !diCTAC            
if(size(patts%patts) < 2) stop "size( patts ) < 2"

! === count the number of polymers of both types ====================
id1 = 1
np1 = 0
np2 = 0
area = cnfg%Lbox(1) * cnfg%Lbox(2)
do while(id1<=cnfg%n)
    if(cnfg%tp(id1) == patts%patts(1)%patt(1)) then
        np1 = np1 + 1.d0
        id1 = id1 + patts%patts(1)%n
    elseif(cnfg%tp(id1) == patts%patts(2)%patt(1)) then
        np2 = np2 + 1
        id1 = id1 + patts%patts(2)%n
    else
        id1 = id1 + 1
        stop ""
    endif
enddo

! === set values of some variables ==================================
PSi = 0.d0
neighborts = 0.d0
dpsi = 2.d0 / dble(nb)

! === the main loop through the trajectories ========================
call reopen_trajectory()
cnfg = get_snap(tFilter)
do while(.not. is_end)

    write(*,*) cur_step, step_num, cnfg%n

    ! === the outer loop of double loop =============================
    id1 = 1
    do while(id1<=cnfg%n-1)
        if(cnfg%tp(id1) == patts%patts(1)%patt(1)) then
            tp1 = 1

            n1 = get_direction(cnfg%r(:,id1), cnfg%r(:,id1+patts%patts(1)%n-1), cnfg%Lbox)

            id1 = id1 + patts%patts(1)%n
        elseif(cnfg%tp(id1) == patts%patts(2)%patt(1)) then
            tp1 = 2

            n1 = get_direction(cnfg%r(:,id1), &
                                (cnfg%r(:,id1+patts%patts(2)%n-1) + &
                                    cnfg%r(:,id1+patts%patts(2)%n-10))*0.5d0, &
                                cnfg%Lbox)

            id1 = id1 + patts%patts(2)%n
        else
            !stop "???"
            id1 = id1 + 1
            cycle
        endif

        ! === the inner loop of double loop =========================
        id2 = id1 + 1
        do while(id2<=cnfg%n)
        
            if(cnfg%tp(id2) == patts%patts(1)%patt(1)) then
                tp2 = 1    

                n2 = get_direction(cnfg%r(:,id2), cnfg%r(:,id2+patts%patts(1)%n-1), cnfg%Lbox)

                id2 = id2 + patts%patts(1)%n
            elseif(cnfg%tp(id2) == patts%patts(2)%patt(1)) then
                tp2 = 2

                n2 = get_direction(cnfg%r(:,id2), &
                                    (cnfg%r(:,id2+patts%patts(2)%n-1) + &
                                        cnfg%r(:,id2+patts%patts(2)%n-10))*0.5d0, &
                                    cnfg%Lbox)

                id2 = id2 + patts%patts(2)%n
            else
                !stop "???"
                id2 = id2 + 1
                cycle
            endif

            rij = cnfg%r(:,id1) - cnfg%r(:,id2)
            rij = rij - cnfg%Lbox * nint(rij/ cnfg%Lbox)
            h = dsqrt(dot_product(rij, rij))
            if(h > rcut) cycle
            si = (n2 - n1) / h

            tpmin = min(tp1, tp2)
            tpmax = max(tp1, tp2)
            neighborts(tpmin, tpmax) = neighborts(tpmin, tpmax) + 1.d0
            i = nint(si(1) / dpsi)
            if (abs(i) <= nb) PSi(tpmin, tpmax, i)        = PSi(tpmin, tpmax, i)  + 1.d0
            i = nint(si(2) / dpsi)
            if (abs(i) <= nb) PSi(tpmin, tpmax, i)        = PSi(tpmin, tpmax, i)  + 1.d0

        enddo
    enddo

    !if(step_num > 10) exit
    cnfg = get_snap(tFilter)
enddo
write(*,'(2a)') "log: <<< end >>>"
call close_trajectory()

! === the normalisation of the results ==============================
neighborts = neighborts / dble(step_num)
PSi = PSi / dble(step_num)

! === make the result directory if doesn't exist ===================
inquire(file=trim(charus),exist=existuje)
if(.not. existuje) then
    call execute_command_line('mkdir '//trim(charus))
endif  

! === print the output files ========================================
open(io,file=trim(charus)//"/neighbots.data")
associate(ngb => neighborts)
write(io,*) "(1,1) (1,2) (2,1) (2,2)"
write(io,*) ngb(1,1), ngb(1,2), ngb(2,1), ngb(2,2)
end associate
close(io)

open(io, file=trim(charus)//"/PSi.data")
write(io,*) "S, PSi(1,1) PSi(1,2) PSi(2,1) PSi(2,2)"
do i = -nb, nb
    write(io,*) i*dpsi, PSi(1,1,i)/sum(PSi(1,1,:)), &
                        PSi(1,2,i)/sum(PSi(1,2,:)), &
                        PSi(2,1,i)/sum(PSi(2,1,:)), &
                        PSi(2,2,i)/sum(PSi(2,2,:))
enddo
close(io)

! === plot and fit the output files by Gnuplot ======================
call gnupl%init()

call gnupl%set_terminal("pngcairo",&
    'enhanced size 800,600 font "Helvetica, 12" linewidth 2')

call gnupl%set_dir(trim(charus))
call gnupl%set("key","autotitle columnhead ")


call gnupl%add_line("set print '@dir/RS_rigidity_fit_gpl.data'")
call gnupl%add_line("k11=10.0; q11=1.0; f11(x)=-k11*x**2 + q11")
call gnupl%add_line("k12=10.0; q12=1.0; f12(x)=-k12*x**2 + q12")
call gnupl%add_line("k21=10.0; q21=1.0; f21(x)=-k21*x**2 + q21")
call gnupl%add_line("k22=10.0; q22=1.0; f22(x)=-k22*x**2 + q22")
write(charus,'(*(g0))') "fit [-0.05:0.05]&
                            & f11(x) '@dir/PSi.data' u 1:(2*log($2)) via k11, q11"
call gnupl%add_line(charus)
write(charus,'(*(g0))') "fit [-0.05:0.05]&
                            & f12(x) '@dir/PSi.data' u 1:(2*log($3)) via k12, q12"
call gnupl%add_line(charus)
write(charus,'(*(g0))') "fit [-0.05:0.05]&
                            & f22(x) '@dir/PSi.data' u 1:(2*log($5)) via k22, q22"
call gnupl%add_line(charus)
write(charus,'(*(g0))') "print 'n1 ",np1,"'"
call gnupl%add_line(charus)
write(charus,'(*(g0))') "print 'n2 ",np2,"'"
call gnupl%add_line(charus)
write(charus,'(*(g0))') "print 'Area ",area,"'"
call gnupl%add_line(charus)
call gnupl%add_line("print '(1,1) ',k11,' ',q11")
call gnupl%add_line("print '(1,2) ',k12,' ',q12")
call gnupl%add_line("print '(2,2) ',k22,' ',q22")
call gnupl%rm_data()
    call gnupl%add('PSi.data',u="1:(2*log($2))", t="(1,1)")
    call gnupl%add('PSi.data',u="1:(2*log($3))", t="(1,2)")
    !call gnupl%add('PSi.data',u="1:4", t="(2,1)")
    call gnupl%add('PSi.data',u="1:(2*log($5))", t="(2,2)")
    call gnupl%add_fce('f11(x)',w="l")
    call gnupl%add_fce('f12(x)',w="l")
    !call gnupl%add_fce('f21(x)',w="l")
    call gnupl%add_fce('f22(x)',w="l")
        call gnupl%plot("PSi_fit.png",xl="S", yl="P(S)", xr="-0.3:0.3")
call gnupl%rm_lines()

call gnupl%rm_data()
call gnupl%add('PSi.data',u="1:2", t="(1,1)")
call gnupl%add('PSi.data',u="1:3", t="(1,2)")
call gnupl%add('PSi.data',u="1:4", t="(2,1)")
call gnupl%add('PSi.data',u="1:5", t="(2,2)")
    call gnupl%plot("PSi.png",xl="S", yl="P(S)")
    
!====================================================================
contains
!--------------------------------------------------------------------
!
! The function gets the unit vector directed from end to end of the polymer.
!
function get_direction(v1, v2, L)
    real(8),intent(in) :: v1(3), v2(3), L(3)
    real(8) :: get_direction(3)
    get_direction = v1 - v2
    get_direction = get_direction - L * nint(get_direction / L)
    get_direction = get_direction/dsqrt(dot_product(get_direction,get_direction))
    !write(*,*) get_direction
end function
end program
    

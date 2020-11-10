!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! The program computes overlap parameter distribution of polymers in monolayers.
!
! inputs:
!         - name of output directory - the first argument in command line
!         - name of trajectory file  - the second argument in command line
!         - setting variables:
!             - tFilter - the array of particle types of bilayer
!             - ddove, ddc, dlz, dd
!
! outputs:
!         - Doverlap_dist.data - distribution of overlap 
!         - Dc_dist.data       - distribution of Dc
!         - Lz_dist.data       - distribution of Lz
!         - Doverlap_dist.png  - plot of distribution of overlap (by Gnuplot)
!         - Dc_dist.png        - plot of distribution of Dc (by Gnuplot)
!         - Lz_dist.png        - plot of distribution of Lz (by Gnuplot)
!
! refs:
!         - bilayer analysis:
!               Venturoli, M., Sperotto, M. M., Kranenburg, M., & Smit, B. (2006).
!               Mesoscopic models of biological membranes. Physics Reports, 437(1-2), 1-54.
!
! Created by: Martin Svoboda (svobodam@icpf.cas.cz )
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program DcLzDoverlap
use ieee_arithmetic
use config_mod
use gnuplot_mod
use lmp_traject_reader_mod
use mol_pattern_mod
implicit none

! === array dimension variables =====================================
integer(4),parameter :: nb=200    ! half the number of bins in distributions
integer(4),parameter :: nbwMax=20 ! size of 2d lattice of results

! === setting =======================================================
integer(4),parameter :: tFilter(4) = [3,4,5,8]         ! CTAC   - trajectory
!integer(4),parameter :: tFilter(6) = [3,4,5,7,8,9] ! diCTAC - trajectory

real(4),parameter :: ddove = 0.05  ! size of bins in overlap distribution
real(4),parameter :: ddc   = 0.05  ! size of bins in Dc distribution
real(4),parameter :: dlz   = 0.025 ! size of bins in Lz distribution
real(4),parameter :: dd    = 2.5   ! the length between neighbours lattice nodes

! === initiation of other variables ================================= 
integer(4) :: n_fol
integer(4) :: n_xtac

integer(4) i,j,k,ix,iy,kx,ky,nbw(2),dir, io

real(4) Lz_FOL(-nbwMax:nbwMax,-nbwMax:nbwMax),Dc_FOL(-nbwMax:nbwMax,-nbwMax:nbwMax),&
&z1t1_FOL(-nbwMax:nbwMax,-nbwMax:nbwMax)
real(4) z2t1_FOL(-nbwMax:nbwMax,-nbwMax:nbwMax),zt1_FOL(2,-nbwMax:nbwMax,-nbwMax:nbwMax),&
&ztn_FOL(2,-nbwMax:nbwMax,-nbwMax:nbwMax)
real(4) Lz_XTAC(-nbwMax:nbwMax,-nbwMax:nbwMax),Dc_XTAC(-nbwMax:nbwMax,-nbwMax:nbwMax),&
&z1t1_XTAC(-nbwMax:nbwMax,-nbwMax:nbwMax)
real(4) z2t1_XTAC(-nbwMax:nbwMax,-nbwMax:nbwMax),zt1_XTAC(2,-nbwMax:nbwMax,-nbwMax:nbwMax),&
&ztn_XTAC(2,-nbwMax:nbwMax,-nbwMax:nbwMax)
real(4) Nz1t1_FOL(-nbwMax:nbwMax,-nbwMax:nbwMax),Nz2t1_FOL(-nbwMax:nbwMax,-nbwMax:nbwMax),&
&Nzt1_FOL(2,-nbwMax:nbwMax,-nbwMax:nbwMax)
real(4) Nztn_FOL(2,-nbwMax:nbwMax,-nbwMax:nbwMax)
real(4) Nz1t1_XTAC(-nbwMax:nbwMax,-nbwMax:nbwMax),Nz2t1_XTAC(-nbwMax:nbwMax,-nbwMax:nbwMax),&
&Nzt1_XTAC(2,-nbwMax:nbwMax,-nbwMax:nbwMax)
real(4) Nztn_XTAC(2,-nbwMax:nbwMax,-nbwMax:nbwMax)

real(8) Dove_FOL_dist(-nb:nb),Dove_XTAC_dist(-nb:nb)
real(8) Lz_FOL_dist(-nb:nb),Lz_XTAC_dist(-nb:nb)
real(8) Dc_FOL_dist(-nb:nb),Dc_XTAC_dist(-nb:nb)

character(len=500) blabla
character(len=100) charus, file_name

integer(4) num_of_args

logical existuje

type(mol_patts) :: patts
type(config) :: cnfg
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
call patts%find_patts(cnfg%tp, [8,4,5], [".","+","."] ) ! CTAC
!call patts%find_patts(cnfg%tp, [7,8,9], [".","+","."] ) !diCTAC tails
if(size(patts%patts) < 2) stop "size( patts ) < 2"

! === set values of some variables =======================================
n_xtac = patts%patts(2)%n
n_fol = patts%patts(1)%n

nbw(1) = int(0.5*cnfg%Lbox(1)/dd)
nbw(2) = int(0.5*cnfg%Lbox(2)/dd)
write(*,'(a,3f15.10)') "Lb ::: ",cnfg%Lbox(:)
write(*,'(a,2i15)') "nbw ::: ",nbw(:)
if(nbw(1)>nbwMax) stop "nbw(1)>nbwMax"
if(nbw(2)>nbwMax) stop "nbw(2)>nbwMax"

Dc_FOL_dist=0.0
Dc_XTAC_dist=0.0

Lz_FOL_dist=0.0
Lz_XTAC_dist=0.0

Dove_FOL_dist=0.0
Dove_XTAC_dist=0.0

! === the main loop through the trajectories ========================
call reopen_trajectory()
cnfg = get_snap(tFilter)
do while(.not. is_end)

    write(*,*) cur_step, step_num, cnfg%n

    ! === zeros accumulated variables ===============================
    z1t1_FOL = 0.0
    Nz1t1_FOL = 0.0
    z2t1_FOL = 0.0
    Nz2t1_FOL = 0.0
    zt1_FOL = 0.0
    Nzt1_FOL = 0.0
    ztn_FOL = 0.0
    Nztn_FOL = 0.0

    z1t1_XTAC = 0.0
    Nz1t1_XTAC = 0.0
    z2t1_XTAC = 0.0
    Nz2t1_XTAC = 0.0
    zt1_XTAC = 0.0
    Nzt1_XTAC = 0.0
    ztn_XTAC = 0.0
    Nztn_XTAC = 0.0
    
    ! === the loop through a one configuration ======================
    i = 1
    do while(i<=cnfg%n)
      ! === if the next sequence starts the same as the pattern(1)
      if(cnfg%tp(i) == patts%patts(1)%patt(1)) then

        j=n_fol-1
          
        dir = get_direction(cnfg%r(3,i+j), cnfg%r(3,i))

        kx=nint(cnfg%r(1,i)/dd)
        ky=nint(cnfg%r(2,i)/dd)
        
        zt1_FOL(dir,kx,ky)=zt1_FOL(dir,kx,ky)+cnfg%r(3,i+j)
        Nzt1_FOL(dir,kx,ky)=Nzt1_FOL(dir,kx,ky)+1.0
        ztn_FOL(dir,kx,ky)=ztn_FOL(dir,kx,ky)+cnfg%r(3,i+1)
        Nztn_FOL(dir,kx,ky)=Nztn_FOL(dir,kx,ky)+1.0
                        
        if(dir==1) then
          z1t1_FOL(kx,ky)=z1t1_FOL(kx,ky)+cnfg%r(3,i+1)
          Nz1t1_FOL(kx,ky)=Nz1t1_FOL(kx,ky)+1.0
        else
          z2t1_FOL(kx,ky)=z2t1_FOL(kx,ky)+cnfg%r(3,i+1)
          Nz2t1_FOL(kx,ky)=Nz2t1_FOL(kx,ky)+1.0
        endif
        
        i=i+n_fol

        i = i + patts%patts(1)%n

      ! === if the next sequence starts the same as the pattern(2)
      elseif(cnfg%tp(i) == patts%patts(2)%patt(1)) then

        j=n_xtac-1    
        
        dir = get_direction(cnfg%r(3,i+j), cnfg%r(3,i))
      
        kx=nint(cnfg%r(1,i)/dd)
        ky=nint(cnfg%r(2,i)/dd)
        
        zt1_XTAC(dir,kx,ky)=zt1_XTAC(dir,kx,ky)+cnfg%r(3,i+j)
        Nzt1_XTAC(dir,kx,ky)=Nzt1_XTAC(dir,kx,ky)+1.0
        ztn_XTAC(dir,kx,ky)=ztn_XTAC(dir,kx,ky)+cnfg%r(3,i+1)
        Nztn_XTAC(dir,kx,ky)=Nztn_XTAC(dir,kx,ky)+1.0
                      
        if(dir==1) then
          z1t1_XTAC(kx,ky)=z1t1_XTAC(kx,ky)+cnfg%r(3,i+1)
          Nz1t1_XTAC(kx,ky)=Nz1t1_XTAC(kx,ky)+1.0
        else
          z2t1_XTAC(kx,ky)=z2t1_XTAC(kx,ky)+cnfg%r(3,i+1)
          Nz2t1_XTAC(kx,ky)=Nz2t1_XTAC(kx,ky)+1.0
        endif
        
        i=i+n_xtac

        i = i + patts%patts(2)%n
      else
        i = i + 1
        cycle
      endif
    enddo
        
    ! === compute average values on lattice =========================
    z1t1_FOL = z1t1_FOL/Nz1t1_FOL
    z2t1_FOL = z2t1_FOL/Nz2t1_FOL
    zt1_FOL = zt1_FOL/Nzt1_FOL
    ztn_FOL = ztn_FOL/Nztn_FOL

    z1t1_XTAC = z1t1_XTAC/Nz1t1_XTAC
    z2t1_XTAC = z2t1_XTAC/Nz2t1_XTAC
    zt1_XTAC = zt1_XTAC/Nzt1_XTAC
    ztn_XTAC = ztn_XTAC/Nztn_XTAC

    Dc_FOL = abs(z1t1_FOL - z2t1_FOL)
    Dc_XTAC = abs(z1t1_XTAC - z2t1_XTAC)

    Lz_FOL = 0.5*(abs(zt1_FOL(1,:,:) - ztn_FOL(1,:,:))+abs(zt1_FOL(2,:,:) - ztn_FOL(2,:,:)))
    Lz_XTAC = 0.5*(abs(zt1_XTAC(1,:,:) - ztn_XTAC(1,:,:))+abs(zt1_XTAC(2,:,:) - ztn_XTAC(2,:,:)))
    
    ! === compute distributions of values on the lattice ============
    do ix=-nbw(1),nbw(1)
      do iy=-nbw(2),nbw(2)

        if((.not. ieee_is_nan(Lz_FOL(ix,iy))) .and. (.not. ieee_is_nan(Dc_FOL(ix,iy)))) then
          k = nint(((2.d0*Lz_FOL(ix,iy)-Dc_FOL(ix,iy))/Lz_FOL(ix,iy))/ddove)
          if(k>=-nb .and. k<=nb) Dove_FOL_dist(k) = Dove_FOL_dist(k) + 1.d0
          k = nint(Lz_FOL(ix,iy)/dlz)
          if(k>=-nb .and. k<=nb) Lz_FOL_dist(k) = Lz_FOL_dist(k) + 1.d0
          k = nint(Dc_FOL(ix,iy)/ddc)
          if(k>=-nb .and. k<=nb) Dc_FOL_dist(k) = Dc_FOL_dist(k) + 1.d0
        endif
        
        if((.not. ieee_is_nan(Lz_XTAC(ix,iy))) .and. (.not. ieee_is_nan(Dc_XTAC(ix,iy)))) then
          k = nint(((2.d0*Lz_XTAC(ix,iy)-Dc_XTAC(ix,iy))/Lz_XTAC(ix,iy))/ddove)
          if(k>=-nb .and. k<=nb) Dove_XTAC_dist(k) = Dove_XTAC_dist(k) + 1.d0
          k = nint(Lz_XTAC(ix,iy)/dlz)
          if(k>=-nb .and. k<=nb) Lz_XTAC_dist(k) = Lz_XTAC_dist(k) + 1.d0
          k = nint(Dc_XTAC(ix,iy)/ddc)
          if(k>=-nb .and. k<=nb) Dc_XTAC_dist(k) = Dc_XTAC_dist(k) + 1.d0
        endif
        
      enddo
    enddo

    !if(step_num > 10) exit
    cnfg = get_snap(tFilter)
enddo

write(*,'(a)') "log: <<< end >>>"
call close_trajectory()

! === the normalisation of the results ==============================
Dove_FOL_dist = Dove_FOL_dist / (sum(Dove_FOL_dist)*ddove)
Dove_XTAC_dist = Dove_XTAC_dist / (sum(Dove_XTAC_dist)*ddove)

Dc_FOL_dist = Dc_FOL_dist / (sum(Dc_FOL_dist)*ddc)
Dc_XTAC_dist = Dc_XTAC_dist / (sum(Dc_XTAC_dist)*ddc)

Lz_FOL_dist = Lz_FOL_dist / (sum(Lz_FOL_dist)*dlz)
Lz_XTAC_dist = Lz_XTAC_dist / (sum(Lz_XTAC_dist)*dlz)

! === make the result directory if doesn't exist ===================
inquire(file=trim(charus),exist=existuje)
  if(.not. existuje) then
    !stop "neexistuje"
    call execute_command_line("mkdir "//trim(charus))
  endif
  
! === print the output files ========================================
open(1,file=trim(charus)//"/Doverlap_dist.data")
write(1,*) "overlap FOL BTAC"
do i=-nb,nb
  write(1,*) real(i)*ddove,Dove_FOL_dist(i),Dove_XTAC_dist(i)
enddo
close(1)

open(1,file=trim(charus)//"/Dc_dist.data")
write(1,*) "Dc FOL BTAC"
do i=-nb,nb
  write(1,*) real(i)*ddc,Dc_FOL_dist(i),Dc_XTAC_dist(i)
enddo
close(1)

open(1,file=trim(charus)//"/Lz_dist.data")
write(1,*) "Lz FOL BTAC"
do i=-nb,nb
  write(1,*) real(i)*dlz,Lz_FOL_dist(i),Lz_XTAC_dist(i)
enddo
close(1)

! === plot the output files by Gnuplot ==============================
call gnupl%init()

call gnupl%set_terminal("pngcairo",&
'enhanced size 800,600 font "Helvetica, 12" linewidth 2')

call gnupl%set_dir(trim(charus))
call gnupl%set("key","autotitle columnhead")

call gnupl%set("ylabel","'dist'")

call gnupl%rm_data()
  call gnupl%add('Doverlap_dist.data',u="1:2")
  call gnupl%add('Doverlap_dist.data',u="1:3")
    call gnupl%plot("Doverlap_dist.png", xl="Doverlap")

call gnupl%rm_data()
  call gnupl%add('Dc_dist.data',u="1:2")
  call gnupl%add('Dc_dist.data',u="1:3")
    call gnupl%plot("Dc_dist.png", xl="Dc")

call gnupl%rm_data()
  call gnupl%add('Lz_dist.data',u="1:2")
  call gnupl%add('Lz_dist.data',u="1:3")
    call gnupl%plot("Lz_dist.png", xl="Lz")


!====================================================================
contains
!--------------------------------------------------------------------
!
! The function gets the direction of the polymer.
! Gets 1 or 2 for up or down. 
!
pure real(8) function get_direction(z1, z2) result(dir)
    real(8),intent(in) :: z1, z2
    if(z2-z1 > 0.d0) then
      dir = 1
    else
      dir = 2
    endif
end function
end program

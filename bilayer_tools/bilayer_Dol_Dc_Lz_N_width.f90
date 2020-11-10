!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! The program computes overlap parameter of polymers in monolayers
! as function of width of bilayer.
!
! inputs:
!         - name of output directory - the first argument in command line
!         - name of trajectory file  - the second argument in command line
!         - setting variables:
!             - tFilter - the array of particle types of bilayer
!             - dwdt
!
! outputs:
!         - Doverlap_width_dist.data - distribution of overlap 
!         - Dc_width_dist.data       - distribution of Dc
!         - Lz_width_dist.data       - distribution of Lz
!         - N_width_dist.data        - the number of cells with a given width
!         - Doverlap_width_dist.png  - plot of distribution of overlap (by Gnuplot)
!         - Dc_width_dist.png        - plot of distribution of Dc (by Gnuplot)
!         - Lz_width_dist.png        - plot of distribution of Lz (by Gnuplot)
!         - N_width_dist.png         - plot of the number of cells with a given width (by Gnuplot)
!
! refs:
!         - ITIM method:
!               Sega, M., Kantorovich, S. S., Jedlovszky, P., & Jorge, M. (2013).
!               The generalized identification of truly interfacial molecules (ITIM) algorithm
!               for nonplanar interfaces. The Journal of chemical physics, 138(4), 044110.

!         - bilayer analysis:
!               Venturoli, M., Sperotto, M. M., Kranenburg, M., & Smit, B. (2006).
!               Mesoscopic models of biological membranes. Physics Reports, 437(1-2), 1-54.
!
! Created by: Martin Svoboda (svobodam@icpf.cas.cz )
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program DolDcLzNWidth
use ieee_arithmetic
use config_mod
use gnuplot_mod
use lmp_traject_reader_mod
use mol_pattern_mod
implicit none

! === array dimension variables =====================================
integer(4),parameter :: nn=180000 ! max number of particles in one configuration 
integer(4),parameter :: nnp=50000 ! max number of particles after filtering
integer(4),parameter :: nb=200    ! half the number of bins in distributions
integer(4),parameter :: nbwMax=20 ! size of 2d lattice of results

! === setting =======================================================
integer(4),parameter :: tFilter(4) = [3,4,5,8]         ! CTAC   - trajectory
!integer(4),parameter :: tFilter(6) = [3,4,5,7,8,9] ! diCTAC - trajectory

real(4),parameter ::    dwdt = 0.25 ! size of bins in distributions

! === initiation of other variables ================================= 
integer(4) :: n_fol
integer(4) :: n_xtac

integer(4) i,j,k,ix,iy,kx,ky,n,np,istat,itimestep,idd(nnp),nbw(2),names(nn),dir,io
integer(4) kxix,kyiy

real(4) r(3,nn),rp(3,nnp),Lb(3),Lbl(3),Lbh(3),xij,yij,rij

real(4) surfw(2,-nbwMax:nbwMax,-nbwMax:nbwMax),zsurf(2),widthDist(0:nb)

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

real(8) Dove_FOL_dist(0:nb),Dove_XTAC_dist(0:nb)
real(8) Lz_FOL_dist(0:nb),Lz_XTAC_dist(0:nb)
real(8) Dc_FOL_dist(0:nb),Dc_XTAC_dist(0:nb)
real(8) N_FOL_dist(0:nb),N_XTAC_dist(0:nb)

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
cnfg = get_snap()
call patts%find_patts(cnfg%tp, [3,4,5], [".","+","."] )  ! FOL
call patts%find_patts(cnfg%tp, [8,4,5], [".","+","."] ) ! CTAC
!call patts%find_patts(cnfg%tp, [7,8,9], [".","+","."] )  ! diCTAC tails
!call patts%find_patts(cnfg%tp, [6], ["."] ) ! diCTAC head
    
if(size(patts%patts) < 2) stop "size( patts ) < 2"

call close_trajectory()

! === set values of some variables ==================================
n_xtac = patts%patts(2)%n
n_fol = patts%patts(1)%n

widthDist = 0.0

N_FOL_dist=0.0
N_XTAC_dist=0.0

Dc_FOL_dist=0.0
Dc_XTAC_dist=0.0

Lz_FOL_dist=0.0
Lz_XTAC_dist=0.0

Dove_FOL_dist=0.0
Dove_XTAC_dist=0.0

! === the main loop through the trajectories ========================
open(newunit=io,file=trim(file_name),action='read',status='old')
istat=0
n=0
itimestep=0
do
  ! === read the head of one configuration ==========================
  read(io,'(a)',iostat=istat) blabla
  if(istat/=0) exit
  read(io,'(a)',iostat=istat) blabla
  if(istat/=0) exit
  read(io,'(a)',iostat=istat) blabla
  if(istat/=0) exit
  read(io,*,iostat=istat) n
  if(istat/=0) exit
  if(n>nn) stop "n>nn"
  read(io,'(a)',iostat=istat) blabla
  if(istat/=0) exit  
  read(io,*,iostat=istat) Lbl(1),Lbh(1)
  read(io,*,iostat=istat) Lbl(2),Lbh(2)
  read(io,*,iostat=istat) Lbl(3),Lbh(3)
  Lb = Lbh - Lbl
  if(itimestep==0) then
    nbw(1) = int(0.5*lb(1))
    nbw(2) = int(0.5*lb(2))
    write(*,'(a,3f15.10)') "Lb ::: ",Lb(:)
    write(*,'(a,2i15)') "nbw ::: ",nbw(:)
    if(nbw(1)>nbwMax) stop "nbw(1)>nbwMax"
    if(nbw(2)>nbwMax) stop "nbw(2)>nbwMax"
  endif
  if(istat/=0) exit  
  read(io,'(a)',iostat=istat) blabla
  if(istat/=0) exit  
  
      ! === read atom positions =====================================
       do i=1,n
         read(io,'(a)',iostat=istat) blabla
         if(istat/=0) exit
         
         read(blabla,*) j
         read(blabla,*) k,names(j),r(:,j)
         
       enddo 
       if(istat/=0) exit
      
        itimestep=itimestep+1
        
        write(*,*) "log: iTimeStep=",itimestep
        
        ! === zeros accumulated variables ===========================
        
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
        
        ! === filtering of atoms ====================================
        np=0
        rp=0.0
        idd=0
        do i=1,n

          if(any(names(i)==tFilter)) then
            np=np+1

            rp(:,np)=r(:,i)
            idd(np)=i
          endif
        enddo
        
        write(*,*) "log: Np=",np
        
        ! === the loop through a one configuration ======================
        surfw(1,:,:) = -huge(1.)
        surfw(2,:,:) = huge(1.)
        i=1
        do while(i<=np)
          ! === if the next sequence starts the same as the pattern(1)
          if(names(idd(i))==patts%patts(1)%patt(1)) then
              
            ! === ITIM method - surface search ==========================
            do j=0,n_fol-1

              kx=nint(rp(1,i+j))
              ky=nint(rp(2,i+j))
              
              do ix=-2,2
                kxix = kx+ix
                if(kxix<-nbw(1)) kxix =  2*nbw(1)+kxix+1
                if(kxix> nbw(1)) kxix = -2*nbw(1)+kxix-1
                do iy=-2,2
                  kyiy = ky+iy
                  if(kyiy<-nbw(2)) kyiy =  2*nbw(2)+kyiy+1
                  if(kyiy> nbw(2)) kyiy = -2*nbw(2)+kyiy-1
                  
                  xij = real(kxix)-rp(1,i+j)
                  yij = real(kyiy)-rp(2,i+j)
                  xij = xij - Lb(1)*nint(xij/Lb(1))
                  yij = yij - Lb(2)*nint(yij/Lb(2))
                  rij = sqrt(xij**2 + yij**2)
                  if(rij>1.d0) cycle
                  rij = sqrt(1.0 - rij**2)
                  zsurf(1) = rp(3,i+j) + rij - 0.5  
                  zsurf(2) = rp(3,i+j) - rij + 0.5 

                  if(zsurf(1) > surfw(1,kxix,kyiy))  surfw(1,kxix,kyiy) = zsurf(1)
                  if(zsurf(2) < surfw(2,kxix,kyiy))  surfw(2,kxix,kyiy) = zsurf(2)
                  
                enddo
              enddo
            enddo

            ! === compute values on the lattice =====================
            j=j-1
            dir = get_direction(rp(3,i+j), rp(3,i))

            kx=nint(rp(1,i))
            ky=nint(rp(2,i))
            if(abs(ix) > nbwMax) stop
            if(abs(iy) > nbwMax) stop
             
              zt1_FOL(dir,kx,ky)=zt1_FOL(dir,kx,ky)+rp(3,i+j)
              Nzt1_FOL(dir,kx,ky)=Nzt1_FOL(dir,kx,ky)+1.0
              ztn_FOL(dir,kx,ky)=ztn_FOL(dir,kx,ky)+rp(3,i+1)
              Nztn_FOL(dir,kx,ky)=Nztn_FOL(dir,kx,ky)+1.0
                            
              if(dir==1) then
                z1t1_FOL(kx,ky)=z1t1_FOL(kx,ky)+rp(3,i+1)
                Nz1t1_FOL(kx,ky)=Nz1t1_FOL(kx,ky)+1.0
              else
                z2t1_FOL(kx,ky)=z2t1_FOL(kx,ky)+rp(3,i+1)
                Nz2t1_FOL(kx,ky)=Nz2t1_FOL(kx,ky)+1.0
              endif
               
            i=i+n_fol
          ! === if the next sequence starts the same as the pattern(2)
          elseif(names(idd(i))==patts%patts(2)%patt(1)) then
            do j=0,n_xtac-1

              kx=nint(rp(1,i+j))
              ky=nint(rp(2,i+j))
    
              do ix=-2,2
                kxix = kx+ix
                if(kxix<-nbw(1)) kxix =  2*nbw(1)+kxix+1
                if(kxix> nbw(1)) kxix = -2*nbw(1)+kxix-1
                do iy=-2,2
                  kyiy = ky+iy
                  if(kyiy<-nbw(2)) kyiy =  2*nbw(2)+kyiy+1
                  if(kyiy> nbw(2)) kyiy = -2*nbw(2)+kyiy-1

                  xij = real(kxix)-rp(1,i+j)
                  yij = real(kyiy)-rp(2,i+j)
                  xij = xij - Lb(1)*nint(xij/Lb(1))
                  yij = yij - Lb(2)*nint(yij/Lb(2))
                  rij = sqrt(xij**2 + yij**2)
                  if(rij>1.d0) cycle
                  rij = sqrt(1.0 - rij**2)
                  zsurf(1) = rp(3,i+j) + rij - 0.5  
                  zsurf(2) = rp(3,i+j) - rij + 0.5 

                  if(zsurf(1) > surfw(1,kxix,kyiy))  surfw(1,kxix,kyiy) = zsurf(1)
                  if(zsurf(2) < surfw(2,kxix,kyiy))  surfw(2,kxix,kyiy) = zsurf(2)
                  
                enddo
              enddo
                            
            enddo

            ! === compute values on the lattice =====================
            j=j-1
            dir = get_direction(rp(3,i+j), rp(3,i))
          
            kx=nint(rp(1,i))
            ky=nint(rp(2,i))
            
            zt1_XTAC(dir,kx,ky)=zt1_XTAC(dir,kx,ky)+rp(3,i+j)
            Nzt1_XTAC(dir,kx,ky)=Nzt1_XTAC(dir,kx,ky)+1.0
            ztn_XTAC(dir,kx,ky)=ztn_XTAC(dir,kx,ky)+rp(3,i+1)
            Nztn_XTAC(dir,kx,ky)=Nztn_XTAC(dir,kx,ky)+1.0
                          
            if(dir==1) then
              z1t1_XTAC(kx,ky)=z1t1_XTAC(kx,ky)+rp(3,i+1)
              Nz1t1_XTAC(kx,ky)=Nz1t1_XTAC(kx,ky)+1.0
            else
              z2t1_XTAC(kx,ky)=z2t1_XTAC(kx,ky)+rp(3,i+1)
              Nz2t1_XTAC(kx,ky)=Nz2t1_XTAC(kx,ky)+1.0
            endif
            
            
            i=i+n_xtac
          else
            i=i+1
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
            
            rij=abs(surfw(1,ix,iy)-surfw(2,ix,iy))
            k = nint(rij/dwdt)
            if(k>nb) cycle
            widthDist(k) = widthDist(k) + 1.0        
             
             if((.not. ieee_is_nan(Lz_FOL(ix,iy))) .and. (.not. ieee_is_nan(Dc_FOL(ix,iy)))) then
               Dove_FOL_dist(k) = Dove_FOL_dist(k) + ((2.d0*Lz_FOL(ix,iy)-Dc_FOL(ix,iy))/Lz_FOL(ix,iy))
               Lz_FOL_dist(k) = Lz_FOL_dist(k) + Lz_FOL(ix,iy)
               Dc_FOL_dist(k) = Dc_FOL_dist(k) + Dc_FOL(ix,iy)
               N_FOL_dist(k) = N_FOL_dist(k) + 1.d0
             endif
             
             if((.not. ieee_is_nan(Lz_XTAC(ix,iy))) .and. (.not. ieee_is_nan(Dc_XTAC(ix,iy)))) then
               Dove_XTAC_dist(k) = Dove_XTAC_dist(k) + ((2.d0*Lz_XTAC(ix,iy)-Dc_XTAC(ix,iy))/Lz_XTAC(ix,iy))
               Lz_XTAC_dist(k) = Lz_XTAC_dist(k) + Lz_XTAC(ix,iy)
               Dc_XTAC_dist(k) = Dc_XTAC_dist(k) + Dc_XTAC(ix,iy)
               N_XTAC_dist(k) = N_XTAC_dist(k) + 1.d0
             endif
            
          enddo
        enddo
      
        !if(itimestep>10) exit
    if(istat/=0) exit
   
enddo

write(*,'(2a)') "log: <<< end >>>"

close(io)

! === the normalisation of the results ==============================
Dc_FOL_dist=Dc_FOL_dist/N_FOL_dist
Dc_XTAC_dist=Dc_XTAC_dist/N_XTAC_dist

Lz_FOL_dist=Lz_FOL_dist/N_FOL_dist
Lz_XTAC_dist=Lz_XTAC_dist/N_XTAC_dist

Dove_FOL_dist=Dove_FOL_dist/N_FOL_dist
Dove_XTAC_dist=Dove_XTAC_dist/N_XTAC_dist

widthDist = widthDist / (sum(widthDist)*dwdt)

! === make the result directory if doesn't exist ===================
inquire(file=trim(charus),exist=existuje)
  if(.not. existuje) then
    !stop "neexistuje"
    call execute_command_line("mkdir "//trim(charus))
  endif
  
! === print the output files ========================================
open(newunit=io,file=trim(charus)//"/N_width_dist.data")
write(io,*) "width N_{FOL8} N_{CTAB}"
do i=0,nb
    write(io,*) real(i)*dwdt,N_FOL_dist(i)/dble(itimestep),N_XTAC_dist(i)/dble(itimestep)
enddo
close(io)

open(newunit=io,file=trim(charus)//"/Doverlap_width_dist.data")
write(io,*) "width Doverlap_{FOL8} Doverlap_{CTAB}"
do i=0,nb
  if(widthDist(i) > 0.01) then
    write(io,*) real(i)*dwdt,Dove_FOL_dist(i),Dove_XTAC_dist(i)
  else
    write(io,*) real(i)*dwdt,"NaN ","NaN ","NaN"
  endif
enddo
close(io)

open(newunit=io,file=trim(charus)//"/Lz_width_dist.data")
write(io,*) "width Lz_{FOL8} Lz_{CTAB}"
do i=0,nb
  if(widthDist(i) > 0.01) then
    write(io,*) real(i)*dwdt,Lz_FOL_dist(i),Lz_XTAC_dist(i)
  else
    write(io,*) real(i)*dwdt,"NaN ","NaN ","NaN"
  endif
enddo
close(io)

open(newunit=io,file=trim(charus)//"/Dc_width_dist.data")
write(io,*) "width Dc_{FOL8} Dc_{CTAB}"
do i=0,nb
  if(widthDist(i) > 0.01) then
    write(io,*) real(i)*dwdt,Dc_FOL_dist(i),Dc_XTAC_dist(i)
  else
    write(io,*) real(i)*dwdt,"NaN ","NaN ","NaN"
  endif
enddo
close(io)

! === plot the output files by Gnuplot ==============================
call gnupl%init()
call gnupl%set_terminal("pngcairo",&
'enhanced size 800,600 font "Helvetica, 12" linewidth 2')

call gnupl%set_dir(trim(charus))
call gnupl%set("key","autotitle columnhead")

call gnupl%set("xlabel","'width'")

call gnupl%rm_data()
call gnupl%add('Doverlap_width_dist.data',u="1:2")
call gnupl%add('Doverlap_width_dist.data',u="1:3")
  call gnupl%plot("Doverlap_width_dist.png", yl="Doverlap")

call gnupl%rm_data()
call gnupl%add('Dc_width_dist.data',u="1:2")
call gnupl%add('Dc_width_dist.data',u="1:3")
  call gnupl%plot("Dc_width_dist.png", yl="Dc")

call gnupl%rm_data()
call gnupl%add('Lz_width_dist.data',u="1:2")
call gnupl%add('Lz_width_dist.data',u="1:3")
  call gnupl%plot("Lz_width_dist.png", yl="Lz")

call gnupl%rm_data()
call gnupl%add('N_width_dist.data',u="1:2")
call gnupl%add('N_width_dist.data',u="1:3")
  call gnupl%plot("N_width_dist.png", yl="N")

!====================================================================
  contains
  !--------------------------------------------------------------------
  !
  ! The function gets the direction of the polymer.
  ! Gets 1 or 2 for up or down. 
  !
  pure real(4) function get_direction(z1, z2) result(dir)
      real(4),intent(in) :: z1, z2
      if(z2-z1 > 0.d0) then
        dir = 1
      else
        dir = 2
      endif
  end function
end program

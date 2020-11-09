program DolDcLzNWidth
use ieee_arithmetic
use config_mod
use gnuplot_mod
use lmp_traject_reader_mod
use mol_pattern_mod
implicit none

integer(4),parameter :: nn=180000
integer(4),parameter :: nnp=50000
integer(4),parameter :: nb=200
integer(4),parameter :: nbwMax=20
integer(4) :: n_fol
integer(4) :: n_xtac

integer(4) i,j,k,ix,iy,kx,ky,n,np,istat,itimestep,idd(nnp),nbw(2),names(nn),dir
integer(4) kxix,kyiy
integer(4),allocatable :: tFilter(:)

real(4) r(3,nn),rp(3,nnp),Lb(3),Lbl(3),Lbh(3),xij,yij,rij

real(4) surfw(2,-nbwMax:nbwMax,-nbwMax:nbwMax),zsurf(2),widthDist(0:nb),dwdt


real(4) Lz_FOL8(-nbwMax:nbwMax,-nbwMax:nbwMax),Dc_FOL8(-nbwMax:nbwMax,-nbwMax:nbwMax),&
&z1t1_FOL8(-nbwMax:nbwMax,-nbwMax:nbwMax)
real(4) z2t1_FOL8(-nbwMax:nbwMax,-nbwMax:nbwMax),zt1_FOL8(2,-nbwMax:nbwMax,-nbwMax:nbwMax),&
&ztn_FOL8(2,-nbwMax:nbwMax,-nbwMax:nbwMax)
real(4) Lz_CTAB(-nbwMax:nbwMax,-nbwMax:nbwMax),Dc_CTAB(-nbwMax:nbwMax,-nbwMax:nbwMax),&
&z1t1_CTAB(-nbwMax:nbwMax,-nbwMax:nbwMax)
real(4) z2t1_CTAB(-nbwMax:nbwMax,-nbwMax:nbwMax),zt1_CTAB(2,-nbwMax:nbwMax,-nbwMax:nbwMax),&
&ztn_CTAB(2,-nbwMax:nbwMax,-nbwMax:nbwMax)
real(4) Nz1t1_FOL8(-nbwMax:nbwMax,-nbwMax:nbwMax),Nz2t1_FOL8(-nbwMax:nbwMax,-nbwMax:nbwMax),&
&Nzt1_FOL8(2,-nbwMax:nbwMax,-nbwMax:nbwMax)
real(4) Nztn_FOL8(2,-nbwMax:nbwMax,-nbwMax:nbwMax)
real(4) Nz1t1_CTAB(-nbwMax:nbwMax,-nbwMax:nbwMax),Nz2t1_CTAB(-nbwMax:nbwMax,-nbwMax:nbwMax),&
&Nzt1_CTAB(2,-nbwMax:nbwMax,-nbwMax:nbwMax)
real(4) Nztn_CTAB(2,-nbwMax:nbwMax,-nbwMax:nbwMax)

real(8) Dove_FOL8_dist(0:nb),Dove_CTAB_dist(0:nb)
real(8) Lz_FOL8_dist(0:nb),Lz_CTAB_dist(0:nb)
real(8) Dc_FOL8_dist(0:nb),Dc_CTAB_dist(0:nb)
real(8) N_FOL8_dist(0:nb),N_CTAB_dist(0:nb)

character(len=500) blabla
character(len=100) charus, file_name

integer(4) num_of_args

logical existuje

type(mol_patts) :: patts
type(config) :: cnfg
type(gnuplot) :: gnupl

num_of_args = command_argument_count()
if(num_of_args < 2) stop "num_of_args < 2"

call get_command_argument(1, charus)
call get_command_argument(2, file_name)

write(*,*) "charus: ",trim(charus)
write(*,*) "file name: ",trim(file_name)

tFilter = [3,4,5,8]    ! CTAC   - trajectory
!tFilter = [3,4,5,7,8,9] ! diCTAC - trajectory

call open_trajectory(trim(file_name))
cnfg = get_snap()
call patts%find_patts(cnfg%tp, [3,4,5], [".","+","."] )  ! FOL
call patts%find_patts(cnfg%tp, [8,4,5], [".","+","."] ) ! CTAC
!call patts%find_patts(cnfg%tp, [7,8,9], [".","+","."] )  ! diCTAC tails
!call patts%find_patts(cnfg%tp, [6], ["."] ) ! diCTAC head
    
if(size(patts%patts) < 2) stop "size( patts ) < 2"

n_xtac = patts%patts(2)%n
n_fol = patts%patts(1)%n

call close_trajectory()

!charus = "_T_0_7_min_stress"

dwdt = 0.25

!open(1,file="Teq0_7_stress_min_half/nvt/lam_vt.cs",action='read',status='old')
open(1,file=trim(file_name),action='read',status='old')
!open(1,file="vtf_overview.txt",action='read',status='old')

widthDist = 0.0

N_FOL8_dist=0.0
N_CTAB_dist=0.0

Dc_FOL8_dist=0.0
Dc_CTAB_dist=0.0

Lz_FOL8_dist=0.0
Lz_CTAB_dist=0.0

Dove_FOL8_dist=0.0
Dove_CTAB_dist=0.0

istat=0
n=0
itimestep=0
do
  read(1,'(a)',iostat=istat) blabla
  if(istat/=0) exit
  read(1,'(a)',iostat=istat) blabla
  if(istat/=0) exit
  read(1,'(a)',iostat=istat) blabla
  if(istat/=0) exit
  read(1,*,iostat=istat) n
  if(istat/=0) exit
  if(n>nn) stop "n>nn"
  read(1,'(a)',iostat=istat) blabla
  if(istat/=0) exit  
  read(1,*,iostat=istat) Lbl(1),Lbh(1)
  read(1,*,iostat=istat) Lbl(2),Lbh(2)
  read(1,*,iostat=istat) Lbl(3),Lbh(3)
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
  read(1,'(a)',iostat=istat) blabla
  if(istat/=0) exit  
  
       do i=1,n
         read(1,'(a)',iostat=istat) blabla
         if(istat/=0) exit
         
         read(blabla,*) j
         read(blabla,*) k,names(j),r(:,j)
         
       enddo 
       if(istat/=0) exit
      
      
        itimestep=itimestep+1
        
        write(*,*) "log: iTimeStep=",itimestep
        
        !analyza ///////////////////////////////////////////////////////
        
        z1t1_FOL8 = 0.0
        Nz1t1_FOL8 = 0.0
        z2t1_FOL8 = 0.0
        Nz2t1_FOL8 = 0.0
        zt1_FOL8 = 0.0
        Nzt1_FOL8 = 0.0
        ztn_FOL8 = 0.0
        Nztn_FOL8 = 0.0

        z1t1_CTAB = 0.0
        Nz1t1_CTAB = 0.0
        z2t1_CTAB = 0.0
        Nz2t1_CTAB = 0.0
        zt1_CTAB = 0.0
        Nzt1_CTAB = 0.0
        ztn_CTAB = 0.0
        Nztn_CTAB = 0.0
        
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
        
        surfw(1,:,:) = -huge(1.)
        surfw(2,:,:) = huge(1.)
        i=1
        do while(i<=np)
          if(names(idd(i))==3) then
              
            do j=0,n_fol-1

              kx=nint(rp(1,i+j))
              ky=nint(rp(2,i+j))
              
              do ix=-2,2
                kxix = kx+ix
                !if(kxix<0 .or. kx+ix>nbwMax) stop "bum"
                if(kxix<-nbw(1)) kxix =  2*nbw(1)+kxix+1
                if(kxix> nbw(1)) kxix = -2*nbw(1)+kxix-1
                do iy=-2,2
                  kyiy = ky+iy
                  !if(ky+iy<0 .or. ky+iy>nbwMax) stop "bum"
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

            j=j-1
            
             if(rp(3,i)-rp(3,i+j) > 0.d0) then
              dir = 1
            else
              dir = 2
            endif

            kx=nint(rp(1,i))
            ky=nint(rp(2,i))
            if(abs(ix) > nbwMax) stop
            if(abs(iy) > nbwMax) stop
             
              zt1_FOL8(dir,kx,ky)=zt1_FOL8(dir,kx,ky)+rp(3,i+j)
              Nzt1_FOL8(dir,kx,ky)=Nzt1_FOL8(dir,kx,ky)+1.0
              ztn_FOL8(dir,kx,ky)=ztn_FOL8(dir,kx,ky)+rp(3,i+1)
              Nztn_FOL8(dir,kx,ky)=Nztn_FOL8(dir,kx,ky)+1.0
                            
              if(dir==1) then
                z1t1_FOL8(kx,ky)=z1t1_FOL8(kx,ky)+rp(3,i+1)
                Nz1t1_FOL8(kx,ky)=Nz1t1_FOL8(kx,ky)+1.0
              else
                z2t1_FOL8(kx,ky)=z2t1_FOL8(kx,ky)+rp(3,i+1)
                Nz2t1_FOL8(kx,ky)=Nz2t1_FOL8(kx,ky)+1.0
              endif
               
            i=i+n_fol
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

            j=j-1
            
           if(rp(3,i)-rp(3,i+j) > 0.d0) then
              dir = 1
            else
              dir = 2
            endif
          
            kx=nint(rp(1,i))
            ky=nint(rp(2,i))
            
            zt1_CTAB(dir,kx,ky)=zt1_CTAB(dir,kx,ky)+rp(3,i+j)
            Nzt1_CTAB(dir,kx,ky)=Nzt1_CTAB(dir,kx,ky)+1.0
            ztn_CTAB(dir,kx,ky)=ztn_CTAB(dir,kx,ky)+rp(3,i+1)
            Nztn_CTAB(dir,kx,ky)=Nztn_CTAB(dir,kx,ky)+1.0
                          
            if(dir==1) then
              z1t1_CTAB(kx,ky)=z1t1_CTAB(kx,ky)+rp(3,i+1)
              Nz1t1_CTAB(kx,ky)=Nz1t1_CTAB(kx,ky)+1.0
            else
              z2t1_CTAB(kx,ky)=z2t1_CTAB(kx,ky)+rp(3,i+1)
              Nz2t1_CTAB(kx,ky)=Nz2t1_CTAB(kx,ky)+1.0
            endif
            
            
            i=i+n_xtac
          else
            i=i+1
          endif
        enddo

        z1t1_FOL8 = z1t1_FOL8/Nz1t1_FOL8
        z2t1_FOL8 = z2t1_FOL8/Nz2t1_FOL8
        zt1_FOL8 = zt1_FOL8/Nzt1_FOL8
        ztn_FOL8 = ztn_FOL8/Nztn_FOL8

        z1t1_CTAB = z1t1_CTAB/Nz1t1_CTAB
        z2t1_CTAB = z2t1_CTAB/Nz2t1_CTAB
        zt1_CTAB = zt1_CTAB/Nzt1_CTAB
        ztn_CTAB = ztn_CTAB/Nztn_CTAB

        Dc_FOL8 = abs(z1t1_FOL8 - z2t1_FOL8)
        Dc_CTAB = abs(z1t1_CTAB - z2t1_CTAB)

        Lz_FOL8 = 0.5*(abs(zt1_FOL8(1,:,:) - ztn_FOL8(1,:,:))+abs(zt1_FOL8(2,:,:) - ztn_FOL8(2,:,:)))
        Lz_CTAB = 0.5*(abs(zt1_CTAB(1,:,:) - ztn_CTAB(1,:,:))+abs(zt1_CTAB(2,:,:) - ztn_CTAB(2,:,:)))
        
        

        do ix=-nbw(1),nbw(1)
          do iy=-nbw(2),nbw(2)
            
            rij=abs(surfw(1,ix,iy)-surfw(2,ix,iy))
            k = nint(rij/dwdt)
            if(k>nb) cycle
            widthDist(k) = widthDist(k) + 1.0        
             
             if((.not. ieee_is_nan(Lz_FOL8(ix,iy))) .and. (.not. ieee_is_nan(Dc_FOL8(ix,iy)))) then
               Dove_FOL8_dist(k) = Dove_FOL8_dist(k) + ((2.d0*Lz_FOL8(ix,iy)-Dc_FOL8(ix,iy))/Lz_FOL8(ix,iy))
               Lz_FOL8_dist(k) = Lz_FOL8_dist(k) + Lz_FOL8(ix,iy)
               Dc_FOL8_dist(k) = Dc_FOL8_dist(k) + Dc_FOL8(ix,iy)
               N_FOL8_dist(k) = N_FOL8_dist(k) + 1.d0
             endif
             
             if((.not. ieee_is_nan(Lz_CTAB(ix,iy))) .and. (.not. ieee_is_nan(Dc_CTAB(ix,iy)))) then
               Dove_CTAB_dist(k) = Dove_CTAB_dist(k) + ((2.d0*Lz_CTAB(ix,iy)-Dc_CTAB(ix,iy))/Lz_CTAB(ix,iy))
               Lz_CTAB_dist(k) = Lz_CTAB_dist(k) + Lz_CTAB(ix,iy)
               Dc_CTAB_dist(k) = Dc_CTAB_dist(k) + Dc_CTAB(ix,iy)
               N_CTAB_dist(k) = N_CTAB_dist(k) + 1.d0
             endif
            
            
          enddo
        enddo
        !analyza ///////////////////////////////////////////////////////
      
        !if(itimestep>10) exit
    if(istat/=0) exit
   
enddo

write(*,'(2a)') "log: <<< end >>>"

close(1)

Dc_FOL8_dist=Dc_FOL8_dist/N_FOL8_dist
Dc_CTAB_dist=Dc_CTAB_dist/N_CTAB_dist

Lz_FOL8_dist=Lz_FOL8_dist/N_FOL8_dist
Lz_CTAB_dist=Lz_CTAB_dist/N_CTAB_dist

Dove_FOL8_dist=Dove_FOL8_dist/N_FOL8_dist
Dove_CTAB_dist=Dove_CTAB_dist/N_CTAB_dist

widthDist = widthDist / (sum(widthDist)*dwdt)

inquire(file=trim(folName(charus)),exist=existuje)
  if(.not. existuje) then
    !stop "neexistuje"
    call execute_command_line("mkdir "//trim(folName(charus)))
  endif
  
  
!///////////////////////////

open(1,file=trim(folName(charus))//"/N_width_dist.data")
write(1,*) "width N_{FOL8} N_{CTAB}"
do i=0,nb
    write(1,*) real(i)*dwdt,N_FOL8_dist(i)/dble(itimestep),N_CTAB_dist(i)/dble(itimestep)
enddo
close(1)

open(1,file=trim(folName(charus))//"/Doverlap_width_dist.data")
write(1,*) "width Doverlap_{FOL8} Doverlap_{CTAB}"
do i=0,nb
  if(widthDist(i) > 0.01) then
    write(1,*) real(i)*dwdt,Dove_FOL8_dist(i),Dove_CTAB_dist(i)
  else
    write(1,*) real(i)*dwdt,"NaN ","NaN ","NaN"
  endif
enddo
close(1)

open(1,file=trim(folName(charus))//"/Lz_width_dist.data")
write(1,*) "width Lz_{FOL8} Lz_{CTAB}"
do i=0,nb
  if(widthDist(i) > 0.01) then
    write(1,*) real(i)*dwdt,Lz_FOL8_dist(i),Lz_CTAB_dist(i)
  else
    write(1,*) real(i)*dwdt,"NaN ","NaN ","NaN"
  endif
enddo
close(1)

open(1,file=trim(folName(charus))//"/Dc_width_dist.data")
write(1,*) "width Dc_{FOL8} Dc_{CTAB}"
do i=0,nb
  if(widthDist(i) > 0.01) then
    write(1,*) real(i)*dwdt,Dc_FOL8_dist(i),Dc_CTAB_dist(i)
  else
    write(1,*) real(i)*dwdt,"NaN ","NaN ","NaN"
  endif
enddo
close(1)


call gnupl%init()
call gnupl%set_terminal("pngcairo",&
'enhanced size 800,600 font "Helvetica, 12" linewidth 2')

call gnupl%set_dir("result"//trim(charus))
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

contains
!----------------------------------------------------------
character(len=50) function folName(charus) result(foldName)

character(len=*),intent(in) :: charus

foldName = "result"//trim(charus)

end function
end program

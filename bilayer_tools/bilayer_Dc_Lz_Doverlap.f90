program DcLzDoverlap
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
integer(4),allocatable :: tFilter(:)

integer(4) i,j,k,ix,iy,kx,ky,n,istat,itimestep,nbw(2),dir, io

real(4) r(3,nn),Lb(3),Lbl(3),Lbh(3)

real(4) dd

real(4) Lz_FOL8(-nbwMax:nbwMax,-nbwMax:nbwMax),Dc_FOL8(-nbwMax:nbwMax,-nbwMax:nbwMax),&
&z1t1_FOL8(-nbwMax:nbwMax,-nbwMax:nbwMax)
real(4) z2t1_FOL8(-nbwMax:nbwMax,-nbwMax:nbwMax),zt1_FOL8(2,-nbwMax:nbwMax,-nbwMax:nbwMax),&
&ztn_FOL8(2,-nbwMax:nbwMax,-nbwMax:nbwMax)
real(4) Lz_BTAC(-nbwMax:nbwMax,-nbwMax:nbwMax),Dc_BTAC(-nbwMax:nbwMax,-nbwMax:nbwMax),&
&z1t1_BTAC(-nbwMax:nbwMax,-nbwMax:nbwMax)
real(4) z2t1_BTAC(-nbwMax:nbwMax,-nbwMax:nbwMax),zt1_BTAC(2,-nbwMax:nbwMax,-nbwMax:nbwMax),&
&ztn_BTAC(2,-nbwMax:nbwMax,-nbwMax:nbwMax)
real(4) Nz1t1_FOL8(-nbwMax:nbwMax,-nbwMax:nbwMax),Nz2t1_FOL8(-nbwMax:nbwMax,-nbwMax:nbwMax),&
&Nzt1_FOL8(2,-nbwMax:nbwMax,-nbwMax:nbwMax)
real(4) Nztn_FOL8(2,-nbwMax:nbwMax,-nbwMax:nbwMax)
real(4) Nz1t1_BTAC(-nbwMax:nbwMax,-nbwMax:nbwMax),Nz2t1_BTAC(-nbwMax:nbwMax,-nbwMax:nbwMax),&
&Nzt1_BTAC(2,-nbwMax:nbwMax,-nbwMax:nbwMax)
real(4) Nztn_BTAC(2,-nbwMax:nbwMax,-nbwMax:nbwMax)

real(8) Dove_FOL8_dist(-nb:nb),Dove_BTAC_dist(-nb:nb)
real(8) Lz_FOL8_dist(-nb:nb),Lz_BTAC_dist(-nb:nb)
real(8) Dc_FOL8_dist(-nb:nb),Dc_BTAC_dist(-nb:nb)
real(4) ddove,ddc,dlz

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
cnfg = get_snap(tFilter)
call patts%find_patts(cnfg%tp, [3,4,5], [".","+","."] ) ! FOL
call patts%find_patts(cnfg%tp, [8,4,5], [".","+","."] ) ! CTAC
!call patts%find_patts(cnfg%tp, [7,8,9], [".","+","."] ) !diCTAC tails
    
if(size(patts%patts) < 2) stop "size( patts ) < 2"

n_xtac = patts%patts(2)%n
n_fol = patts%patts(1)%n

ddove = 0.05
ddc = 0.05
dlz = 0.025

dd = 2.5


Lb = cnfg%Lbox

nbw(1) = int(0.5*lb(1)/dd)
nbw(2) = int(0.5*lb(2)/dd)
write(*,'(a,3f15.10)') "Lb ::: ",Lb(:)
write(*,'(a,2i15)') "nbw ::: ",nbw(:)
if(nbw(1)>nbwMax) stop "nbw(1)>nbwMax"
if(nbw(2)>nbwMax) stop "nbw(2)>nbwMax"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call reopen_trajectory()

Dc_FOL8_dist=0.0
Dc_BTAC_dist=0.0

Lz_FOL8_dist=0.0
Lz_BTAC_dist=0.0

Dove_FOL8_dist=0.0
Dove_BTAC_dist=0.0

cnfg = get_snap(tFilter)
do while(.not. is_end)

    write(*,*) cur_step, step_num, cnfg%n

    z1t1_FOL8 = 0.0
    Nz1t1_FOL8 = 0.0
    z2t1_FOL8 = 0.0
    Nz2t1_FOL8 = 0.0
    zt1_FOL8 = 0.0
    Nzt1_FOL8 = 0.0
    ztn_FOL8 = 0.0
    Nztn_FOL8 = 0.0

    z1t1_BTAC = 0.0
    Nz1t1_BTAC = 0.0
    z2t1_BTAC = 0.0
    Nz2t1_BTAC = 0.0
    zt1_BTAC = 0.0
    Nzt1_BTAC = 0.0
    ztn_BTAC = 0.0
    Nztn_BTAC = 0.0

    i = 1
    do while(i<=cnfg%n)
      if(cnfg%tp(i) == patts%patts(1)%patt(1)) then

          j=n_fol-1
            
          if(cnfg%r(3,i)-cnfg%r(3,i+j) > 0.d0) then
            dir = 1
          else
            dir = 2
          endif

          kx=nint(cnfg%r(1,i)/dd)
          ky=nint(cnfg%r(2,i)/dd)
          
          zt1_FOL8(dir,kx,ky)=zt1_FOL8(dir,kx,ky)+cnfg%r(3,i+j)
          Nzt1_FOL8(dir,kx,ky)=Nzt1_FOL8(dir,kx,ky)+1.0
          ztn_FOL8(dir,kx,ky)=ztn_FOL8(dir,kx,ky)+cnfg%r(3,i+1)
          Nztn_FOL8(dir,kx,ky)=Nztn_FOL8(dir,kx,ky)+1.0
                          
          if(dir==1) then
            z1t1_FOL8(kx,ky)=z1t1_FOL8(kx,ky)+cnfg%r(3,i+1)
            Nz1t1_FOL8(kx,ky)=Nz1t1_FOL8(kx,ky)+1.0
          else
            z2t1_FOL8(kx,ky)=z2t1_FOL8(kx,ky)+cnfg%r(3,i+1)
            Nz2t1_FOL8(kx,ky)=Nz2t1_FOL8(kx,ky)+1.0
          endif
          
          i=i+n_fol

          i = i + patts%patts(1)%n
      elseif(cnfg%tp(i) == patts%patts(2)%patt(1)) then

          j=n_xtac-1    
          
          if(cnfg%r(3,i)-cnfg%r(3,i+j) > 0.d0) then
            dir = 1
          else
            dir = 2
          endif
        
          kx=nint(cnfg%r(1,i)/dd)
          ky=nint(cnfg%r(2,i)/dd)
          
          zt1_BTAC(dir,kx,ky)=zt1_BTAC(dir,kx,ky)+cnfg%r(3,i+j)
          Nzt1_BTAC(dir,kx,ky)=Nzt1_BTAC(dir,kx,ky)+1.0
          ztn_BTAC(dir,kx,ky)=ztn_BTAC(dir,kx,ky)+cnfg%r(3,i+1)
          Nztn_BTAC(dir,kx,ky)=Nztn_BTAC(dir,kx,ky)+1.0
                        
          if(dir==1) then
            z1t1_BTAC(kx,ky)=z1t1_BTAC(kx,ky)+cnfg%r(3,i+1)
            Nz1t1_BTAC(kx,ky)=Nz1t1_BTAC(kx,ky)+1.0
          else
            z2t1_BTAC(kx,ky)=z2t1_BTAC(kx,ky)+cnfg%r(3,i+1)
            Nz2t1_BTAC(kx,ky)=Nz2t1_BTAC(kx,ky)+1.0
          endif
          
          i=i+n_xtac

          i = i + patts%patts(2)%n
      else
          !stop "???"
          i = i + 1
          cycle
      endif
    enddo
        
        
    z1t1_FOL8 = z1t1_FOL8/Nz1t1_FOL8
    z2t1_FOL8 = z2t1_FOL8/Nz2t1_FOL8
    zt1_FOL8 = zt1_FOL8/Nzt1_FOL8
    ztn_FOL8 = ztn_FOL8/Nztn_FOL8

    z1t1_BTAC = z1t1_BTAC/Nz1t1_BTAC
    z2t1_BTAC = z2t1_BTAC/Nz2t1_BTAC
    zt1_BTAC = zt1_BTAC/Nzt1_BTAC
    ztn_BTAC = ztn_BTAC/Nztn_BTAC

    Dc_FOL8 = abs(z1t1_FOL8 - z2t1_FOL8)
    Dc_BTAC = abs(z1t1_BTAC - z2t1_BTAC)

    Lz_FOL8 = 0.5*(abs(zt1_FOL8(1,:,:) - ztn_FOL8(1,:,:))+abs(zt1_FOL8(2,:,:) - ztn_FOL8(2,:,:)))
    Lz_BTAC = 0.5*(abs(zt1_BTAC(1,:,:) - ztn_BTAC(1,:,:))+abs(zt1_BTAC(2,:,:) - ztn_BTAC(2,:,:)))
      
    do ix=-nbw(1),nbw(1)
      do iy=-nbw(2),nbw(2)
          
          !------------------------------------

          if((.not. ieee_is_nan(Lz_FOL8(ix,iy))) .and. (.not. ieee_is_nan(Dc_FOL8(ix,iy)))) then
            k = nint(((2.d0*Lz_FOL8(ix,iy)-Dc_FOL8(ix,iy))/Lz_FOL8(ix,iy))/ddove)
            !if(k<0)  write(*,*) "FOL8",k,Lz_FOL8(ix,iy),Dc_FOL8(ix,iy),Nz1t1_FOL8(ix,iy),&
            !Nz2t1_FOL8(ix,iy),Nzt1_FOL8(:,ix,iy),Nztn_FOL8(:,ix,iy)
            if(k>=-nb .and. k<=nb) Dove_FOL8_dist(k) = Dove_FOL8_dist(k) + 1.d0
            k = nint(Lz_FOL8(ix,iy)/dlz)
            if(k>=-nb .and. k<=nb) Lz_FOL8_dist(k) = Lz_FOL8_dist(k) + 1.d0
            k = nint(Dc_FOL8(ix,iy)/ddc)
            if(k>=-nb .and. k<=nb) Dc_FOL8_dist(k) = Dc_FOL8_dist(k) + 1.d0
          endif
          
          if((.not. ieee_is_nan(Lz_BTAC(ix,iy))) .and. (.not. ieee_is_nan(Dc_BTAC(ix,iy)))) then
            k = nint(((2.d0*Lz_BTAC(ix,iy)-Dc_BTAC(ix,iy))/Lz_BTAC(ix,iy))/ddove)
            !if(k<0) write(*,*) "BTAC",k,Lz_BTAC(ix,iy),Dc_BTAC(ix,iy),Nz1t1_BTAC(ix,iy),&
            !Nz2t1_BTAC(ix,iy),Nzt1_BTAC(:,ix,iy),Nztn_BTAC(:,ix,iy)
            if(k>=-nb .and. k<=nb) Dove_BTAC_dist(k) = Dove_BTAC_dist(k) + 1.d0
            k = nint(Lz_BTAC(ix,iy)/dlz)
            if(k>=-nb .and. k<=nb) Lz_BTAC_dist(k) = Lz_BTAC_dist(k) + 1.d0
            k = nint(Dc_BTAC(ix,iy)/ddc)
            if(k>=-nb .and. k<=nb) Dc_BTAC_dist(k) = Dc_BTAC_dist(k) + 1.d0
          endif

          !------------------------------------
        
      enddo
    enddo

    !if(step_num > 10) exit
    cnfg = get_snap(tFilter)
enddo

write(*,'(2a)') "log: <<< end >>>"

call close_trajectory()


Dove_FOL8_dist = Dove_FOL8_dist / (sum(Dove_FOL8_dist)*ddove)
Dove_BTAC_dist = Dove_BTAC_dist / (sum(Dove_BTAC_dist)*ddove)

Dc_FOL8_dist = Dc_FOL8_dist / (sum(Dc_FOL8_dist)*ddc)
Dc_BTAC_dist = Dc_BTAC_dist / (sum(Dc_BTAC_dist)*ddc)

Lz_FOL8_dist = Lz_FOL8_dist / (sum(Lz_FOL8_dist)*dlz)
Lz_BTAC_dist = Lz_BTAC_dist / (sum(Lz_BTAC_dist)*dlz)

inquire(file=trim(folName(charus)),exist=existuje)
  if(.not. existuje) then
    !stop "neexistuje"
    call execute_command_line("mkdir "//trim(folName(charus)))
  endif
  
  

  
!///////////////////////////
open(1,file=trim(folName(charus))//"/Doverlap_dist.data")
write(1,*) "overlap FOL BTAC"
do i=-nb,nb
  write(1,*) real(i)*ddove,Dove_FOL8_dist(i),Dove_BTAC_dist(i)
enddo
close(1)

open(1,file=trim(folName(charus))//"/Dc_dist.data")
write(1,*) "Dc FOL BTAC"
do i=-nb,nb
  write(1,*) real(i)*ddc,Dc_FOL8_dist(i),Dc_BTAC_dist(i)
enddo
close(1)

open(1,file=trim(folName(charus))//"/Lz_dist.data")
write(1,*) "Lz FOL BTAC"
do i=-nb,nb
  write(1,*) real(i)*dlz,Lz_FOL8_dist(i),Lz_BTAC_dist(i)
enddo
close(1)



call gnupl%init()

call gnupl%set_terminal("pngcairo",&
'enhanced size 800,600 font "Helvetica, 12" linewidth 2')

call gnupl%set_dir("result"//trim(charus))
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


contains
!----------------------------------------------------------
character(len=50) function folName(charus) result(foldName)

character(len=*),intent(in) :: charus

foldName = "result"//trim(charus)

end function
end program

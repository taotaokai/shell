! Haskell matrix of P-SV motion in plane-layered isotropic elastic medium 
!
!--THEORY
!
!   A*(dV/dt) = dV/dz + F
!   , where V = [Vr, Vz, Tr, Tz] is the velocity-traction vector radial(r), vertical-down(z)
!
!   In frequency domain: 1i*w*A*V = dV/dz + F
!
!   The Jordan decomposition gives A = M * J * Minv
!
!   M: mode matrix (4 by 4 matrix),
!       M(:,1|3): down-going P|SV, M(:,2|4): up-going P|SV;
!
!   Minv: inverse mode (projection) matrix (4 by 4 matrix); 
!       Minv(1|3,:): down-going P|SV, Minv(2|4,:): up-going P|SV.
!
!   J: diagonal matrix of vertical slownesses, J = diag([-qp,qp,-qs,qs]);
!
!   In one homogeneous layer: V(z) = M * exp(1i*w*J*z) * Minv * V(0)
!
!   The Haskell matrix is defined as H(z;0) = M * exp(1i*w*J*z) * Minv
!
!--NOTE 
!
! 1) for near critical rayp (qs or qp = 0), 1/0 will occur in Minv.
!	However this will be circumvented by combining the phase delay term
!	in the propagation matrix which will lead to the terms like: sin(qs*tau)/qs
!	and no singularity occurs at qs=0;
!
! 2) vertical direction points downward to the earth's center;
!
! 3) Polarization of down- and up-going waves
! 
!                     Down-going |  Up-going
!  (tangential) T.-------------------------------> R(radial)
!                |      SV       |       P
!                |     /         |      /
!                | SH *          |  SH * 
!                |     \         |      \
!                |       P       |       SV
!                V   
!                Z(down)  
!
! 4) Inverse Fourier transformation is defined as: F(t) = Int[F(w)*exp(1i*w*t),{w,-Inf,Inf}]
!
!--HISTORY
!
! [2012-05-27] created
! [2012-06-09] modified: change the polarization
! [2012-11-17] modified: change the polarization
! [2013-03-16] export to fortran 
! [2017-06-17] add RF_PSV_ISO

!
!-------------------------------------------------------------
!

subroutine RF_PSV_ISO(nz,z,vp,vs,rho,nw,w,np,p,itype, Vr,Vu)
!-- P-SV receiver site response (free surface, i.e. zero traction )
!
!-- INPUT
!
!   nz: number of layers (including sub-sediment layer)
!   z(nz): array of layer thickness
!   vp(nz): array of P-wave velocity
!   vs(nz): array of S-wave velocity
!   rho(nz): array of density
!
!   w(nw): array of angular frequencies
!
!   np: number of ray parameters
!   p(np): array of ray parameters
!
!   itype: incident wave type, 0=P, 1=SV
!
!-- OUTPUT
!
!   Vr(nw,np), Vu(nw,np): radial/vertical(upward positive) particle velocities at the surface
!
implicit none

complex*16, parameter :: XJ = (0,1)
complex*16, parameter :: XZERO = (0,0)
complex*16, parameter :: XONE = (1,0)

integer, intent(in) :: nz, nw, np
!f2py integer intent(hide),depend(z) :: nz = len(z)
!f2py integer intent(hide),depend(w) :: nw = len(w)
!f2py integer intent(hide),depend(p) :: np = len(p)
real*8, dimension(nz), intent(in) :: z, vp, vs, rho
complex*16, dimension(nw), intent(in) :: w
real*8, dimension(np), intent(in) :: p
integer, intent(in) :: itype

complex*16, dimension(nw,np), intent(out) :: Vr
complex*16, dimension(nw,np), intent(out) :: Vu

! local variables
complex*16, dimension(4,4,nw) :: H
real*8, dimension(4,4) :: M, Minv, Q
complex*16, dimension(4,4) :: MinvH
complex*16 :: a,b,c,d,det
integer :: iw, ip

!------ loop each event
do ip = 1,np

  ! Mode matrix in the last layer
  call Mode_PSV_ISO(vp(nz),vs(nz),rho(nz),p(ip), M,Minv,Q)

  ! Haskell matrix
  call Haskell_PSV_ISO(nz,z,vp,vs,rho,nw,w,p(ip), H)

  ! apply boundary conditions
  ! Minv*H * [Vr, Vz, 0, 0] = [?, 1, ?, 0] for P-wave incidence
  ! Minv*H * [Vr, Vz, 0, 0] = [?, 0, ?, 1] for SV-wave incidence
  do iw = 1,nw

    MinvH = matmul(Minv, H(:,:,iw))

    ! sub-matrix of MinvH
    a = MinvH(2,1)
    b = MinvH(2,2)
    c = MinvH(4,1)
    d = MinvH(4,2)
    det = a*d-b*c ! determinant

    if (abs(det) < 1e-8) then
      write(*,*) "[WARN] very small determinant of Minv*H: ", abs(det)
    endif

    ! inverse of sub-matrix of MinvH
    if (itype == 0) then
      Vr(iw,ip) = d/det
      Vu(iw,ip) = -1.0*(-c/det) ! change sign such that up direction is positive
    else if (itype == 1) then
      Vr(iw,ip) = -b/det
      Vu(iw,ip) = -1.0*(a/det) ! change sign such that up direction is positive 
    else
      write(*,*) "[ERROR] Wrong input of incident wave type (itype=0/1 for P/SV): ", itype
      stop
    endif

  end do

  ! phase shift

end do ! loop each event

end subroutine

!
!-------------------------------------------------------------
!

subroutine DC_PSV_ISO(nz,z,vp,vs,rho,nw,w,np,p,V0, V1)
!-- P-SV surface wavefield downward continuation and decomposition
!
!--INPUT
!
!   nz: number of layers (including sub-sediment layer)
!   z(nz): array of layer thickness
!   vp(nz): array of P-wave velocity
!   vs(nz): array of S-wave velocity
!   rho(nz): array of density
!
!   nw: number of frequency opoints
!   w(nw): angular frequency array 
!
!   np: number of ray parameters
!   p(np): array of ray parameters
!
!   V0(4,nw,np): surface velocity-stree vector, dim-1: vr,vz(downward positive),tr,tz,
!                dim-2: frequency samples, dim-3: ray parameters
!
!--OUTPUT
!
!   V1(4,nw,np): mode vector of downward continuated wavefield, dim-1: down-P,up-P,down-S,up-S
!
!--NOTE
!
!   1. The sign for vz is downward positive.
!   2. The wavefield is downward extrapolated to the bottom of the last layer
!   and decomposed.

implicit none

complex*16, parameter :: XJ = (0,1)
complex*16, parameter :: XZERO = (0,0)
complex*16, parameter :: XONE = (1,0)

integer, intent(in) :: nz, nw, np
!f2py integer intent(hide),depend(p) :: np = len(p)
!f2py integer intent(hide),depend(w) :: nw = len(w)
!f2py integer intent(hide),depend(z) :: nz = len(z)
real*8, dimension(nz), intent(in) :: z, vp, vs, rho
complex*16, dimension(nw), intent(in) :: w
real*8, dimension(np), intent(in) :: p
complex*16, dimension(4,nw,np), intent(in) :: V0

complex*16, dimension(4,nw,np), intent(out) :: V1

! local variables
complex*16, dimension(4,4,nw) :: H
real*8, dimension(4,4) :: M, Minv, Q
integer :: iw, ip

do ip = 1,np ! loop each event

  ! Mode matrix in the last layer
  call Mode_PSV_ISO(vp(nz),vs(nz),rho(nz),p(ip), M,Minv,Q)

  ! Haskell matrix
  call Haskell_PSV_ISO(nz,z(1:nz),vp(1:nz),vs(1:nz),rho(1:nz), nw,w, p(ip), H)
  do iw = 1,nw
    V1(:,iw,ip) = matmul(Minv, matmul(H(:,:,iw), V0(:,iw,ip)) )
  end do

  !if (nz > 1) then
  !  call Haskell_PSV_ISO(nz-1,z(1:nz-1),vp(1:nz-1),vs(1:nz-1),rho(1:nz-1), nw,w, p(ip), H)
  !  do iw = 1,nw
  !    V1(:,iw,ip) = matmul(Minv, matmul(H(:,:,iw), V0(:,iw,ip)) )
  !  end do
  !else
  !  do iw = 1,nw
  !    V1(:,iw,ip) = matmul(Minv, V0(:,iw,ip))
  !  end do
  !end if

end do ! loop each event

end subroutine DC_PSV_ISO

!
!--------------------------------------------------------
!

subroutine Haskell_PSV_ISO(nz,z,vp,vs,rho,nw,w,p, H)
! P-SV Haskell matrix for layered isotropic elastic media
!
!--INPUT
!
!   real*8 :: z(nz): layer thicknesses
!   real*8 :: vp(nz), vs(nz), rho(nz): Vp(km/s), Vs(km/s) and density(kg/m^3) for each layer
!
!   complex*16 :: w(nw): angular frequencies
!   real*8 :: p: ray parameter (s/km)
!
!--OUTPUT
!
!   complex*16 :: H(4,4,nw): Haskell matrix of P-SV motion for each frequency sample
!
!--NOTE
!
! inverse Fourier transformation is defined as: F(t) = Int[H(w)*exp(1i*w*t),{w,-Inf,Inf}]
!
!   In one homogeneous layer:
!
!   1i*w*A*v = dv/dz, A = M*Q*Minv, v(z) = M*exp(1i*w*Q*z)*Minv*v(0)
!   H(z;0) = M*exp(1i*w*Q*z)*Minv

implicit none

complex*16, parameter :: XJ = (0,1) 
complex*16, parameter :: XZERO = (0,0) 
complex*16, parameter :: XONE = (1,0) 

! inputs
integer, intent(in) :: nz, nw
!f2py integer intent(hide),depend(z) :: nz = len(z)
!f2py integer intent(hide),depend(w) :: nw = len(w)
real*8, dimension(nz), intent(in) :: z, vp, vs, rho
complex*16, dimension(nw), intent(in) :: w
real*8, intent(in) :: p

! output
complex*16, dimension(4,4,nw), intent(out) :: H

! local variables
complex*16, dimension(4,4) :: Jm
real*8, dimension(4,4) :: M, Minv
real*8, dimension(4) :: Q
integer :: i, iz, iw

! initialize H to identity matrix for each frequency sample
H = XZERO
do i = 1,4
  H(i,i,:) = XONE
enddo

! cumulative product of Haskell matrices of each layer
do iz = 1,nz

  call Mode_PSV_ISO(vp(iz),vs(iz),rho(iz),p, M,Minv,Q)

  do iw = 1,nw
    ! phase matrix of up-/down-going P and SV waves from the top to the bottom in each layer
    Jm = XZERO
    do i = 1,4
      Jm(i,i) = exp(XJ*w(iw)*Q(i)*z(iz))
    enddo
    H(:,:,iw) = matmul( matmul(M, matmul(Jm, Minv)), H(:,:,iw))
  end do

end do

end subroutine

!
!--------------------------------------------------------
!

subroutine Mode_PSV_ISO(vp,vs,rho,p, M,Minv,Q)
!
!--INPUT
!   real*8 :: vp,vs,rho: Vp, Vs(km/s) and density(g/cm^3)
!   real*8 :: p: ray parameter(s/km)
!
!--OUTPUT
!   real*8 :: M(4,4), Minv(4,4), mode and inverse mode matrix for P-SV motion
!   real*8 :: Q(4), vertical slownesses

implicit none

real*8, intent(in) :: vp, vs, rho, p

real*8, dimension(4,4), intent(out) :: M, Minv
real*8, dimension(4), intent(out) :: Q

real*8 :: p2, vs2, mu, qp, qs, qss

p2 = p**2
vs2 = vs**2
mu = rho*vs2
qp = sqrt(1/vp**2-p2)
qs = sqrt(1/vs2-p2)
qss = 1/vs2-2*p2

! Mode matrix
!               down P              up P      down SV            up SV
M(1,:) = ( (/         p*vp,          p*vp,         qs*vs,         qs*vs /) ) ! Vr
M(2,:) = ( (/        qp*vp,        -qp*vp,         -p*vs,          p*vs /) ) ! Vz
M(3,:) = ( (/-2*mu*p*qp*vp,  2*mu*p*qp*vp,    -mu*qss*vs,     mu*qss*vs /) ) ! Tr
M(4,:) = ( (/   -mu*qss*vp,    -mu*qss*vp,  2*mu*p*qs*vs,  2*mu*p*qs*vs /) ) ! Tz

! inverse of the Mode matrix
!                  Vr                Vz                 Tr            Tz         
Minv(1,:) = ( (/        p*mu/vp,   mu*qss/qp/2/vp,  -p/qp/2/vp,    -0.5/vp /) ) ! down P
Minv(2,:) = ( (/        p*mu/vp,  -mu*qss/qp/2/vp,   p/qp/2/vp,    -0.5/vp /) ) ! up P
Minv(3,:) = ( (/ mu*qss/qs/2/vs,         -p*mu/vs,     -0.5/vs,  p/qs/2/vs /) ) ! down SV
Minv(4,:) = ( (/ mu*qss/qs/2/vs,          p*mu/vs,      0.5/vs,  p/qs/2/vs /) ) ! up SV
Minv = Minv/rho

! vertical slowness 
!     down  up  down up
Q = ((/-qp, qp, -qs, qs/))

end subroutine

!
!------------------------------------------------
!

program test_Haskell

implicit none

integer, parameter :: nz = 2
double precision, dimension(nz) :: vp, vs, rho, z
DATA vp/4.d0, 6.0d0/ vs/2.d0, 4.0d0/ rho/2.7d0,3.0d0/ z/1.d0,4.d0/

integer, parameter :: nw = 2
complex*16, dimension(nw) :: w

double precision, parameter :: p = 0.05 

!double precision, dimension(4,4) :: M, Minv, Q 

complex*16, dimension(4,4,nw) :: H

integer :: i

w(1) = (1.0, 0.0)
w(2) = (2.0, 0.0)

call Haskell_PSV_ISO(nz,z,vp,vs,rho, nw,w, p, H)
  
do i = 1, nw
  write(*,*) "iw = ", i
  write(*,*) H(:,:,i)
enddo

end program test_Haskell

  program hartree
  implicit none
!=================================================================================================================================
! Contains the following procedures:
!
!                                   -------------------------
!                                   |   F U N C T I O N S   |
!                                   -------------------------
!    i. Fact(n) -- outputs the factorial for a given integer value
!
!   ii. Gam(n,l) -- outputs a gamma function, gamma(n + l + 3.0/2.0) for quantum numbers n and l
!
!  iii. Laguerre(l, max_n, b, r) -- outputs the Hermite polynomial associated with wavefunction solution for a harmonic oscillator
!
!==================================================================================================================================
!                                   -------------------------
!                                   | S U B R O U T I N E S |
!                                   -------------------------
!    i. Subroutine setpsi(l, n, max_n, nr, nsq, dr, psi)--sets harmonic oscillator states
!
!   ii. Subroutine density(l, n, max_n, nr, nsq, dr, Lag, psi, rho, rho_tot)-- uses psi(nr, nsq, 2) to calculate the density
!       rho(nr, nsq, 2) for neutrons and protons individual as well as total density rho_tot(nr,nsq).  These arrays are integrated
!
!  iii. Subroutine damp(r, dr, l, T, nr)-- Damping operator applied to harmonic oscillator states when calculating psitmp(nr, nsq, 2)
!
!   vi. Subroutine Coulomb(nr, nsq, dr, phi_c, rho)-- solves the Poisson equation to find the Coulomb
!       contribution to the potential.
!
!    v. Subroutine yukawa(nr, nsq, pi, phi_y, rho_tot, dr)-- solves the Helmhotz equation to find
!       contribution from the Yukawa potential (both proton and neutron).
!
!   iv. Subroutine BKN(phi_c, phi_y, nr, nsq, dr, rho, rho_tot, UBKN)-- uses density calculated
!       from the wavefunction, and both Coulomb and Yukawa to calculate total potential.
!
!   vi. Subroutine spet(nr, nsq, l, n, dr, psi, UBKN, spe, tpe)-- finds eigenvalue spe(ns, iq) and kinetic
!       energy
!
!  vii. Subroutine grstep(T, psi, psitmp, spe, l, dr, nr, nsq, UBKN)-- calculates new wavefunction psitmp using damping operator
!       and single particle energies
!
! viii. Subroutine SHF(l, nr, nsq, dr, rho, rho_tot, spe, tpe, EHF)
!---------------------------------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------------------
!   P a r a m e t e r s:
!----------------------------------------------------------------------------------------------------------------------------------------------
      integer, parameter                   :: wp = kind(1.0d0)
      real(wp), parameter                  :: pi = acos(-1.0d0)
      integer, parameter                   :: nr = 200               ! mesh points, number of proton or neutron states
!----------------------------------------------------------------------------------------------------------------------------------------------
!   V a r i a b l e s:
!----------------------------------------------------------------------------------------------------------------------------------------------
      integer                              :: i, j, is, iq, il       ! radial index, proton/neutron states, type fermion
      integer                              :: max_n, m               ! principal quantum #, major #, angular momentum
      integer                              :: blank
      integer                              :: A, nsq                 ! atomic mass, and maximun number of n,l values
      integer                              :: nsl                    ! maximum number of groups of l values
      integer                              :: lcount(2)              ! counts groups lstates
      integer                              :: max_shell(2)
!
      real(wp)                             :: sum1, sum2             ! sum1 to normalize chi, sum2 to calculate # protons and neutrons
      real(wp)                             :: rmax, dr               ! rmax in femtometers, rstep
      real(wp)                             :: sum_n, sum_p, sum_t    ! total # neutrons, protons (and both combined)
!----------------------------------------------------------------------------------------------------------------------------------------------
!   A r r a y s:
!----------------------------------------------------------------------------------------------------------------------------------------------
      real(wp), dimension(nr)              :: r                      ! radial distance array
      real(wp), allocatable                :: psi(:,:,:)             ! factor functions for wave function
      real(wp), allocatable                :: psitmp(:,:,:)          ! temporary wavefunction in gradient iteration subroutine
      integer, allocatable                 :: n(:,:), l(:,:)         ! quantum numbers array
!
      real(wp), dimension(nr,2)            :: rho                    ! rho(1:nr,1) = neutron density, rho(1:nr,2) proton density
      real(wp), dimension(nr)              :: rho_tot                ! total density
!
      real(wp), allocatable                :: T(:,:,:,:)             ! Damping operator
      real(wp)                             :: Q(nr,nr)
!
      real(wp), dimension(nr)              :: phi_c                  ! Coulomb potential
!
      real(wp), dimension(nr,2)            :: phi_y                  ! Yukawa potential
!
      real(wp), dimension(nr)              :: Un                     ! BKN potential for neutron
      real(wp), dimension(nr)              :: Up                     ! BKN potential for proton
      real(wp), dimension(nr,2)            :: UBKN                   ! BKN for both neutron and proton
!
      real(wp), allocatable                :: spe(:,:)               ! single particle total energy
      real(wp), allocatable                :: tpe(:,:)               ! single particle kinetic energy
      real(wp), allocatable                :: occ(:,:)
      integer,  allocatable                :: lst(:,:)
!
      integer, dimension(2)                :: num                    ! num(1) = no. neutrons, num(2) = no. protons
      integer, dimension(2)                :: num_states             ! num_states(1) = no. available neutron states, num_states(2) = no. available proton states
      integer, dimension(2)                :: tot_lvl                ! tot_lvl(1) = no. occupied neutron levels, tot_lvl(2) = no. occupied proton levels
!
      real(wp)                             :: EHF
      real(wp)                             :: tolerance
!
      character (len=16), dimension(2)     :: file_name = [ 'n_l_neutron.txt', 'n_l_proton.txt ']   ! Files contain N, l, n values
      character (len=16), dimension(2)     :: lgroups   = [ 'lg_neutrons.txt', 'lg_protons.txt ']   ! File that groups l values in ascending order
      character (len=12), dimension(2)     :: FERM(2)   = ['FOR NEUTRONS', 'FOR PROTONS ']
!
!----------------------------------------------------------------------------------------------------------------------------
! Initialize:
!----------------------------------------------------------------------------------------------------------------------------
!
! set parameters for wave function expressions and rho runction--fixed for now, but eventually we'll program a loop:
! -integers
    num_states = 0
!
! -floats
    rmax = 20.0_wp                 ! max radial distance in femtometers
    dr   = rmax/nr
    r    = 0.0_wp
!
! Prompt user for neutron, proton number:
      print *, " Enter no. neutron and protons:"
      read *, num(1), num(2)
      print *, " "
      print *, "Atomic mass = ", sum(num)
      print *, " "
!
    A      = sum(num)
!
    r = [(i * dr, i = 1, nr)]
!
    call find_states(num, num_states, tot_lvl, nsq, file_name, max_shell)
!
    call occupation(num, num_states, nsq, tot_lvl, occ, l, n, file_name, lgroups, lst, nsl, max_shell)
!
    call setpsi(r, n, l, nr, nsq, dr, psi, tot_lvl, A)
!
    call damp(T, Q, r, nr, nsq, dr, l, tot_lvl)
!
     do i = 1, 130
!
    call density(r, l, n, nsq, tot_lvl, nr, dr, psi, rho, rho_tot, sum_n, sum_p, sum_t, occ)
!
    call Coulomb(r, nr, dr, phi_c, rho, sum_p)
!
    call yukawa(r, nr, pi, phi_y, rho, dr)
!
    call BKN(r, phi_c, phi_y, nr, dr, rho, rho_tot, UBKN)
!
    call spet(r, nr, nsq, l, dr, tot_lvl, psi, UBKN, spe, tpe)
!
    call grstep(r, l, tot_lvl, T, psi, psitmp, spe, dr, nr, nsq, UBKN, lst, nsl, il, m, blank)
!
    call SHF(r, l, tot_lvl,nr, nsq, dr, rho, rho_tot, spe, tpe, EHF)
!
    end do
!
!
        print *, "========================================================================="
        print *, "EHF = " , Ehf, "MeV"
        print *, "========================================================================="
        do iq = 1, 2
          print *, '   ', FERM(iq), ':'
          do is = 1, tot_lvl(iq)
            print *, "l =", l(is,iq),"n =", n(is,iq), "      spe =", spe(is,iq)
          end do
        end do
!
  contains
!===================================================================================================================================================
!===================================================================================================================================================
	subroutine find_states(num, num_states, tot_lvl, nsq, file_name, max_shell)
!
! Local variables:
	integer, parameter   :: wp = kind(1.0d0)
	integer              :: sum1, iq
	integer              :: il, in, is
	integer              :: Big_N, lil_n
	character (len=16), dimension(2), intent(inout)   :: file_name
!
! Global variables:
    integer, intent(out)  :: nsq
	integer, intent(in)   :: num(2)
	integer, intent(out)  :: num_states(2)
	integer, intent(out)  :: tot_lvl(2)
	integer, intent(out)  :: max_shell(2)
!
	tot_lvl = 0
	num_states = 0
!
! Determine the number of states for neutrons/protons:
!
	do iq = 1, 2
!
		open (unit = 11, file = file_name(iq))
! 		print *, '   ', FERM(iq), ':'
!
		do Big_N = 0, 20
!
		if (mod(Big_N , 2) .eq. 0) then
			lil_n = Big_N - Big_N/2
		else if (mod(Big_N , 2) .eq. 1) then
			lil_n = Big_N - (Big_N + 1)/2
		end if
!
		if (num_states(iq) .ge. num(iq)) exit
!
! 		print *, "N = ", Big_N
! 		print *, "========================================="
!
    do in = lil_n, 0, -1
!
        il = Big_N - 2 * in
!
!             print *, "n = ", in, " l = ", l
            tot_lvl(iq) = tot_lvl(iq) + 1
!
            num_states(iq) = num_states(iq) + 2 * ( 2 * il + 1)
!
            write (11, fmt = 100) Big_N, il, in
            100 format (2x, i2, 2x, i2, 2x, i2)
!
			end do
! 			print *, "          "
		end do
		max_shell(iq) = Big_N - 1
		close(11)
	end do
!
    nsq = maxval(tot_lvl)
!
	end subroutine find_states
!=========================================================================================================================================
!=========================================================================================================================================
    subroutine occupation(num, num_states, nsq, tot_lvl, occ, l, n, file_name, lgroups, lst, nsl, max_shell)
    implicit none
!
!  Local variables:
    integer, parameter                             :: wp = kind(1.0d0)
    integer                                        :: i, j, is, iq, sum_s
!
!  Local array:
    real(wp)                                       :: frac(2)
    integer                                        :: lcount(2)
    integer                                        :: A(2)
!
!  Global variables:
    integer, intent(in)                            :: nsq
    integer, intent(out)                           :: nsl
!
!  Global array:
    integer, intent(in)                            :: num(2)
    integer, intent(in)                            :: num_states(2)
    integer, intent(in)                            :: tot_lvl(2)
    integer, intent(in)                            :: max_shell(2)
    integer, allocatable, intent(out)              :: l(:,:), n(:,:)
    integer, allocatable, intent(out)              :: lst(:,:)
    real(wp), allocatable, intent(out)             :: occ(:,:)
	character(len=16), intent(in)                  :: file_name(2)
    character(len=16), dimension(2), intent(inout) :: lgroups
!
	allocate (occ(nsq,2))
	allocate (l(nsq,2))
	allocate (n(nsq,2))
!
    occ    = 0.0_wp
	l      = 0
	n      = 0
	lcount = 0
	A      = 0

!
	do iq = 1, 2
      open (unit = 11, file = file_name(iq))
		do is = 1, tot_lvl(iq)
			read (11, fmt = 110) l(is,iq), n(is,iq)
			110 format (6x, i2, 2x, i2)
			occ(is, iq) = 2.0_wp * (2.0_wp * l(is,iq) + 1.0_wp)
		end do
      close(11)
	end do
!
    do iq = 1, 2
      open (unit = 12, file = lgroups(iq))
        do m = 0, 10
            do is = 1, tot_lvl(iq)
               i = l(is,iq)
               j = n(is,iq)
            if ( i .eq. m ) write (12, fmt = 120) i, j, is
               120 format (2x, i2, 2x, i2, 2x, i2)
            end do
            if (count (l(1:nsq,iq) == m) .gt. 0) lcount(iq) = lcount(iq) + 1
        end do
      close(12)
	end do
!
	nsl = maxval(lcount) - 1
!
	allocate ( lst(0:nsl,2) )
	lst = 0
!
	do iq = 1, 2
        do is = 0, (lcount(iq) - 1)
          lst(is,iq) = count( l(1:tot_lvl(iq),iq) == is )
      end do
    end do

    do iq = 1, 2
      sum_s = 0
      do i  = 0, max_shell(iq)
        A(iq) = sum_s
        frac(iq) = (num(iq) - A(iq)) / ( sum(occ(:,iq)) - A(iq))
        sum_s =  sum_s + (i + 1) * (i + 2)
        if (sum_s .gt. num(iq)) exit
      end do
    end do

    do iq = 1, 2
      is = 1
      do while (sum(occ(1:is,iq)) .le. A(iq))
      is = is + 1
      end do
      occ(is:nsq,iq) = occ(is:nsq,iq) * frac(iq)
    end do
!
    end subroutine occupation
!=========================================================================================================================================================
!=========================================================================================================================================================
  Subroutine setpsi(r, n, l, nr, nsq, dr, psi, tot_lvl, A)
  implicit none

    integer, parameter   :: wp = kind(1.0d0)

! Local variables:
    integer  :: i, is, iq                                     ! radial index, proton/neutron states index, type fermion index
    real(wp) :: b                                             ! sqrt of h-bar/(mass * omega_w)
    real(wp) :: sum1, sum2                                    ! sum1 to normalize chi,

! Local arrays:
    real(wp), allocatable                                :: Lag(:,:,:)
    real(wp), allocatable                                :: summ(:,:)

! Global variables:
    integer, intent(in)  :: nr, nsq                           ! mesh points, number of proton or neutron states
    real(wp), intent(in) :: dr                                ! rstep
    integer, intent(in)  :: A                                 ! atomic mass

! Global arrays:
    real(wp), dimension(nr), intent(in)                  :: r
    integer,  dimension(2), intent(in)                   :: tot_lvl
    integer,  intent(in)                                 :: l(:,:), n(:,:)
    real(wp), allocatable, intent(out)                   :: psi(:,:,:)



!  Initialize:
    b   = 1.010_wp * (A**(1.0_wp/6.0_wp))           ! value for Oxygen 16
    Lag = 0.0_wp

!--------------------------------------------------------------------------------------------------------------------------
!                     S e t  p s i (i, is, iq)
!--------------------------------------------------------------------------------------------------------------------------
! For a harmonic oscillator function, the original expression for psi is:
!
! psi(i, is, iq) = (b**(-3/2)) * (((2.0_wp*fact(n))/(gamma(n + l + 3.0_wp/2.0_wp))**(1.0_wp/2.0_wp)))&
!   * ((r(i)/b)**l) * exp((-1.0_wp/2.0_wp)*((r(i)/b)**2))*Lag(i).
!--------------------------------------------------------------------------------------------------------------------------

! Set sizes of allocatable arrays:
  allocate ( psi(0:nr+1,nsq,2) )
  allocate ( Lag(nr,nsq,2) )
  allocate ( summ(nsq,2) )

! Fill Laguerre array:
      do iq = 1, 2
        do is = 1, tot_lvl(iq)
            if (n(is,iq) .eq. 0)then
              do i = 1, nr
                Lag(i,is,iq) = 1.0_wp
              end do
           else
              do i = 1, nr
                Lag(i,is,iq) = Laguerre(l(is,iq), n(is,iq), b, r(i))
              end do
           end if
          end do
        end do

! Fill psi(i, is, iq):
      do iq = 1, 2
        do is = 1, tot_lvl(iq)
            do i = 1, nr
              psi(i,is,iq) = r(i) * ( (b**(-3.0_wp/2.0_wp)) * (2.0_wp * fact(n(is,iq)) / gam(n(is,iq),l(is,iq)))**(1.0_wp/2.0_wp)&
                              * (r(i)/b)**l(is,iq) * exp( (-1.0_wp/2.0_wp ) * (r(i)/b)**2) * Lag(i,is,iq) )
            end do
        end do
      end do

! Normalize psi:
      do iq = 1, 2
        do is = 1, tot_lvl(iq)
          summ(is,iq) = sum(psi(1:nr,is,iq)**2) * dr
          psi(:,is,iq) = psi(:,is,iq)/sqrt(summ(is,iq))
        end do
      end do

  end subroutine setpsi

  pure FUNCTION Fact(n)

IMPLICIT NONE
INTEGER :: Fact, i
INTEGER, INTENT(IN) :: n

IF (n == 0) THEN
Fact = 1
ELSE
Fact = 1
Do i = 1, n
  fact = i * fact
End Do
END IF

END FUNCTION Fact
!===========================================================================================================================================
  pure function gam(n, l)

    integer, parameter  :: wp = kind(1.0d0)
    integer, intent(in) :: n, l
    real(wp)            :: gam

    gam = gamma(n + l + 3.0_wp/2.0_wp)

    end function gam

!===========================================================================================================================================
  pure elemental function Laguerre(l, max_n, b, r)

!----------------------------------------------------------------------------------------------------------------------------------------------------------------------
!       alpha : alpha in the recurrence relation = l + 1/2
!       max_n : order of the polynomial desired
!       x     : the argument of the function
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------
        implicit none
        integer, parameter :: wp = kind(1.0d0)
!
        real(wp) :: alpha, Laguerre, LNM1, LN, LNP1
        integer  :: order, n, i
        integer, parameter   :: nr = 50
        integer, intent(in)  :: l, max_n
        real(wp), intent(in) :: b
        real(wp), intent(in) :: r



        order = max_n
        alpha =  l + 1.0_wp/2.0_wp


        LNM1 = 1.0_wp
        LN   = -(r/b)**2 + alpha + 1.0_wp
        do n = 1, (max_n - 1)
           LNP1 = (LN*(2.0_wp*n + alpha + 1.0_wp - (r/b)**2) - LNM1*(n + alpha))/(n + 1.0_wp)
           LNM1 = LN
           LN = LNP1
        end do
        Laguerre = LN
!
        end function Laguerre
!=====================================================================================================================================================================
!=====================================================================================================================================================================
  Subroutine damp(T, Q, r, nr, nsq, dr, l, tot_lvl)
  implicit none

    integer, parameter        :: wp = kind(1.0d0)

! Global variables:
    integer, intent(in)       :: nr, nsq
    real(wp), intent(in)      :: dr
!
! Global arrays:
    real(wp), intent(in)                   :: r(nr)
    integer, intent(in)                    :: l(nsq,2)
    integer, intent(in)                    :: tot_lvl(2)
    real(wp), allocatable, intent(inout)   :: T(:,:,:,:)

! Local variables:
    integer       :: i, j                 ! indexes for arrays
    real(wp)      :: beta                 ! average of hbar**2/(2*mass) for masses of proton and neutron
    real(wp)      :: E0                   ! E nought
    real(wp)      :: b, c                 ! diagonal elements of T

! Local arrays:
    real(wp)      :: a(nr,nsq,2)          ! diagonal elements of T
    real(wp), intent(inout)      :: Q(nr,nr)             ! matrix to test inversion of T
!
! Initialize:
    allocate ( T(nr,nr,nsq,2) )

    beta = 20.726_wp
    E0   = 5.0_wp
    a    = 0.0_wp
    T    = 0.0_wp
!
! Loop for diagonal elements a(i):
    do iq = 1, 2
      do is = 1, tot_lvl(iq)
        do i = 1, nr
          a(i, is, iq) =  (beta/E0) * ( (l(is,iq)*l(is,iq) + l(is,iq))/r(i)**2 + 2.0_wp/dr**2) + 1.0_wp
        end do
      end do
    end do
!
! Off-diagonal elements b, c
    b = -beta/(dr**2 * E0)
    c =  b
!
! fill T(i,j) matrix array:
    do iq = 1, 2
      do is = 1, tot_lvl(iq)
        do j = 1, nr
          do i = 1, nr
            if( i .eq. j) then
              T(i,j,is,iq) = a(i,is,iq)
            elseif( i .eq. j - 1) then
              T(i,j,is,iq) = b
            elseif( i .eq. j + 1) then
              T(i,j,is,iq) = c
            else
              T(i,j,is,iq) = 0.0_wp
            end if
          end do
        end do
      end do
    end do

!
  do iq = 1, 2
      do is = 1, tot_lvl(iq)
        Q = T(1:nr,1:nr,is,iq)
        call matinv(Q, nr)     ! inverts T(i,j) matrix
        T(1:nr,1:nr,is,iq) = Q
      end do
  end do
!
  end subroutine damp
!
! !===============================================================================================================================

      subroutine matinv(a, ndm)
      implicit none
      integer, parameter  ::  wp = kind(1.0d0)
! ! !-----------------------------------------------
! ! !   a r g u m e n t s
! ! !   n=ndm for square matrices
! ! !-----------------------------------------------
      integer , intent(in) :: ndm
      real(wp), intent(inout) :: a(ndm,ndm)
! ! !-----------------------------------------------
! ! !   l o c a l   v a r i a b l e s
! ! !-----------------------------------------------
      integer , dimension(1000)     :: index
      integer                       :: i, j, k, irow, icolum, i1, n
      real(wp), dimension(1000) :: pivot
      real(wp)                  :: amax, temp, swap, t
! ! !--------------------------------------------------------------------
! ! !        subroutine for inversion of a general-real matrix.
! ! !        modified by a.s. umar
! ! !
! ! !        a - on output contains the inverse matrix
! ! !        ndm - the maximum dimension of a in the calling routine
! ! !        n - the actual dimension used in calculations (n<=ndm)
! ! !---------------------------------------------------------------------
! ! !
! ! !        initialize pivot element array
! ! !
      n = ndm
      pivot(:n) = 0.0_wp
      index(:n) = 0
! ! !
! ! !        perform successive pivot operations (grand loop)
! ! !
      do i = 1, n
! ! !
! ! !        search for pivot element
! ! !
         amax = 0.0_wp
         do j = 1, n
            if (pivot(j) /= 0.0_wp) cycle
            do k = 1, n
               if (pivot(k) /= 0.0_wp) cycle
               temp = dabs(a(j,k))
               if (temp < amax) cycle
               irow = j
               icolum = k
               amax = temp
            end do
         end do
         index(i) = 4096*irow + icolum
         j = irow
         amax = a(j,icolum)
         pivot(icolum) = amax
!
! !        interchange rows to put pivot element on diagonal
! !
         if (irow /= icolum) then
            do k = 1, n
               swap = a(j,k)
               a(j,k) = a(icolum,k)
               a(icolum,k) = swap
            end do
         endif
! !
! !        divide pivot row by pivot element
! !
         k = icolum
         a(icolum,k) = 1.0_wp
         a(icolum,:n) = a(icolum,:n)/amax
! !
! !        reduce non-pivot rows
! !
         do j = 1, n
            if (j == icolum) cycle
            t = a(j,icolum)
            a(j,icolum) = 0.0_wp
            a(j,:n) = a(j,:n) - a(icolum,:n)*t
         end do
      end do
! !
! !     interchange columns after all pivot operations have been performed
! !
      do i = 1, n
         i1 = n + 1 - i
         k = index(i1)/4096
         icolum = index(i1) - 4096*k
         if (k == icolum) cycle
         do j = 1, n
            swap = a(j,k)
            a(j,k) = a(j,icolum)
            a(j,icolum) = swap
         end do
      end do
      return
      end subroutine matinv
!
!===============================================================================================================================================
!=====================================================================================================================================================================
 subroutine density(r, l, n, nsq, tot_lvl, nr, dr, psi, rho, rho_tot, sum_n, sum_p, sum_t, occ)
  implicit none
!
! Local variables:
    integer, parameter       :: wp = kind(1.0d0)
    real(wp), parameter      :: pi = acos(-1.d0)
    integer                  :: i, is, iq
!
! Global variables:
    integer, intent(in)      :: nr, nsq
    real(wp), intent(in)     :: dr
    real(wp), intent(out)    :: sum_n                ! to integrate rho_p (finds # of protons)
    real(wp), intent(out)    :: sum_p                ! to integrate rho_n (finds # of neutrons)
    real(wp), intent(out)    :: sum_t                ! to normalize rho_tot
!
! Local Arrays:
    real(wp)                 :: A(nr,2)
!
! Global Arrays:
    real(wp), intent(in), dimension(nr)               :: r
    real(wp), intent(in), dimension(nsq,2)            :: occ(nsq,2)
    integer, intent(in), dimension(2)                 :: tot_lvl
    integer, intent(in), dimension(:,:)               :: l, n
    real(wp), intent(in), dimension(0:nr+1,nsq,2)     :: psi
    real(wp), intent(out), dimension(nr,2)            :: rho                  ! density array
    real(wp), intent(out), dimension(nr)              :: rho_tot              ! total density array
!
! Initialize:
    rho     = 0.0_wp
    rho_tot = 0.0_wp
    A       = 0.0_wp
!
! Loop to calculate density, rho(i, iq):
    do iq = 1, 2
        do is = 1, tot_lvl(iq)
            do i = 1, nr
            A(i,iq) = (occ(is,iq)/(4.0_wp * pi * r(i)**2) ) * abs(psi(i,is,iq))**2
            rho(i,iq) = rho(i,iq) + A(i,iq)
          end do
        end do
    end do
!
! Total density rho_tot:
    rho_tot    = rho(1:nr,1) + rho(1:nr,2)
!
! Integrate rho(i,1), rho(i,2) to find neutron, proton number, and calculate rho total, rho_tot:
!
    sum_n      = 4.0_wp * pi * sum(rho(:,1) * r**2) * dr
!
    sum_p      = 4.0_wp * pi * sum(rho(:,2) * r**2) * dr
!
    sum_t      = 4.0_wp * pi * sum(rho_tot * r**2) * dr

!
  end subroutine density
!===============================================================================================================================================
! ===============================================================================================================================================
  subroutine Coulomb(r, nr, dr, phi_c, rho, sum_p)
  implicit none

! Local variables:
    integer, parameter      :: wp = kind(1.0d0)
    real(wp), parameter     :: pi = acos(-1.0d0)
    real(wp), parameter     :: e_2  = 1.4399643929         ! charge of a proton**2 w/units MeV*fm
    integer                 :: i, is, iq

! Global variables:
    integer, intent(in)     :: nr
    real(wp), intent(in)    :: dr
    real(wp), intent(in)    :: sum_p                       ! number of protons Z

! Local arrays
    real(wp)                :: F(0:nr+1), E(0:nr+1), Wc(0:nr+1)
    real(wp)                :: di(nr)

! Global arrays:
    real(wp), intent(in)    :: r(nr)
    real(wp), intent(in)    :: rho(nr,2)
    real(wp), intent(out)   :: phi_c(nr)

!
!
! Initialize and establish boundary conditions:
    di   = 0.0_wp
    E    = 0.0_wp
    F    = 0.0_wp
    Wc   = 0.0_wp
    Wc(nr+1) = sum_p * e_2

! Loop for arrays, E, F, and di:
      do is = 1, nr
        do i = 1, nr
            di(i) = -4.0_wp * pi * e_2 * rho(i,2) * r(i) * dr**2
            E(i) = -1.0_wp / (E(i-1) - 2.0_wp)
            F(i) = (di(i) - F(i-1)) / (E(i-1) - 2.0_wp)
        end do
      end do

! Backwards loop to fill array Wc using arrays E and F:
      do i = nr, 1, -1
        Wc(i) = E(i) * Wc(i+1) + F(i)
      end do

! Loop to calculate phi_c (Coulomb potential) using Wc and r:

      do i = 1, nr
        phi_c(i) = Wc(i) / r(i) - ( ((3.0_wp/pi)**(1.0_wp/3.0_wp)) * e_2 * rho(i,2)**(1.0_wp/3.0_wp))
      end do

  end subroutine Coulomb
!==========================================================================================================================================================
!==========================================================================================================================================================
subroutine yukawa(r, nr, pi, phi_y, rho, dr)
implicit none
!---------------------------------------------------------------------------------------------------------
! Calculation of the Yukawa potential, in a spherically symmetrical system,
! accomplished by a modified solution to the Helmholtz equation:
! (Laplacian - mu**2)*Phi(r) = -4*pi*Vy*rho(r). In this case we'll set Vy = 1.
!-----------------------------------------------------------------------------------------------------------

! Local variables:
    integer , parameter     :: wp = kind(1.0d0)
    integer                 :: i, is, iq
    real(wp)                :: bi, mu, Vy

! Global variables:
    integer, intent(in)     :: nr
    real(wp), intent(in)    :: pi
    real(wp), intent(in)    :: dr


! Internal arrays:
    real(wp)                :: Ui(0:nr,2), di(nr,2)
    real(wp)                :: Ei(0:nr+1,2), Fi(0:nr+1,2)

! Global arrays:
    real(wp), intent(in)    :: r(nr)
    real(wp), intent(out)   :: phi_y(nr,2)
    real(wp), intent(in)    :: rho(nr,2)

! Initialize and set boundary conditions:
    mu  =  2.174906_wp
    bi  = -mu**2 * dr**2 - 2.0_wp
    Vy  = -166.9239_wp
!
    Ei    = 0.0_wp
    Fi    = 0.0_wp
    Ui    = 0.0_wp
    phi_y = 0.0_wp

! Loop for array di:
      do iq = 1, 2
        do i = 1, nr
          di(i, iq) = -4.0_wp * pi * Vy * rho(i, iq) * r(i) * (dr**2)
        end do
      end do

! Backwards loop for arrays Ei , Fi:
      do iq = 1, 2
        do i = nr, 1, -1
          Ei(i, iq) = -1.0_wp/(bi + Ei(i+1, iq))
          Fi(i, iq) = (di(i, iq) - Fi(i+1, iq))/(bi + Ei(i+1, iq))
        end do
      end do

! Loop to calculate Ui:
      do iq = 1, 2
        do i = 1, nr
          Ui(i, iq) = (Ei(i, iq) * Ui(i-1, iq)) + Fi(i, iq)
        end do
      end do

! Calculate contribution from Yukawa potential, phi_y, using Ui and r:
      do iq = 1, 2
        do i = 1, nr
          phi_y(i, iq) = Ui(i, iq)/r(i)
        end do
      end do

end subroutine yukawa
!=========================================================================================================================================================
!=========================================================================================================================================================
   subroutine BKN(r, phi_c, phi_y, nr, dr, rho, rho_tot, UBKN)
    implicit NONE

!  Local variables:
    integer, parameter       :: wp = kind(1.0d0)
    real(wp), parameter      :: pi = acos(1.d0)
    real(wp), parameter      :: t0 = -497.726_wp, t3 = 17270.0_wp
    integer                  :: i, is, iq

!  Global variables:
    integer,  intent(in)     :: nr                    ! number of points
    real(wp), intent(in)     :: dr                    ! rstep

!  Global arrays:
    real(wp), dimension(nr), intent(in)          :: r              ! radial array
    real(wp), dimension(nr), intent(in)          :: rho_tot        ! proton + neutron density
    real(wp), dimension(nr), intent(in)          :: phi_c          ! coulomb potential
    real(wp), dimension(nr,2), intent(in)        :: phi_y          ! yukawa potential
    real(wp), dimension(nr,2), intent(in)        :: rho            ! neutron or proton density
    real(wp), dimension(nr,2), intent(out)       :: UBKN           ! U(nr,nsq,1)= neutron, U(nr,nsq,2) = proton
!
! Initialize:
      UBKN = 0.0_wp
!
! Loop for UBKN(nr,2):
      do i = 1, nr
!     for neutrons:
        UBKN(i, 1) = (3.0_wp/2.0_wp) * t0 * rho(i, 2) + (1.0_wp/4.0_wp) * t3 * (rho_tot(i)**2 - rho(i, 1)**2) &
          + phi_y(i, 1) + phi_y(i, 2)
!     for protons:
        UBKN(i, 2) = (3.0_wp/2.0_wp) * t0 * rho(i, 1) + (1.0_wp/4.0_wp) * t3 * (rho_tot(i)**2 - rho(i, 2)**2) &
          + phi_y(i, 2) + phi_y(i, 1) + phi_c(i)
      end do
!
  end subroutine BKN
!======================================================================================================================================================
!======================================================================================================================================================
  subroutine spet(r, nr, nsq, l, dr, tot_lvl, psi, UBKN, spe, tpe)
    implicit none
!
! Local variables:
    integer, parameter    :: wp   = kind(1.0d0)
    real(wp), parameter   :: beta = 20.726_wp                          ! represents hbar*2/(2 * mass_iq)--beta(1) for neutron, beta(2) for proton
    integer               :: i, is, iq
    real(wp)              :: sum1, sum2
!
! Global variables:
    integer, intent(in)   :: nr, nsq
    real(wp), intent(in)  :: dr
!
! Global arrays:
    real(wp), intent(in)                            :: r(nr)
    integer, intent(in)                             :: l(nsq,2)
    integer , intent(in)                            :: tot_lvl(2)
    real(wp), intent(in)                            :: psi(0:nr+1,nsq,2)
    real(wp), intent(in)                            :: UBKN(nr,2)
    real(wp), allocatable, intent(out)              :: spe(:,:)          ! spe(is,1) = neutron energy & spe(is,2) = proton energy
    real(wp), allocatable, intent(out)              :: tpe(:,:)          ! tpe(is,1) = neutron kinetic energy & tepe(is,2) = proton kinetic energy
!
! Local array:
    real(wp), dimension(nr, nsq, 2)                 :: D, F        ! Stand in arrays for when I need to sum for spe(ns,2) or tpe(ns,2)
!
    allocate( spe(nsq,2) )
    allocate( tpe(nsq,2) )
!
    spe  = 0.0_wp
    tpe  = 0.0_wp
    sum1 = 0.0_wp
    sum2 = 0.0_wp
!
!                                C A L C U L A T E  S I N G L E  P A R T I C L E  E N E R G I E S :
!
! Total energy spe(nsq,2):
    do iq = 1, 2
          do is = 1, tot_lvl(iq)
            do i = 1, nr

            D(i, is, iq) =  ( (-1.0_wp) * beta/dr**2 * psi(i,is,iq) * psi(i - 1, is, iq) + &
                        ((beta/(r(i)**2) * (l(is,iq)**2 + l(is,iq))) + beta * (2.0_wp/dr**2) + UBKN(i,iq)) * psi(i,is,iq)**2  &
                        -  (beta/dr**2 * psi(i,is,iq) * psi(i + 1, is, iq) ))

                                  spe(is, iq) = spe(is, iq) + dr * D(i, is, iq)    ! single particle energy array

            end do
          end do
      end do

       do iq = 1, 2
         do is = 1, tot_lvl(iq)
           do i = 1, nr

             !  Kinetic energy tpe(nsq,2):
              F(i, is, iq) =  ( (-1.0_wp) * beta/dr**2 * psi(i,is,iq) * psi(i - 1, is, iq) + &
                                   ((beta/(r(i)**2) * (l(is,iq)**2 + l(is,iq))) + beta * (2.0_wp/dr**2)) * psi(i,is,iq)**2 &
                                    -  (beta/dr**2 * psi(i,is,iq) * psi(i + 1, is, iq) ))

                                   tpe(is, iq) = tpe(is, iq) + dr * F(i, is, iq)     ! single particle kinetic energy array
            end do
         end do
      end do

end subroutine spet
!=========================================================================================================================================================
!=========================================================================================================================================================
    subroutine grstep(r, l, tot_lvl, T, psi, psitmp, spe, dr, nr, nsq, UBKN, lst, nsl, il, m, blank)
    implicit none

! Local variables:
    integer, parameter    :: wp = kind(1.0d0)
    integer               :: i, j, is
    integer, intent(out)  :: il, m
    real(wp), parameter   :: X0 = 0.02_wp
    real(wp), parameter   :: beta = 20.726                          ! represents hbar*2/(2 * mass_iq)--beta(1) for neutron, beta(2) for proton

! Global variables:
    integer, intent(in)   :: nr, nsq, nsl
    real(wp), intent(in)  :: dr
    integer , intent(out) :: blank

! Local arrays:
    real(wp), dimension(nr,nsq,2) :: C
    real(wp), dimension(nr,nsq,2) :: P
    real(wp), allocatable :: summ(:,:)

! Global arrays:
    real(wp), intent(in)    :: r(nr)
    integer, intent(inout)  :: l(nsq,2)
    integer, intent(in)     :: tot_lvl(2)
    integer, intent(in)     :: lst(0:nsl,2)
    real(wp), intent(inout) :: psi(0:nr+1,nsq,2)
    real(wp), intent(in)    :: spe(nsq,2)
    real(wp), intent(in)    :: T(nr,nr,nsq,2)
    real(wp), intent(in)    :: UBKN(nr,2)
    real(wp), allocatable, intent(out) :: psitmp(:,:,:)

!
! Initialize and set bundary conditions psitmp(0,:,:) = psitmp(nr+1,:,:) = 0.0 :
    allocate ( psitmp(0:nr+1, nsq, 2) )
    allocate ( summ(nsq, 2) )

    C      = 0.0_wp
    M      = 0.0_wp
    P      = 0.0_wp
    psitmp = 0.0_wp
    summ   = 0.0_wp

    do iq = 1, 2
      do is = 1, tot_lvl(iq)
        do i = 1, nr
          C(i,is,iq) = ( (-1.0_wp) * beta/dr**2 * psi(i - 1, is, iq) + &
                       ( (beta/(r(i)**2) * (l(is,iq)**2 + l(is,iq))) + beta * (2.0_wp/dr**2) + UBKN(i,iq)) * psi(i,is,iq)  &
                      -( beta/dr**2 * psi(i + 1, is, iq) )) &
                      - spe(is,iq)*psi(i,is,iq)
        end do
     end do
    end do
!
    do iq = 1, 2
      do is = 1, tot_lvl(iq)
          psitmp(1:nr, is, iq) =  psi(1:nr, is, iq) - X0 * matmul(T(:,:,is,iq),C(:,is,iq))
      end do
    end do

        do iq = 1, 2                                        ! Normalize the wave function
      do is = 1, tot_lvl(iq)                            ! for different states/fermions
        summ(is,iq) = sum(psitmp(1:nr,is,iq)**2) * dr
        psitmp(1:nr,is,iq) = psitmp(1:nr,is,iq) / sqrt(summ(is,iq))
      end do
    end do


!
!
    psi = psitmp

    do il = 1, 2
    blank = 0
      do i = 0, maxval(l(:,il))
              m = lst(i,il)
              if (m .gt. 1) then
              call gramschmid(m, psi, psitmp, nr, nsq, il, lgroups, blank)
              end if
              if (m .gt. 1) blank = blank + m
          end do
      end do

      do iq = 1, 2                                        ! Normalize the wave function
      do is = 1, tot_lvl(iq)                            ! for different states/fermions
        summ(is,iq) = sum(psi(1:nr,is,iq)**2) * dr
        psi(1:nr,is,iq) = psi(1:nr,is,iq) / sqrt(summ(is,iq))
      end do
    end do

    end subroutine grstep
!=======================================================================================================================================================================
!=======================================================================================================================================================================
    subroutine gramschmid(m, psi, psitmp, nr, nsq, il, lgroups, blank)
    implicit none

    integer, parameter :: wp = kind(1.0d0)

! Local variables::
    integer :: i, j, k
    integer :: nsi

! Local arrays:
    integer, allocatable  :: is(:,:)
    real(wp), allocatable :: lower_product(:,:)
    real(wp), allocatable :: upper_product(:,:,:)

! Global variables:
    integer, intent(in) :: nr, nsq
    integer, intent(in) :: m
    integer, intent(in) :: il
    character (len=16), dimension(2), intent(in)     :: lgroups

! Global arrays:
    integer, intent(in)     :: blank
    real(wp), intent(inout) :: psi(0:nr+1,nsq,2), psitmp(0:nr+1,nsq,2)
    real(wp), allocatable   :: summ(:,:)

    allocate ( is(m, 2) )
    allocate ( lower_product(m,2) )
    allocate ( upper_product(m,m,2) )
    allocate ( summ(m,2))

    lower_product = 0.0_wp
    upper_product = 0.0_wp
    is = 0
    summ = 0

!       do j = 1, 2
      open (unit = 12, file = lgroups(il))
      do k = 1, blank
      read(12, fmt = '(a)')
      end do
      do i = 1, m
      read(12,fmt = 130)is(i,il)
      130 format (10x, i2)
      end do
      close(12)
!
    nsi = m
        do j = m, 1, -1
          do i = 1, m
            if ((j - i) .eq. 0 ) exit
            psi(1:nr,is(j,il), il) = psitmp(1:nr,is(j,il),il) - (dot_product(psi(1:nr,is(j,il),il), psi(1:nr,is(j-i,il),il)*dr ) &
                                   / (norm2(psi(1:nr,is(j-i,il),il))**2 *dr ) ) * psi(1:nr,is(j-i,il),il)
            psitmp(1:nr,is(j,il),il) = psi(1:nr,is(j,il),il)
            end do
        end do
!
        do i = 1, m
          summ(i,il) = sqrt(sum(psi(1:nr, is(i,il), il)**2) * dr)
          psi(1:nr,is(i,il),il) = psi(1:nr,is(i,il),il)/summ(i,il)

        end do
!
    deallocate (is)
    deallocate (lower_product)
    deallocate (upper_product)
    deallocate (summ)

    end subroutine gramschmid
!============================================================================================================================================================================
!============================================================================================================================================================================

  subroutine SHF(r, l, tot_lvl,nr, nsq, dr, rho, rho_tot, spe, tpe, EHF)
  implicit none

! Local variables:
  integer, parameter    :: wp = kind(1.0d0)
  real(wp), parameter   :: pi = acos(-1.0d0)
  real(wp), parameter   :: t3 = 17270.0_wp
  real(wp)              :: summ
  integer               :: i, is, iq
  integer               :: k
  real(wp)              :: C3


! Global variables:
  integer, intent(in)   :: nr, nsq
  real(wp), intent(in)  :: dr

! Global arrays:
  real(wp), intent(in)  :: r(nr)
  integer, intent(in)   :: tot_lvl(2)
  integer, intent(in)   :: l(nsq,2)
  real(wp), intent(in)  :: rho(nr,2), rho_tot(nr)
  real(wp), intent(in)  :: spe(nsq,2), tpe(nsq,2)
  real(wp), intent(out) :: EHF

! Local arrays:
  real(wp)              :: D(nr,2), P(nsq,2)

  summ = 0.0_wp
  C3   = 0.0_wp
  EHF  = 0.0_wp

  C3 = sum(r**2 * rho(1:nr,1) * rho(1:nr,2) * rho_tot * dr * 4.0_wp * pi * (1.0_wp/8.0_wp) * t3)

  do iq = 1, 2
    do is = 1, tot_lvl(iq)
      summ = summ + (occ(is,iq)/2.0_wp) * (spe(is,iq) + tpe(is,iq))
    end do
  end do


  EHF = summ - C3

end subroutine SHF
!==========================================================================================================================================================
!==========================================================================================================================================================
    end program

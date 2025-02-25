!--------------------------------------------------------------
module constants
!--------------------------------------------------------------
    ! Universal constants (in SI units)
    double precision, parameter :: G = 6.674d-11
    double precision, parameter :: pi = 4.d0 * atan(1.d0)
    double precision, parameter :: c = 2.998d8
    double precision, parameter :: h = 6.626d-34
    double precision, parameter :: sigma = 5.670d-8
    double precision, parameter :: kb = 1.381d-23
    
    ! Conversion factors (SI to given units)
    double precision, parameter :: day = 3600*24
    double precision, parameter :: Msun = 1.988d30
    double precision, parameter :: erg = 1d7
    double precision, parameter :: lday = 2.998d8 * 3600*24
    double precision, parameter :: Lsun = 3.828d26
    double precision, parameter :: Ang = 1d-10
    double precision, parameter :: year = 3600 * 24 * 365.25
    double precision, parameter :: parsec = 3.256 * c * year
    double precision, parameter :: erg_mJy = 1d26
end module constants


!--------------------------------------------------------------
function Planck(freq, temp)
! Computes black body intensity (J/m^2/s/Hz/ster) at given
! temperature and wavelength
! Parameters - 
! freq (double precision) : frequency (Hz)
! temp (double precision) : temperature (K)
! Returns - 
! Plank (double precision) : black body intensity
!--------------------------------------------------------------
    use constants
    implicit none
    save
    double precision, intent(in) :: freq
    double precision, intent(in) :: temp
    double precision :: x, c1, c2
    logical :: firstcall = .true.
    double precision :: Planck

    if (firstcall) then
        c1 = (2 * kb**3)/(c**2 * h**2) ! (2kb^3)/(c^2h^2)
        c2 = h / kb ! h/kb
        firstcall = .false.
    end if

    x = c2 * freq / temp

    if (x .lt. 1.d-5) then ! Rayleigh-Jeans tail
        Planck = c1*(x**2)*(temp**3)
    else
        Planck = c1*(x**3)*(temp**3) / (exp(x)-1)
    end if
end function


!--------------------------------------------------------------
subroutine SED(nr, r_r, dr_r, T_r, nw, wl, incl, Lnu)
! Computes the Disc SED (with self-occultation)
! Parameters - 
! nr (integer) : number of annular elements
! r_r(nr) (double precision) : radial elements (in m)
! h_r(nr) (double precision) : height of the disc (not used here)
! dr_r(nr) (double precision) : radial thickness of the annulus (in m)
! T_r(nr) (double precision) : Temperature of annulus (in K)
! nw (integer) : number of wavelenghts
! wl(nw) (double precision) : wavelength range
! incl (double precision) : inclination of accretion disc
! Returns -
! Lnu(nw) (double precision) : Lnu (erg/s/Hz)
!--------------------------------------------------------------
    use constants
    implicit none
    integer, intent(in) :: nr, nw
    double precision, intent(in) :: r_r(nr), dr_r(nr), T_r(nr)
    double precision, intent(in) :: wl(nw)
    double precision, intent(in) :: incl
    double precision, intent(out) :: Lnu(nw)

    double precision, external :: Planck

    integer :: ir, iaz, iw
    integer :: n_az = 360
    double precision :: dr, dh, dx, dA, daz
    double precision :: rmid, sum
    double precision :: ex, ez
    double precision :: Ax, Az!, Ay
    double precision :: dot, sintilt, costilt
    double precision :: freq
    double precision :: Bnu

    double precision :: h_mid, tan_incl
    double precision :: phi, delta_H, delta_R
    logical :: vis


    if (incl .eq. 0d0) n_az = 1 
    daz = 2*pi/n_az

    ! Unit vector towards Earth
    ex = sin(incl)
    ez = sqrt(1 - ex*ex)
    
    tan_incl = ex/ez

    ! Iterate through wavelengths
    do iw=1, nw
        sum = 0.d0
        ! Iterate through annuli
        do ir=1, nr
            ! dr = r_r(ir+1) - r_r(ir)
            ! rmid = (r_r(ir+1) + r_r(ir))/2
            rmid = r_r(ir)
            ! dh = h_r(ir+1) - h_r(ir)

            ! dx = sqrt(dr**2 + dh**2)
            dx = dr_r(ir)
            dA = dx * rmid * daz
            
            ! sintilt = dh/dx
            ! costilt = dr/dx
            sintilt = 0d0
            costilt = 1d0

            freq = c/wl(iw)
            Bnu = Planck(freq, T_r(ir))

            ! h_mid = (h_r(ir+1) + h_r(ir)) / 2
            ! delta_H = h_r(nr) - h_mid
            ! delta_R = delta_H * tan_incl

            ! Iterate through azimuthal angles
            do iaz=1, n_az

                phi = iaz*daz
                ! call self_occulation_check(rmid, phi, delta_R, r_r(nr), vis)

                ! if (.not. vis) cycle

                ! Normal to surface element
                Ax = -cos(iaz*daz) * sintilt
                !Ay = sin(iaz*daz) * sintilt
                Az = costilt

                dot = Ax*ex + Az*ez ! + Ay*ey (=0)

                if (dot .gt. 0d0) then
                    sum = sum + Bnu * dA * dot * erg
                end if
            end do
        end do
        Lnu(iw) = sum
        ! Note - output is in cgs, so to convert to flux at earth, do
        !       Lnu = Lnu / distance / distance

        ! Also, distance must be in cm
    end do
end subroutine SED
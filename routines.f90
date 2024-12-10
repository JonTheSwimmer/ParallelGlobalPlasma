!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    module routines
    use parameters_base
    implicit none

    contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine integrate_dlsode(Neq,X,tstart,tstop,derivs)
    implicit none
    integer, intent(in) :: neq
    real*8, intent(inout) :: X(neq)
    real*8, intent(in) :: tstart, tstop
    integer, parameter :: nmax   = 2000
    integer, parameter :: neqmax =  30
    real*8, parameter :: prec = tolerance
    integer ::  itol, itask, istate, iopt, lrw, iwork(nmax), liw, mf
    real*8 ::  t1, t2, rtol, atol, rwork(nmax)!, y(neqmax)
    external :: derivs

    !solver parameters
    !iopt = 1 means were are inputting expanded maximum stepsize
    itol=1; rtol=prec; atol=rtol; itask=1; istate=1; iopt=1
    lrw=nmax; liw=nmax; mf=10

    rwork(6) = 0d0; rwork(5) = 0d0; rwork(7) = 0d0; rwork(8) = 0d0; &
    rwork(9) = 0d0; rwork(10) = 0d0
    iwork(6) = 30000000; iwork(5) = 0; iwork(7) = 1; iwork(8) = 0; &
    iwork(9) = 0; iwork(10) = 0

    t1 = tstart;  t2 = tstop
    call dlsode(derivs,neq,X,t1,t2,itol,rtol,atol,itask, &
                  istate,iopt,rwork,lrw,iwork,liw,jacobian,mf)

    return
    end subroutine integrate_dlsode

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Attempt to do an RK4 integrator so I can parallelize things
!   Uses the RKF4(5) method to determine the adaptive step size
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine integrate_ode(Neq,X,tstart,tstop,derivs)
    implicit none
    integer, intent(in) :: neq
    real*8, intent(inout) :: X(neq)
    real*8, intent(in) :: tstart, tstop
    external :: derivs
    !local variables
    real*8, dimension(neq) :: yn, yn2, yerr, k1, k2, k3, k4, k5, k6
    real*8 :: h, h_init, err, tn
    real*8 :: a3, a4
    real*8 :: b31, b32, b41, b42, b43, b51, b52, b53, b54, b61, b62, b63, b64, b65
    real*8 :: c11, c13, c14, c15, c16, c21, c23, c24, c25, c26
    logical :: cont
    real*8 :: prec_arr(neq), h_scale
    tn = tstart; yn = X
    h_init = (tstop - tstart)/8
    a3 = 3.d0/8.d0; a4 = 1.2d1/1.3d1
    b31 = 3.d0/3.2d1; b32 = 9.d0/3.2d1
    b41 = 1.932d3/2.197d3; b42 = -7.2d3/2.197d3; b43 = 7.296d3/2.197d3
    b51 = 4.39d2/2.16d2; b52 = -8.d0; b53 = 3.68d3/5.13d2; b54 = -8.45d2/4.104d3
    b61 = -8.d0/2.7d1; b62 = 2.d0; b63 = -3.544d3/2.565d3; b64 = 1.859d3/4.104d3; b65 = -1.1d1/4.d1
    c11 = 1.6d1/1.35d2; c13 = 6.656d3/1.2825d4
    c14 = 2.8561d4/5.643d4; c15 = -9.d0/5.d1; c16 = 2.d0/5.5d1
    c21 = -1.d0/3.6d2; c23 = 1.28d2/4.275d3; c24 = 2.197d3/7.524d4
    c25 = -1.d0/5.d1; c26 = -2.d0/5.5d1
    do while (tn .lt. tstop)
        h = min(h_init, tstop - tn)
        cont = .true.
        do while (cont)
            call derivs(neq, tn, yn, k1)
            call derivs(neq, tn+ (h/4.d0), yn + (h/4.d0) * k1, k2)
            call derivs(neq, tn+ h*a3, yn + h*(b31*k1 + b32*k2), k3)
            call derivs(neq, tn + h*a4, yn + h*(b41*k1 + b42*k2 + b43*k3), k4)
            call derivs(neq, tn + h, yn + &
                h*(b51*k1 + b52*k2 + b53*k3 + b54*k4), k5)
            call derivs(neq, tn + 0.5d0*h, yn + &
                h*(b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5), k6)
            yn2 = yn + h*(c11*k1 + c13*k3 + c14*k4 + c15*k5 + c16*k6)
            yerr = h*(c21*k1 + c23*k3 + c24*k4 + c25*k5 + c26*k6)
            !err = norm2(yerr)
            prec_arr = tolerance * (1.d0 + abs(yn2))
            h_scale = maxval(abs(yerr/prec_arr))
            if (h_scale .le. 1.d0) then
                cont = .false.
            else
                h = 0.9d0 * h * (h_scale)**(-0.2d0)
            end if
        end do
        !h_init = min((tstop - tstart)/8, h * (h_scale)**(0.2d0))
        tn = tn + h; yn = yn2
        !print *, tstop - tn, maxval(abs(yerr)), yerr(9), h_scale, h
    X = yn
    end do
    end subroutine integrate_ode
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine jacobian
    implicit none

    return
    end subroutine jacobian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine spline(x,y,n,yp1,ypn,y2)
    implicit none
    integer :: n
    real*8 :: yp1,ypn,x(n),y(n),y2(n)
    integer, parameter :: NMAX = 100000
    integer :: i,k
    real*8 :: p,qn,sig,un,u(NMAX)

    if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
    else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    endif
    do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    end do
    if (ypn.gt..99e30) then
        qn=0.
        un=0.
    else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    endif
    y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
    do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
    end do

    return
    end subroutine spline

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine splint(xa,ya,y2a,n,x,y,abscissa_increasing)
    implicit none
    integer :: n
    real*8 :: x,y,xa(n),y2a(n),ya(n)
    logical :: abscissa_increasing
    integer :: k,khi,klo
    real*8 :: a,b,h

    klo=1
    khi=n
1   if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
            if (abscissa_increasing) then
                khi=k
            else
                klo=k
            endif
        else
            if (abscissa_increasing) then
                klo=k
            else
                khi=k
            endif
        endif
    goto 1
    endif
    h=xa(khi)-xa(klo)
    if (h.eq.0.) print *, xa(khi), khi, klo
    if (h.eq.0.) pause 'bad xa input in splint'
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
!      print *, a, b, xa(khi), xa(klo), x, y

    return
    end subroutine splint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine splint_arr(xa,ya,y2a,n,x,y,nx)
    !spline procedure for array inputs
    !xa must be a linearly uniform array and monotonically increasing!
    !otherwise use splint
    implicit none
    integer :: n, nx
    real*8 :: xa(n), ya(n), y2a(n), x(nx), y(nx)
    integer, dimension(nx) :: khi, klo
    real*8, dimension(nx) :: a, b
    real*8 :: step, tol
    integer :: i
    tol = 1.d-8
    step = xa(2) - xa(1)
    khi = max(2, ceiling((x - xa(1)) / step))
    klo = khi - 1
    a = (xa(khi) - x) / step
    b = 1 - a
    y =a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(step**2)/6.
    return
    end subroutine splint_arr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function ran2(idum_in)
    integer, intent(inout) :: idum_in
    INTEGER :: IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
    REAL :: ran2,AM,EPS,RNMX
    PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
    IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
    NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
    INTEGER :: idum2,j,k,iv(NTAB),iy
    SAVE :: iv,iy,idum2
    DATA idum2/123456789/, iv/NTAB*0/, iy/0/

    if (idum_in <= 0) then
        idum_in=max(-idum_in,1)
        idum2=idum_in
        do j=NTAB+8,1,-1
            k=idum_in/IQ1
            idum_in=IA1*(idum_in-k*IQ1)-k*IR1
            if (idum_in < 0) idum_in=idum_in+IM1
            if (j <= NTAB) iv(j)=idum_in
        end do
        iy=iv(1)
    endif
    k=idum_in/IQ1
    idum_in=IA1*(idum_in-k*IQ1)-k*IR1
    if (idum_in < 0) idum_in=idum_in+IM1
    k=idum2/IQ2
    idum2=IA2*(idum2-k*IQ2)-k*IR2
    if (idum2 < 0) idum2=idum2+IM2
    j=1+iy/NDIV
    iy=iv(j)-idum2
    iv(j)=idum_in
    if(iy < 1)iy=iy+IMM1
    ran2=min(AM*iy,RNMX)

    return
    end function ran2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! for a monotonic array, return the index of the last element that is less than/greater than the point
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function bisect(arr, point, length, array_increasing)
    integer, intent(in) :: length
    real*8, dimension(length), intent(in) :: arr
    real*8, intent(in) :: point
    logical, intent(in) :: array_increasing
    integer :: ind, bisect
    integer :: upper, lower
    logical :: cont
    upper = length + 1
    lower = 1
    cont = .true.
    ind = floor(upper*1.d0/2)
    if (array_increasing) then
        do while (cont)
            if (arr(ind) < point) then
                lower = ind
            else if (arr(ind) > point) then
                upper = ind
            else
                bisect = ind
                return
            end if
            ind = floor((upper + lower)*1.d0 / 2)
            if ((upper - lower) == 1) then
                cont = .false.
            end if
        end do
    else
        do while (cont)
            if (arr(ind) > point) then
                lower = ind
            else if (arr(ind) < point) then
                upper = ind
            else
                bisect = ind
                return
            end if
            ind = floor((upper + lower)*1.d0 / 2)
            if ((upper - lower) == 1) then
                cont = .false.
            end if
        end do
    endif
    bisect = ind
    return
    end function bisect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   cross product
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine cross(a, b, c)
    real*8, dimension(3), intent(in) :: a, b
    real*8, dimension(3), intent(out) :: c

    c(1) = a(2)*b(3) - b(2)*a(3)
    c(2) = a(3)*b(1) - b(3)*a(1)
    c(3) = a(1)*b(2) - b(1)*a(2)

    return
    end subroutine cross

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   convert spherical to cartesian
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine cartesian(R, theta, phi, x, y, z)
    real*8, intent(in) :: R, theta, phi
    real*8, intent(out) :: x, y, z
    x = R * sin(theta) * cos(phi)
    y = R * sin(theta) * sin(phi)
    z = R * cos(theta)
    return
    end subroutine cartesian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   convert cartesian to spherical
    subroutine spherical(x, y, z, R, theta, phi)
    real*8, intent(in) :: x, y, z
    real*8 :: rho
    real*8, intent(out) :: R, theta, phi
    R = (x**2 + y**2 + z**2)**0.5
    rho = (x**2 + y**2)**0.5
    if (rho == 0) then
        if (z >= 0) then
            theta = 0
            phi = 0
            return
        else
            theta = pi
            phi = 0
            return
        end if
    else
        theta = acos(z/R)
        if (y == 0) then
            if (x >= 0) then
                phi = 0
            else
                phi = pi
            end if
        else
            phi = acos(x/rho)
            if (y < 0) then
               phi = twopi - phi
            end if
        end if
    end if

    return
    end subroutine spherical


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Shorthand for opening and checking if a file already exists
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine file_open(unitnum, filename)
    implicit none
    integer, intent(in) :: unitnum
    character(len=*), intent(in) :: filename
    logical :: exist
    inquire(file=filename, exist = exist)
    if (exist) then
        open(unit = unitnum, file = filename,Status='old', position="append", action="write")
    else
        open(unit = unitnum, file = filename, status = 'new')
    endif
    return
    end subroutine file_open
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end module routines

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

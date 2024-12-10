      program main
      implicit none

      integer, parameter :: nmax   = 200
      integer, parameter :: neqmax =  10 
      real(8), parameter :: prec   = 1d-12 

      integer  itol,itask,istate,iopt,lrw,iwork(nmax),liw,mf,neq
      real(8)  xstart,xstop,rtol,atol,rwork(nmax),y(neqmax)
      external derivs, jacobian

      !solver parameters
      itol=1; rtol=prec; atol=rtol; itask=1; istate=1; iopt=0
      lrw=nmax; liw=nmax; mf=10
   
      !integrate cos(x) from x = 0 to 
      xstart = 0. 
      xstop  = acos(-1d0)/2d0
      neq    = 1
      y(1)   = 0.

      call dlsode(derivs,neq,y,xstart,xstop,itol,rtol,atol,itask, 
     &              istate,iopt,rwork,lrw,iwork,liw,jacobian,mf)

      print *, 'y(pi/2) =', y(1)    !answer should be 1
 
      end program main

      !-------------------------
      subroutine derivs(nvar,x,y,dydx)
      implicit none
      
      integer nvar
      real(8) x, y(nvar), dydx(nvar)

      dydx(1) = cos(x)

      return
      end subroutine derivs

      !--------------------------
      subroutine jacobian
      implicit none

      return
      end subroutine jacobian

       module mod_functions
!===================================================================
!  This module contains all the functions of the code
!
!  Last revision: 9/April/2019
!===================================================================
!
!--Common variables definitions for the module
!
       use mod_parameters, only : pi
!
       contains
!
       real function wk(r)
!=============================================================
!      KERNEL function
!=============================================================
!
!--Force to declare EVERYTHinG
!
       implicit none
!
!--I/O variables
!
       real, intent(in) :: r
!
!--Local variables
!
       real :: u, sigma
!
       u=sqrt(r)
!
       if (u >= 0.0d0 .and. u < 1.0d0) then
          wk=1.0d0 - 1.5d0*u**2.0d0 + 0.75d0*u**3  
       elseif (u >= 1.0d0 .and. u < 2.0d0) then
          wk=0.25d0*(2.0d0-u)**3.0d0
       elseif (u >= 2.0d0) then
          wk=0.0d0
       endif
!
!--Include normalization factor
!
       wk = wk/pi
!
       end function wk
!
       real function dk(r)
!=======================================================================
!      KERNEL DERIVATIVE function
!
!      Remember that the kernel gradient can been expresed as:
!
!      grad W(r) = (dW/dr_x,dW/dr_y,dW/dr_z) = (r_x/r,r_y/r,r_z/r)*dW/dr
!                = (r_x,r_y,r_z)*[(1/r)*(dW/dr)]= (r_x,r_y,r_z)*F
!
!      Function dk(r) is nothing but function F !!!!
!=======================================================================
!
!--Force do declare EVERYTHinG
!
       implicit none
!
!--I/O variables
!
       real, intent(in) :: r    
!
!--Local variables
!
       real :: u, sigma
!
       u=sqrt(r)
!
       if (u >= 0.0d0.and. u < 1.0d0) then
          dk=-3.0d0 + 2.25d0*u
       elseif (u >= 1.0d0 .and. u < 2.0d0) then
          dk=-(0.750d0/u)*(2.0d0-u)**2
       elseif (u >= 2.0d0) then
          dk=0.0d0
       endif
!
!--Include normalization factor
!
       dk = dk/pi
!
       end function dk   
!
       real function wk_Q5(r)
!=============================================================
!      QUinTIC KERNEL function
!=============================================================
!
!--Force to declare EVERYTHinG
!
       implicit none
!
!--I/O variables
!
       real, intent(in) :: r
!
!--Local variables
!
       real :: u, sigma
!
       u=sqrt(r)
!
       if (u >= 0.0d0 .and. u < 1.0d0) then
          wk_Q5=(3.0-u)**5 - 6.0*(2.0-u)**5 + 15.0*(1.0-u)**5
       elseif (u >= 1.0d0 .and. u < 2.0d0) then
          wk_Q5=(3.0-u)**5 - 6.0*(2.0-u)**5
       elseif (u >= 2.0d0 .and. u < 3.0d0) then
          wk_Q5=(3.0-u)**5
       elseif (u >= 3.0d0) then
          wk_Q5=0.0d0
       endif
!
!--Include normalization factor
!
       wk_Q5 = (1./(120.*pi))*wk_Q5
!
       end function wk_Q5
!
       real function dk_Q5(r)
!=======================================================================
!      QUinTIC KERNEL DERIVATIVE function
!
!      Remember that the kernel gradient can been expresed as:
!
!      grad W(r) = (dW/dr_x,dW/dr_y,dW/dr_z) = (r_x/r,r_y/r,r_z/r)*dW/dr
!                = (r_x,r_y,r_z)*[(1/r)*(dW/dr)]= (r_x,r_y,r_z)*F
!
!      Function dk(r) is nothing but function F !!!!
!========================================================================
!
!--Force to declare EVERYTHinG
!
       implicit none
!
!--I/O variables
!
       real, intent(in) :: r
!
!--Local variables
!
       real :: u, sigma
!
       u=sqrt(r)
!
       if (u == 0.0) then
          dk_Q5=0.0
       elseif (u > 0.0d0 .and. u < 1.0d0) then
          dk_Q5=-(5./u)*(3.0-u)**4 + (30.0/u)*(2.0-u)**4 - (75.0/u)*(1.0-u)**4
       elseif (u >= 1.0d0 .and. u < 2.0d0) then
          dk_Q5=-(5./u)*(3.0-u)**4 + (30.0/u)*(2.0-u)**4
       elseif (u >= 2.0d0 .and. u < 3.0d0) then
          dk_Q5=-(5./u)*(3.0-u)**4
       elseif (u >= 3.0d0) then
          dk_Q5=0.0d0
       endif
!
!--Include normalization factor
!
       dk_Q5 = (1.0/(120.*pi))*dk_Q5
!
       end function dk_Q5
!
       end module mod_functions

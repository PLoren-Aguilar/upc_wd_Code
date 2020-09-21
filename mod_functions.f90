       MODULE mod_functions
!===================================================================
!
!  This module contains all the functions of the code
!
!  Last revision: 15/March/2015
!
!===================================================================

!
!--Common variables definitions for the module
!
       USE mod_parameters, ONLY : pi
!
       CONTAINS
!
       REAL FUNCTION wk(r)
!=============================================================
!      KERNEL FUNCTION
!=============================================================
!
!--Force to declare EVERYTHING
!
       IMPLICIT NONE
!
!--I/O variables
!
       REAL, INTENT(IN) :: r
!
!--Local variables
!
       REAL :: u, sigma
!
       u=SQRT(r)
!
       IF (u >= 0.0d0 .AND. u < 1.0d0) THEN
          wk=1.0d0 - 1.5d0*u**2.0d0 + 0.75d0*u**3  
       ELSEIF (u >= 1.0d0 .AND. u < 2.0d0) THEN
          wk=0.25d0*(2.0d0-u)**3.0d0
       ELSEIF (u >= 2.0d0) THEN
          wk=0.0d0
       ENDIF
!
!--Include normalization factor
!
       wk = wk/pi
!
       END FUNCTION wk
!
       REAL FUNCTION dk(r)
!=======================================================================
!      KERNEL DERIVATIVE FUNCTION
!
!      Remember that the kernel gradient can been expresed as:
!
!      grad W(r) = (dW/dr_x,dW/dr_y,dW/dr_z) = (r_x/r,r_y/r,r_z/r)*dW/dr
!                = (r_x,r_y,r_z)*[(1/r)*(dW/dr)]= (r_x,r_y,r_z)*F
!
!      Function dk(r) is nothing but function F !!!!
!=======================================================================
!
!--Force do declare EVERYTHING
!
       IMPLICIT NONE
!
!--I/O variables
!
       REAL, INTENT(IN) :: r    
!
!--Local variables
!
       REAL :: u, sigma
!
       u=SQRT(r)
!
       IF (u >= 0.0d0.AND. u < 1.0d0) THEN
          dk=-3.0d0 + 2.25d0*u
       ELSEIF (u >= 1.0d0 .AND. u < 2.0d0) THEN
          dk=-(0.750d0/u)*(2.0d0-u)**2
       ELSEIF (u >= 2.0d0) THEN
          dk=0.0d0
       ENDIF
!
!--Include normalization factor
!
       dk = dk/pi
!
       END FUNCTION dk   
!
       REAL FUNCTION wk_Q5(r)
!=============================================================
!      KERNEL FUNCTION
!=============================================================
!
!--Force to declare EVERYTHING
!
       IMPLICIT NONE
!
!--I/O variables
!
       REAL, INTENT(IN) :: r
!
!--Local variables
!
       REAL :: u, sigma
!
       u=SQRT(r)
!
       IF (u >= 0.0d0 .AND. u < 1.0d0) THEN
          wk_Q5=(3.0-u)**5 - 6.0*(2.0-u)**5 + 15.0*(1.0-u)**5
       ELSEIF (u >= 1.0d0 .AND. u < 2.0d0) THEN
          wk_Q5=(3.0-u)**5 - 6.0*(2.0-u)**5
       ELSEIF (u >= 2.0d0 .AND. u < 3.0d0) THEN
          wk_Q5=(3.0-u)**5
       ELSEIF (u >= 3.0d0) THEN
          wk_Q5=0.0d0
       ENDIF
!
!--Include normalization factor
!
       wk_Q5 = 3.0*wk_Q5/(359.0*pi)
!
       END FUNCTION wk_Q5
!
       REAL FUNCTION dk_h(r)
!==========================================================================
!      KERNEL DERIVATIVE FUNCTION WITH RESPECT TO h
!
!      This function is going to be used to include the grad-h terms
!      so in this case it is really dW/dh
!===========================================================================
!
!--Force to declare EVERYTHING
!
       IMPLICIT NONE
!
!--I/O variables
!
       REAL, INTENT(IN) :: r
!
!--Local variables
!
       REAL :: u, sigma
!
       u=SQRT(r)
!
       IF (u >= 0.0d0.AND. u < 1.0d0) THEN
          dk_h=-3.0d0+7.5d0*(u**2)-4.5d0*(u**3)
       ELSEIF (u >= 1.0d0 .AND. u < 2.0d0) THEN
          dk_h=0.75d0*(-((2.0d0-u)**3)+u*((2.0d0-u)**2))
       ELSEIF (u >= 2.0d0) THEN
          dk_h=0.0d0
       ENDIF
!
!--Include normalization factor
!
       dk_h = dk_h/pi
!
       END FUNCTION dk_h
!
       END MODULE mod_functions

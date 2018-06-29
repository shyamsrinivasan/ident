      subroutine all_flux(y, p, flux)

       implicit none
       real, dimension(0:1), intent(in) :: y
       real, dimension(0:3), intent(in) :: p
       real :: alpha, beta, delta, gamma
       real, dimension(0:3), intent(out) :: flux

       alpha = p(0)
       beta = p(1)
       delta = p(2)
       gamma = p(3)

       flux(0) = alpha*y(0)
       flux(1) = beta*y(0)*y(1)
       flux(2) = delta*y(0)*y(1)
       flux(3) = gamma*y(1)

      end subroutine all_flux

      subroutine ode(flux, rhs)
       
       implicit none
       real, dimension(0:3), intent(in) :: flux
       real, dimension(0:1), intent(out) :: rhs       

       rhs(0) = flux(0)-flux(1)
       rhs(1) = flux(2)-flux(3)

      end subroutine ode
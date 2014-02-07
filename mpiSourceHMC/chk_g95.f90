module chk_g95_m
   implicit none
   private
   public is_g95
   type t
      integer x
   end type t
   contains
      function is_g95()
         logical is_g95
         interface assignment(=)
            module procedure assign_t
         end interface assignment(=)
         type(t) y(3)

         y%x = (/1,2,3/)
         y = y((/2,3,1/))
         is_g95 = y(3)%x == 1
      end function is_g95

      elemental subroutine assign_t(lhs,rhs)
         type(t), intent(in) :: rhs
         type(t), intent(out) :: lhs

         lhs%x = rhs%x
      end subroutine assign_t
end module chk_g95_m
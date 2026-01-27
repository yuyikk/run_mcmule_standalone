function fortran_read_int() bind(C, name="fortran_read_int")
  use iso_c_binding
  integer(c_int) :: fortran_read_int
  read(*,*) fortran_read_int
end function fortran_read_int

function fortran_read_real() bind(C, name="fortran_read_real")
  use iso_c_binding
  real(c_double) :: fortran_read_real
  read(*,*) fortran_read_real
end function fortran_read_real

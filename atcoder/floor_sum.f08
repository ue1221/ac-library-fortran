module floor_sum
  use, intrinsic :: iso_fortran_env
  implicit none
  integer, private, parameter :: intkind = int64
contains
  integer(intkind) function floor_linear_sum(n_, m_, a_, b_) result(res)
    integer(intkind), intent(in) :: n_, m_, a_, b_
    integer(intkind) :: n, m, a, b, x, y
    res = 0
    n = n_
    m = m_
    a = a_
    b = b_
    do
      if (a >= m) then
        res = res + (n - 1) * n * (a / m) / 2
        a = mod(a, m)
      end if
      if (b >= m) then
        res = res + n * (b / m)
        b = mod(b, m)
      end if
      y = (a * n + b) / m
      if (y == 0) return
      x = y * m - b
      res = res + (n - (x + a - 1) / a) * y
      n = y
      b = modulo(-x, a)
      y = m
      m = a
      a = y
    end do
  end
end module floor_sum
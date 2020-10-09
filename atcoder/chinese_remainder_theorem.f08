module chinese_remainder_theorem
  use, intrinsic :: iso_fortran_env
  use internal
  implicit none
  integer, private, parameter :: intkind = int64
contains
  function crt(remainders, moduli) result(res)
    integer(intkind), intent(in) :: remainders(:), moduli(:)
    integer(intkind) :: res(2), r0, m0, r1, m1, g, im, tmp(2), u, x
    integer(int32) :: i
    res = 0
    r0 = 0
    m0 = 1
    do i = 1, size(remainders)
      r1 = modulo(remainders(i), moduli(i))
      m1 = moduli(i)
      if (m0 < m1) then
        call swap(r0, r1)
        call swap(m0, m1)
      end if
      if (mod(m0, m1) == 0) then
        if (mod(r0, m1) /= r1) return
        cycle
      end if
      tmp = inv_gcd(m0, m1)
      g = tmp(1)
      im = tmp(2)
      u = m1 / g
      if (mod(r1 - r0, g) /= 0) return
      x = mod(mod((r1 - r0) / g, u) * im, u)
      r0 = r0 + x * m0
      m0 = m0 * u
      if (r0 < 0) r0 = r0 + m0
    end do
    res = (/r0, m0/)
  end
end module chinese_remainder_theorem
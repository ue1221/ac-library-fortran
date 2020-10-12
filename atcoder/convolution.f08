module module_convolution
  use, intrinsic :: iso_fortran_env
  use internal
  implicit none
  integer(int64), private, parameter :: mod1 = 167772161, root1 = 3
  integer(int64), private, parameter :: mod2 = 469762049, root2 = 3
  integer(int64), private, parameter :: mod3 = 1224736769, root3 = 3
  integer(int64), private, parameter :: mod12 = mod1 * mod2
  integer(int64), private, parameter :: inv12 = 104391568 ! inv_mod(mod1, mod2)
  integer(int64), private, parameter :: inv123 = 721017874 ! inv_mod(mod12, mod3)
  private :: pow_mod, swap
contains
  integer(int32) function find_pow2(n) result(res)
    integer(int32), intent(in) :: n
    res = n - 1
    res = or(res, rshift(res, 1))
    res = or(res, rshift(res, 2))
    res = or(res, rshift(res, 4))
    res = or(res, rshift(res, 8))
    res = or(res, rshift(res, 16))
    res = res + 1
  end
  subroutine ntt(n, a, modulus, root, inverse)
    integer(int32), intent(in) :: n
    integer(int64), intent(inout) :: a(0:n - 1)
    integer(int64), intent(in) :: modulus, root
    logical, optional, intent(in) :: inverse
    logical :: inv
    integer(int32) :: i, j, k, l
    integer(int64) :: h, b, w, u, d
    inv = merge(inverse, .false., present(inverse))
    h = pow_mod(root, (modulus - 1) / n, modulus)
    if (inv) h = inv_mod(h, modulus)
    i = 0
    do j = 1, n - 2
      k = rshift(n, 1)
      i = xor(i, k)
      do while (k > i)
        k = rshift(k, 1)
        i = xor(i, k)
      end do
      if (j < i) call swap(a(i), a(j))
    end do
    i = 1
    do while (i < n)
      k = i * 2
      b = pow_mod(h, int(n / k, int64), modulus)
      w = 1
      do j = 0, i - 1
        do l = j, n - 1, k
          u = a(l)
          d = mod(a(l + i) * w, modulus)
          a(l) = u + d
          if (a(l) > modulus) a(l) = a(l) - modulus
          a(l + i) = u - d
          if (a(l + i) < 0) a(l + i) = a(l + i) + modulus
        end do
        w = mod(w * b, modulus)
      end do
      i = i * 2
    end do
    do i = 0, n - 1
      if (a(i) < 0) a(i) = a(i) + modulus
    end do
    if (inv) then
      a = mod(a * inv_mod(int(n, int64), modulus), modulus)
    end if
  end
  function convolution(na_, a, nb_, b, modulus, root) result(res)
    integer(int32), intent(in) :: na_, nb_
    integer(int64), intent(in) :: a(0:na_ - 1), b(0:nb_ - 1), modulus, root
    integer(int32) :: na, nb, n
    integer(int64) :: res(0:na_ + nb_ - 2)
    integer(int64), allocatable :: x(:), y(:), z(:)
    res = 0
    na = na_
    do while (na > 0)
      if (a(na - 1) /= 0) exit
      na = na - 1
    end do
    if (na == 0) return
    nb = nb_
    do while (nb > 0)
      if (b(nb - 1) /= 0) exit
      nb = nb - 1
    end do
    if (nb == 0) return
    n = find_pow2(na + nb - 1)
    allocate(x(0:n - 1), y(0:n - 1), z(0:n - 1))
    x(0:na - 1) = a(0:na - 1)
    x(na:n - 1) = 0
    call ntt(n, x, modulus, root, .false.)
    y(0:nb - 1) = b(0:nb - 1)
    y(nb:n - 1) = 0
    call ntt(n, y, modulus, root, .false.)
    z = mod(x * y, modulus)
    call ntt(n, z, modulus, root, .true.)
    res(0:na + nb - 2) = z(0:na + nb - 2)
  end
  function convolution_ll(na, a, nb, b, modulus) result(res)
    integer(int32), intent(in) :: na, nb
    integer(int64), intent(in) :: a(0:na - 1), b(0:nb - 1), modulus
    integer(int64) :: p(0:na - 1), q(0:nb - 1), modm, v1, v2
    integer(int64), dimension(0:na + nb - 2) :: res, x, y, z
    integer(int32) :: i
    p = modulo(a, modulus)
    q = modulo(b, modulus)
    x = convolution(na, p, nb, q, mod1, root1)
    y = convolution(na, p, nb, q, mod2, root2)
    z = convolution(na, p, nb, q, mod3, root3)
    modm = mod(mod12, modulus)
    do i = 0, na + nb - 2
      v1 = modulo((y(i) - x(i)) * inv12, mod2)
      v2 = modulo((z(i) - mod(x(i) + mod1 * v1, mod3)) * inv123, mod3)
      res(i) = modulo(x(i) + mod1 * v1 + modm * v2, modulus)
    end do
  end
end module
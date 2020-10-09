module internal
  use, intrinsic :: iso_fortran_env
  implicit none
contains

  ! swap `a` and `b`
  subroutine swap(a, b)
    integer(int64), intent(inout) :: a, b
    a = xor(a, b)
    b = xor(a, b)
    a = xor(a, b)
  end

  ! @param n `0 <= n`
  ! @return `i >>> n` in Java
  integer(int32) function urshift(i, n) result(res)
    integer(int32), intent(in) :: i
    integer(int32), intent(in) :: n
    integer(int32) :: m
    m = and(n, 31)
    if (m == 0) then
      res = i
    else
      res = rshift(ibclr(rshift(i, 1), 31), m - 1)
    end if
  end

  ! @param n `0 <= n`
  ! @return minimum non-negative `x` s.t. `n <= 2**x`
  integer(int32) function ceil_pow2(n) result(res)
    integer(int32), intent(in) :: n
    integer(int32) :: i
    i = n
    res = 0
    if (i == 0) return
    res = 1
    if (urshift(i, 16) == 0) then
      res = res + 16
      i = lshift(i, 16)
    end if
    if (urshift(i, 24) == 0) then
      res = res + 8
      i = lshift(i, 8)
    end if
    if (urshift(i, 28) == 0) then
      res = res + 4
      i = lshift(i, 4)
    end if
    if (urshift(i, 30) == 0) then
      res = res + 2
      i = lshift(i, 2)
    end if
    res = 32 - (res - urshift(i, 31))
  end

  ! @param n `0 <= n`
  ! @param m `1 <= m`
  ! @return `(x ** n) % m`
  pure elemental integer(int64) function pow_mod(x, n, m) result(res)
    integer(int64), intent(in) :: x, n, m
    integer(int64) :: y, k
    res = 0
    if (m == 1) return
    res = 1
    y = modulo(x, m)
    k = n
    do while (k > 0)
      if (btest(k, 0)) res = mod(res * y, m)
      y = mod(y * y, m)
      k = rshift(k, 1)
    end do
  end

  ! Reference:
  ! M. Forisek and J. Jancina,
  ! Fast Primality Testing for Integers That Fit into a Machine Word
  ! @param n `0 <= n`
  pure elemental logical function is_prime(n) result(res)
    integer(int64), intent(in) :: n
    integer(int64), parameter :: p(3) = (/2, 7, 61/)
    integer(int64) :: d, a, t, y
    integer(int32) :: i
    res = .false.
    if (n < 2) return
    res = .true.
    if (any(n == p)) return
    res = .false.
    if (.not.btest(n, 0)) return
    d = n - 1
    do while (.not.btest(d, 0))
      d = rshift(d, 1)
    end do
    do i = 1, 3
      a = p(i)
      t = d
      y = pow_mod(a, t, n)
      do while (t /= n - 1 .and. y /= 1 .and. y /= n - 1)
        y = mod(y * y, n)
        t = lshift(t, 1)
      end do
      if (y /= n - 1 .and. .not.btest(t, 0)) return
    end do
    res = .true.
  end

  ! @param b `1 <= b`
  ! @return pair(g, x) s.t. g = gcd(a, b), xa = g (mod b), 0 <= x < b / g
  function inv_gcd(a, b) result(res)
    integer(int64), intent(in) :: a, b
    integer(int64) :: res(2), s, t, u, m0, m1
    t = modulo(a, b)
    res = (/b, 0_8/)
    if (t == 0) return
    s = b
    m0 = 0
    m1 = 1
    do while (t > 0)
      u = s / t
      s = s - t * u
      m0 = m0 - m1 * u
      call swap(s, t)
      call swap(m0, m1)
    end do
    if (m0 < 0) m0 = m0 + b / s
    res = (/s, m0/)
  end

  ! Compile time primitive root
  ! @param m must be prime
  ! @return primitive root (and minimum in now)
  pure elemental integer(int64) function primitive_root(m) result(res)
    integer(int64), intent(in) :: m
    integer(int64) :: divs(20), x, i
    integer(int32) :: cnt
    res = 2
    select case (m)
    case (2)
      res = 1
    case (167772161)
      res = 3
    case (469762049)
      res = 3
    case (754974721)
      res = 11
    case (998244353)
      res = 3
    end select
    if (res /= 2) return
    divs = 0
    divs(1) = 2
    cnt = 1
    x = (m - 1) / 2
    do while (.not.btest(x, 0))
      x = rshift(x, 1)
    end do
    do i = 3, x
      if (i * i > x) exit
      if (mod(x, i) == 0) then
        cnt = cnt + 1
        divs(cnt) = i
        do while (mod(x, i) == 0)
          x = x / i
        end do
      end if
    end do
    if (x > 1) then
      cnt = cnt + 1
      divs(cnt) = x
    end if
    do
      if (all(pow_mod(res, (m - 1) / divs(1:cnt), m) /= 1)) return
      res = res + 1
    end do
  end
end module internal
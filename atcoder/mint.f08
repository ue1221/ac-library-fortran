module module_mint
  use, intrinsic :: iso_fortran_env
  implicit none
  integer(int64), parameter :: modulus = 998244353_int64
  type mint
    integer(int64), private :: val = 0_int64
  contains
    procedure, public :: get => get
    procedure, public :: inv => inv
  end type
  interface mint
    module procedure :: new_mint_int64, new_mint_int32, new_mint_int16, new_mint_int8
  end interface mint
  interface assignment(=)
    module procedure :: asgn_int64, asgn_int32, asgn_int16, asgn_int8
  end interface
  interface operator(+)
    module procedure :: add, add_int64, add_int64_rev, add_int32, add_int32_rev, add_int16, add_int16_rev, add_int8, add_int8_rev
  end interface
  interface operator(-)
    module procedure :: sub, sub_int64, sub_int64_rev, sub_int32, sub_int32_rev, sub_int16, sub_int16_rev, sub_int8, sub_int8_rev
  end interface
  interface operator(*)
    module procedure :: mul, mul_int64, mul_int64_rev, mul_int32, mul_int32_rev, mul_int16, mul_int16_rev, mul_int8, mul_int8_rev
  end interface
  interface operator(/)
    module procedure :: div, div_int64, div_int64_rev, div_int32, div_int32_rev, div_int16, div_int16_rev, div_int8, div_int8_rev
  end interface
  interface operator(**)
    module procedure :: pow_int64, pow_int32, pow_int16, pow_int8
  end interface
  interface dot_product
    module procedure :: dot_prod
  end interface dot_product
  interface matmul
    module procedure :: matmul1, matmul2, matmul3
  end interface matmul
contains
  pure elemental integer(int64) function get(this) result(res)
    class(mint), intent(in) :: this
    res = this%val
  end
  pure elemental subroutine asgn_int64(this, other)
    type(mint), intent(inout) :: this
    integer(int64), intent(in) :: other
    this%val = merge(other, modulo(other, modulus), 0_int64 <= other .and. other < modulus)
  end
  pure elemental subroutine asgn_int32(this, other)
    type(mint), intent(inout) :: this
    integer(int32), intent(in) :: other
    call asgn_int64(this, int(other, int64))
  end
  pure elemental subroutine asgn_int16(this, other)
    type(mint), intent(inout) :: this
    integer(int16), intent(in) :: other
    call asgn_int64(this, int(other, int64))
  end
  pure elemental subroutine asgn_int8(this, other)
    type(mint), intent(inout) :: this
    integer(int8), intent(in) :: other
    call asgn_int64(this, int(other, int64))
  end
  pure elemental type(mint) function add(this, other) result(res)
    type(mint), intent(in) :: this
    type(mint), intent(in) :: other
    integer(int64) :: tmp
    tmp = this%val + other%val
    res%val = merge(tmp, tmp - modulus, 0_int64 <= tmp .and. tmp < modulus)
  end
  pure elemental type(mint) function sub(this, other) result(res)
    type(mint), intent(in) :: this
    type(mint), intent(in) :: other
    integer(int64) :: tmp
    tmp = this%val - other%val
    res%val = merge(tmp, tmp + modulus, 0_int64 <= tmp .and. tmp < modulus)
  end
  pure elemental type(mint) function mul(this, other) result(res)
    type(mint), intent(in) :: this
    type(mint), intent(in) :: other
    integer(int64) :: tmp
    tmp = this%val * other%val
    res%val = merge(tmp, modulo(tmp, modulus), 0_int64 <= tmp .and. tmp < modulus)
  end
  pure elemental type(mint) function div(this, other) result(res)
    type(mint), intent(in) :: this
    type(mint), intent(in) :: other
    res = mul(this, inv(other))
  end
  pure elemental type(mint) function inv(this) result(res)
    class(mint), intent(in) :: this
    integer(int64) :: a, b, c, n
    integer(int64) :: s, t
    a = this%val
    b = modulus
    c = 0_int64
    n = 1_int64
    do while (b /= 0_int64)
      s = mod(a, b)
      t = n - a / b * c
      a = b
      b = s
      n = c
      c = t
    end do
    res%val = modulo(n, modulus)
  end
  pure elemental type(mint) function pow_int64(this, other) result(res)
    type(mint), intent(in) :: this
    integer(int64), intent(in) :: other
    integer(int64) :: a, k, ak
    a = get(this)
    k = other
    ak = 1_int64
    do while (k > 0_int64)
      if (btest(k, 0)) then
        ak = ak * a
        if (ak >= modulus) ak = mod(ak, modulus)
      end if
      a = a * a
      if (a >= modulus) a = mod(a, modulus)
      k = rshift(k, 1)
    end do
    res%val = modulo(ak, modulus)
  end
  pure elemental type(mint) function pow_int32(this, other) result(res)
    type(mint), intent(in) :: this
    integer(int32), intent(in) :: other
    res = pow_int64(this, int(other, int64))
  end
  pure elemental type(mint) function pow_int16(this, other) result(res)
    type(mint), intent(in) :: this
    integer(int16), intent(in) :: other
    res = pow_int64(this, int(other, int64))
  end
  pure elemental type(mint) function pow_int8(this, other) result(res)
    type(mint), intent(in) :: this
    integer(int8), intent(in) :: other
    res = pow_int64(this, int(other, int64))
  end
  pure type(mint) function dot_prod(x, y) result(res)
    type(mint), intent(in) :: x(:), y(:)
    integer(int32) :: l, k
    l = size(x, 1)
    do k = 1, l
      res%val = mod(res%val + x(k)%val * y(k)%val, modulus)
    end do
  end
  pure function matmul1(x, y) result(res)
    type(mint), intent(in) :: x(:, :), y(:)
    type(mint) :: res(size(x, 1))
    integer(int32) :: n, l, i, k
    integer(int64) :: tmp
    n = size(x, 1)
    l = size(x, 2)
    do k = 1, l
      tmp = y(k)%val
      do i = 1, n
        res(i)%val = mod(res(i)%val + x(i, k)%val * tmp, modulus)
      end do
    end do
  end
  pure function matmul2(x, y) result(res)
    type(mint), intent(in) :: x(:), y(:, :)
    type(mint) :: res(size(y, 2))
    integer(int32) :: m, l, j, k
    m = size(y, 2)
    l = size(y, 1)
    do j = 1, m
      do k = 1, l
        res(j)%val = mod(res(j)%val + x(k)%val * y(k, j)%val, modulus)
      end do
    end do
  end
  pure function matmul3(x, y) result(res)
    type(mint), intent(in) :: x(:, :), y(:, :)
    type(mint) :: res(size(x, 1), size(y, 2))
    integer(int32) :: n, m, l, i, j, k
    integer(int64) :: tmp
    n = size(x, 1)
    m = size(y, 2)
    l = size(x, 2)
    do j = 1, m
      do k = 1, l
        tmp = y(k, j)%val
        do i = 1, n
          res(i, j)%val = mod(res(i, j)%val + x(i, k)%val * tmp, modulus)
        end do
      end do
    end do
  end
  pure function matpow(x, n_) result(res)
    type(mint), intent(in) :: x(:, :)
    integer(int64), intent(in) :: n_
    type(mint) :: res(size(x, 1), size(x, 2)), tmp(size(x, 1), size(x, 2))
    integer(int32) :: i, n
    do i = 1, size(x, 1)
      res(i, i)%val = 1_int64
    end do
    tmp = x
    n = n_
    do while (n > 0_int64)
      if (btest(n, 0)) res = matmul3(res, tmp)
      tmp = matmul3(tmp, tmp)
      n = rshift(n, 1)
    end do
  end
  pure elemental type(mint) function new_mint_int64(val) result(res)
    integer(int64), intent(in) :: val
    call asgn_int64(res, val)
  end
  pure elemental type(mint) function add_int64(this, other) result(res)
    type(mint), intent(in) :: this
    integer(int64), intent(in) :: other
    res = add(this, new_mint_int64(other))
  end
  pure elemental type(mint) function sub_int64(this, other) result(res)
    type(mint), intent(in) :: this
    integer(int64), intent(in) :: other
    res = sub(this, new_mint_int64(other))
  end
  pure elemental type(mint) function mul_int64(this, other) result(res)
    type(mint), intent(in) :: this
    integer(int64), intent(in) :: other
    res = mul(this, new_mint_int64(other))
  end
  pure elemental type(mint) function div_int64(this, other) result(res)
    type(mint), intent(in) :: this
    integer(int64), intent(in) :: other
    res = div(this, new_mint_int64(other))
  end
  pure elemental type(mint) function add_int64_rev(other, this) result(res)
    type(mint), intent(in) :: this
    integer(int64), intent(in) :: other
    res = add(new_mint_int64(other), this)
  end
  pure elemental type(mint) function sub_int64_rev(other, this) result(res)
    type(mint), intent(in) :: this
    integer(int64), intent(in) :: other
    res = sub(new_mint_int64(other), this)
  end
  pure elemental type(mint) function mul_int64_rev(other, this) result(res)
    type(mint), intent(in) :: this
    integer(int64), intent(in) :: other
    res = mul(new_mint_int64(other), this)
  end
  pure elemental type(mint) function div_int64_rev(other, this) result(res)
    type(mint), intent(in) :: this
    integer(int64), intent(in) :: other
    res = div(new_mint_int64(other), this)
  end
  pure elemental type(mint) function new_mint_int32(val) result(res)
    integer(int32), intent(in) :: val
    call asgn_int32(res, val)
  end
  pure elemental type(mint) function add_int32(this, other) result(res)
    type(mint), intent(in) :: this
    integer(int32), intent(in) :: other
    res = add(this, new_mint_int32(other))
  end
  pure elemental type(mint) function sub_int32(this, other) result(res)
    type(mint), intent(in) :: this
    integer(int32), intent(in) :: other
    res = sub(this, new_mint_int32(other))
  end
  pure elemental type(mint) function mul_int32(this, other) result(res)
    type(mint), intent(in) :: this
    integer(int32), intent(in) :: other
    res = mul(this, new_mint_int32(other))
  end
  pure elemental type(mint) function div_int32(this, other) result(res)
    type(mint), intent(in) :: this
    integer(int32), intent(in) :: other
    res = div(this, new_mint_int32(other))
  end
  pure elemental type(mint) function add_int32_rev(other, this) result(res)
    type(mint), intent(in) :: this
    integer(int32), intent(in) :: other
    res = add(new_mint_int32(other), this)
  end
  pure elemental type(mint) function sub_int32_rev(other, this) result(res)
    type(mint), intent(in) :: this
    integer(int32), intent(in) :: other
    res = sub(new_mint_int32(other), this)
  end
  pure elemental type(mint) function mul_int32_rev(other, this) result(res)
    type(mint), intent(in) :: this
    integer(int32), intent(in) :: other
    res = mul(new_mint_int32(other), this)
  end
  pure elemental type(mint) function div_int32_rev(other, this) result(res)
    type(mint), intent(in) :: this
    integer(int32), intent(in) :: other
    res = div(new_mint_int32(other), this)
  end
  pure elemental type(mint) function new_mint_int16(val) result(res)
    integer(int16), intent(in) :: val
    call asgn_int16(res, val)
  end
  pure elemental type(mint) function add_int16(this, other) result(res)
    type(mint), intent(in) :: this
    integer(int16), intent(in) :: other
    res = add(this, new_mint_int16(other))
  end
  pure elemental type(mint) function sub_int16(this, other) result(res)
    type(mint), intent(in) :: this
    integer(int16), intent(in) :: other
    res = sub(this, new_mint_int16(other))
  end
  pure elemental type(mint) function mul_int16(this, other) result(res)
    type(mint), intent(in) :: this
    integer(int16), intent(in) :: other
    res = mul(this, new_mint_int16(other))
  end
  pure elemental type(mint) function div_int16(this, other) result(res)
    type(mint), intent(in) :: this
    integer(int16), intent(in) :: other
    res = div(this, new_mint_int16(other))
  end
  pure elemental type(mint) function add_int16_rev(other, this) result(res)
    type(mint), intent(in) :: this
    integer(int16), intent(in) :: other
    res = add(new_mint_int16(other), this)
  end
  pure elemental type(mint) function sub_int16_rev(other, this) result(res)
    type(mint), intent(in) :: this
    integer(int16), intent(in) :: other
    res = sub(new_mint_int16(other), this)
  end
  pure elemental type(mint) function mul_int16_rev(other, this) result(res)
    type(mint), intent(in) :: this
    integer(int16), intent(in) :: other
    res = mul(new_mint_int16(other), this)
  end
  pure elemental type(mint) function div_int16_rev(other, this) result(res)
    type(mint), intent(in) :: this
    integer(int16), intent(in) :: other
    res = div(new_mint_int16(other), this)
  end
  pure elemental type(mint) function new_mint_int8(val) result(res)
    integer(int8), intent(in) :: val
    call asgn_int8(res, val)
  end
  pure elemental type(mint) function add_int8(this, other) result(res)
    type(mint), intent(in) :: this
    integer(int8), intent(in) :: other
    res = add(this, new_mint_int8(other))
  end
  pure elemental type(mint) function sub_int8(this, other) result(res)
    type(mint), intent(in) :: this
    integer(int8), intent(in) :: other
    res = sub(this, new_mint_int8(other))
  end
  pure elemental type(mint) function mul_int8(this, other) result(res)
    type(mint), intent(in) :: this
    integer(int8), intent(in) :: other
    res = mul(this, new_mint_int8(other))
  end
  pure elemental type(mint) function div_int8(this, other) result(res)
    type(mint), intent(in) :: this
    integer(int8), intent(in) :: other
    res = div(this, new_mint_int8(other))
  end
  pure elemental type(mint) function add_int8_rev(other, this) result(res)
    type(mint), intent(in) :: this
    integer(int8), intent(in) :: other
    res = add(new_mint_int8(other), this)
  end
  pure elemental type(mint) function sub_int8_rev(other, this) result(res)
    type(mint), intent(in) :: this
    integer(int8), intent(in) :: other
    res = sub(new_mint_int8(other), this)
  end
  pure elemental type(mint) function mul_int8_rev(other, this) result(res)
    type(mint), intent(in) :: this
    integer(int8), intent(in) :: other
    res = mul(new_mint_int8(other), this)
  end
  pure elemental type(mint) function div_int8_rev(other, this) result(res)
    type(mint), intent(in) :: this
    integer(int8), intent(in) :: other
    res = div(new_mint_int8(other), this)
  end
end module
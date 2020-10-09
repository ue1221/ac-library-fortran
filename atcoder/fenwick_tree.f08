module mod_fenwick_tree
  use, intrinsic :: iso_fortran_env
  implicit none
  integer, private, parameter :: intkind = int32
  type fenwick_tree
    private
    integer(int32) :: n, p
    integer(intkind), allocatable :: body(:)
  contains
    procedure :: init => init
    procedure :: add => add
    procedure :: sum0 => sum0
    procedure :: sum1 => sum1
    generic :: sum => sum0, sum1
    procedure :: get => get
    procedure :: lower_bound => lower_bound
    procedure :: update => update
    procedure :: max => max0
  end type
  interface fenwick_tree
    module procedure :: newft0
  end interface
contains
  type(fenwick_tree) function newft0(n) result(res)
    integer(int32), intent(in) :: n
    integer(int32) :: p
    p = 1
    do while (lshift(p, 1) <= n)
      p = lshift(p, 1)
    end do
    res%n = n
    res%p = p
    allocate(res%body(n))
    res%body = 0
  end
  subroutine init(this, n)
    class(fenwick_tree), intent(out) :: this
    integer(int32), intent(in) :: n
    integer(int32) :: p
    p = 1
    do while (lshift(p, 1) <= n)
      p = lshift(p, 1)
    end do
    this%n = n
    this%p = p
    allocate(this%body(n))
    this%body = 0
  end
  subroutine add(this, i, v)
    class(fenwick_tree), intent(inout) :: this
    integer(int32), intent(in) :: i
    integer(intkind), intent(in) :: v
    integer(int32) :: x
    x = i
    do while (x <= this%n)
      this%body(x) = this%body(x) + v
      x = x + and(x, -x)
    end do
  end
  integer(intkind) function sum0(this, i) result(res)
    class(fenwick_tree), intent(in) :: this
    integer(int32), intent(in) :: i
    integer(int32) :: x
    res = 0
    x = i
    do while (x > 0)
      res = res + this%body(x)
      x = x - and(x, -x)
    end do
  end
  integer(intkind) function sum1(this, l, r) result(res)
    class(fenwick_tree), intent(in) :: this
    integer(int32), intent(in) :: l, r
    integer(int32) :: i, j
    res = 0
    i = l - 1
    j = r
    do while (j > i)
      res = res + this%body(j)
      j = j - and(j, -j)
    end do
    do while (i > j)
      res = res - this%body(i)
      i = i - and(i, -i)
    end do
  end
  integer(intkind) function get(this, i) result(res)
    class(fenwick_tree), intent(in) :: this
    integer(int32), intent(in) :: i
    integer(int32) :: x
    res = sum1(this, i, i)
  end
  integer(int32) function lower_bound(this, v) result(res)
    class(fenwick_tree), intent(in) :: this
    integer(intkind), intent(in) :: v
    integer(int32) :: k
    integer(intkind) :: w
    res = 0
    if (v <= 0) return
    k = this%p
    w = v
    do while (k > 0)
      if (res + k <= this%n .and. this%body(res + k) < w) then
        w = w - this%body(res + k)
        res = res + k
      end if
      k = rshift(k, 1)
    end do
    res = res + 1
  end
  subroutine update(this, i, v)
    class(fenwick_tree), intent(inout) :: this
    integer(int32), intent(in) :: i
    integer(intkind), intent(in) :: v
    integer(int32) :: x
    x = i
    do while (x <= this%n)
      if (this%body(x) < v) this%body(x) = v
      x = x + and(x, -x)
    end do
  end
  integer(intkind) function max0(this, i) result(res)
    class(fenwick_tree), intent(in) :: this
    integer(int32), intent(in) :: i
    integer(int32) :: x
    x = i
    res = 0
    do while (x > 0)
      res = max(res, this%body(x))
      x = x - and(x, -x)
    end do
  end
  integer(intkind) function lis_length(a) result(res)
    integer, intent(in) :: a(:)
    type(fenwick_tree) :: bit
    integer(intkind) :: tmp
    integer(int32) :: n, i
    n = size(a)
    call init(bit, n)
    res = 0
    do i = 1, n
      tmp = max0(bit, a(i) - 1) + 1
      res = max(res, tmp)
      call update(bit, a(i), tmp)
    end do
  end
  integer(intkind) function inversion_num(a) result(res)
    integer(int32), intent(in) :: a(:)
    type(fenwick_tree) :: bit
    integer(int32) :: n, i
    n = size(a)
    call init(bit, n)
    res = 0
    do i = 1, n
      call add(bit, a(i), 1_intkind)
      res = res + i - sum0(bit, a(i))
    end do
  end
end module mod_fenwick_tree
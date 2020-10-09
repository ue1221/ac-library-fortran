module module_segtree
  use, intrinsic :: iso_fortran_env
  implicit none
  type monoid
    ! Define by yourself when you use this.
  end type
  type segtree
    procedure(operation), pointer, nopass :: op => null()
    procedure(identity), pointer, nopass :: e => null()
    type(monoid), allocatable :: array(:)
    integer(int32) :: n, p
  contains
    procedure :: set => set
    procedure :: get => get
    procedure :: prod => prod
    procedure :: all_prod => all_prod
    procedure :: max_right => max_right
    procedure :: min_left => min_left
  end type
  abstract interface
    type(monoid) function operation(a, b) result(res)
      import monoid
      type(monoid), intent(in) :: a, b
    end
  end interface
  abstract interface
    type(monoid) function identity() result(res)
      import monoid
    end
  end interface
  interface segtree
    module procedure :: segtree0, segtree1, segtree2
  end interface
  private :: update
contains
  type(segtree) function segtree0(op, e) result(res)
    procedure(operation) :: op
    procedure(identity) :: e
    res = segtree1(1, op, e)
  end
  type(segtree) function segtree1(n, op, e) result(res)
    integer(int32), intent(in) :: n
    procedure(operation) :: op
    procedure(identity) :: e
    type(monoid), allocatable :: array(:)
    integer(int32) :: i
    allocate(array(n))
    do i = 1, n
      array(i) = e()
    end do
    res = segtree2(array, op, e)
  end
  type(segtree) function segtree2(array, op, e) result(res)
    type(monoid), intent(in) :: array(:)
    procedure(operation) :: op
    procedure(identity) :: e
    integer(int32) :: n, p, i
    n = size(array)
    p = 1
    do while (p < n)
      p = lshift(p, 1)
    end do
    res%e => e
    res%op => op
    allocate(res%array(2 * p - 1))
    res%n = n
    res%p = p
    do i = p, n + 1, -1
      res%array(i + p - 1) = e()
    end do
    do i = n, 1, -1
      res%array(i + p - 1) = array(i)
    end do
    do i = p - 1, 1, -1
      call update(res, i)
    end do
  end
  subroutine update(this, i)
    type(segtree), intent(inout) :: this
    integer(int32), intent(in) :: i
    this%array(i) = this%op(this%array(2 * i), this%array(2 * i + 1))
  end
  subroutine set(this, i, x)
    class(segtree), intent(inout) :: this
    integer(int32), intent(in) :: i
    type(monoid), intent(in) :: x
    integer(int32) :: k
    k = i + this%p - 1
    this%array(k) = x
    do while (k > 1)
      k = k / 2
      call update(this, k)
    end do
  end
  type(monoid) function get(this, i) result(res)
    class(segtree), intent(in) :: this
    integer(int32), intent(in) :: i
    res = this%array(i + this%p - 1)
  end
  type(monoid) function prod(this, a, b) result(res)
    class(segtree), intent(in) :: this
    integer(int32), intent(in) :: a, b
    type(monoid) :: x, y
    integer(int32) :: l, r
    x = this%e()
    y = this%e()
    l = a + this%p - 1
    r = b + this%p
    do while (l < r)
      if (btest(l, 0)) then
        x = this%op(x, this%array(l))
        l = l + 1
      end if
      if (btest(r, 0)) then
        r = r - 1
        y = this%op(this%array(r), y)
      end if
      l = rshift(l, 1)
      r = rshift(r, 1)
    end do
    res = this%op(x, y)
  end
  type(monoid) function all_prod(this) result(res)
    class(segtree), intent(in) :: this
    res = this%array(1)
  end
  integer(int32) function max_right(this, l, f) result(res)
    class(segtree), intent(in) :: this
    integer(int32), intent(in) :: l
    type(monoid) :: x, y
    interface
      logical function f(x) result(res)
        import monoid
        type(monoid), intent(in) :: x
      end
    end interface
    res = this%n + 1
    if (l > this%n) return
    x = this%e()
    res = l + this%p - 1
    do
      do while (.not.btest(res, 0))
        res = rshift(res, 1)
      end do
      if (.not.f(this%op(x, this%array(res)))) then
        do while (res < this%p)
          res = lshift(res, 1)
          y = this%op(x, this%array(res))
          if (f(y)) then
            x = y
            res = res + 1
          end if
        end do
        res = res - this%p + 1
        return
      end if
      x = this%op(x, this%array(res))
      res = res + 1
      if (and(res, -res) == res) exit
    end do
    res = this%n + 1
  end
  integer(int32) function min_left(this, r, f) result(res)
    class(segtree), intent(in) :: this
    integer(int32), intent(in) :: r
    type(monoid) :: x, y
    interface
      logical function f(x) result(res)
        import monoid
        type(monoid), intent(in) :: x
      end
    end interface
    res = 0
    if (r < 1) return
    x = this%e()
    res = r + this%p
    do
      res = res - 1
      do while (res > 1 .and. btest(res, 0))
        res = rshift(res, 1)
      end do
      if (.not.f(this%op(this%array(res), x))) then
        do while (res < this%p)
          res = lshift(res, 1) + 1
          y = this%op(this%array(res), x)
          if (f(y)) then
            x = y
            res = res - 1
          end if
        end do
        res = res - this%p
        return
      end if
      x = this%op(this%array(res), x)
      if (and(res, -res) == res) exit
    end do
    res = 0
  end
end module module_segtree
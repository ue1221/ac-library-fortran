module module_lazy_segtree
  use, intrinsic :: iso_fortran_env
  implicit none
  type monoid
    ! Define by yourself when you use this.
  end type
  type func
    ! Define by yourself when you use this.
  end type
  type lazy_segtree
    procedure(operation), pointer, nopass :: op => null()
    procedure(identity_monoid), pointer, nopass :: e => null()
    procedure(mapping), pointer, nopass :: map => null()
    procedure(composition), pointer, nopass :: compose => null()
    procedure(identity_func), pointer, nopass :: id => null()
    type(monoid), allocatable :: array(:)
    type(func), allocatable :: lazy(:)
    integer(int32) :: n, logp, p
  contains
    procedure :: set => set
    procedure :: get => get
    procedure :: prod => prod
    procedure :: all_prod => all_prod
    procedure, private :: apply_at => apply_at
    procedure, private :: apply_range => apply_range
    generic :: apply => apply_at, apply_range
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
    type(monoid) function identity_monoid() result(res)
      import monoid
    end
  end interface
  abstract interface
    type(monoid) function mapping(f, s) result(res)
      import monoid, func
      type(func), intent(in) :: f
      type(monoid), intent(in) :: s
    end
  end interface
  abstract interface
    type(func) function composition(f, g) result(res)
      import func
      type(func), intent(in) :: f, g
    end
  end interface
  abstract interface
    type(func) function identity_func() result(res)
      import func
    end
  end interface
  interface lazy_segtree
    module procedure :: lazy_segtree0, lazy_segtree1, lazy_segtree2
  end interface
  private :: update, all_apply, push
contains
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
  type(lazy_segtree) function lazy_segtree0(op, e, map, compose, id) result(res)
    procedure(operation) :: op
    procedure(identity_monoid) :: e
    procedure(mapping) :: map
    procedure(composition) :: compose
    procedure(identity_func) :: id
    res = lazy_segtree1(1, op, e, map, compose, id)
  end
  type(lazy_segtree) function lazy_segtree1(n, op, e, map, compose, id) result(res)
    integer(int32), intent(in) :: n
    procedure(operation) :: op
    procedure(identity_monoid) :: e
    procedure(mapping) :: map
    procedure(composition) :: compose
    procedure(identity_func) :: id
    type(monoid), allocatable :: array(:)
    integer(int32) :: i
    allocate(array(n))
    do i = 1, n
      array(i) = e()
    end do
    res = lazy_segtree2(array, op, e, map, compose, id)
  end
  type(lazy_segtree) function lazy_segtree2(array, op, e, map, compose, id) result(res)
    type(monoid), intent(in) :: array(:)
    procedure(operation) :: op
    procedure(identity_monoid) :: e
    procedure(mapping) :: map
    procedure(composition) :: compose
    procedure(identity_func) :: id
    integer(int32) :: n, logp, p, i
    n = size(array)
    logp = ceil_pow2(n)
    p = lshift(1, logp)
    res%op => op
    res%e => e
    res%map => map
    res%compose => compose
    res%id => id
    allocate(res%array(2 * p - 1), res%lazy(p))
    res%n = n
    res%logp = logp
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
    do i = 1, p
      res%lazy(i) = id()
    end do
  end
  subroutine update(this, i)
    type(lazy_segtree), intent(inout) :: this
    integer(int32), intent(in) :: i
    this%array(i) = this%op(this%array(2 * i), this%array(2 * i + 1))
  end
  subroutine all_apply(this, i, f)
    type(lazy_segtree), intent(inout) :: this
    integer(int32), intent(in) :: i
    type(func), intent(in) :: f
    this%array(i) = this%map(f, this%array(i))
    if (i < this%p) this%lazy(i) = this%compose(f, this%lazy(i))
  end
  subroutine push(this, i)
    type(lazy_segtree), intent(inout) :: this
    integer(int32), intent(in) :: i
    call all_apply(this, 2 * i, this%lazy(i))
    call all_apply(this, 2 * i + 1, this%lazy(i))
    this%lazy(i) = this%id()
  end
  subroutine set(this, i, x)
    class(lazy_segtree), intent(inout) :: this
    integer(int32), intent(in) :: i
    type(monoid) :: x
    integer(int32) :: p, k
    p = i + this%p - 1
    do k = this%logp, 1, -1
      call push(this, rshift(p, k))
    end do
    this%array(p) = x
    do k = 1, this%logp
      call update(this, rshift(p, k))
    end do
  end
  type(monoid) function get(this, i) result(res)
    class(lazy_segtree), intent(inout) :: this
    integer(int32), intent(in) :: i
    integer(int32) :: p, k
    p = i + this%p - 1
    do k = this%logp, 1, -1
      call push(this, rshift(p, k))
    end do
    res = this%array(p)
  end
  type(monoid) function prod(this, a, b) result(res)
    class(lazy_segtree), intent(inout) :: this
    integer(int32), intent(in) :: a, b
    type(monoid) :: x, y
    integer(int32) :: l, r, i
    l = a + this%p - 1
    r = b + this%p
    if (l >= r) then
      res = this%e()
      return
    end if
    do i = this%logp, 1, -1
      if (lshift(rshift(l, i), i) /= l) call push(this, rshift(l, i))
      if (lshift(rshift(r, i), i) /= r) call push(this, rshift(r, i))
    end do
    x = this%e()
    y = this%e()
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
    class(lazy_segtree), intent(in) :: this
    res = this%array(1)
  end
  subroutine apply_at(this, i, f)
    class(lazy_segtree), intent(inout) :: this
    integer(int32), intent(in) :: i
    type(func), intent(in) :: f
    integer(int32) :: p, k
    p = i + this%p - 1
    do k = this%logp, 1, -1
      call push(this, rshift(p, k))
    end do
    this%array(p) = this%map(f, this%array(p))
    do k = 1, this%logp
      call update(this, rshift(p, k))
    end do
  end
  subroutine apply_range(this, a, b, f)
    class(lazy_segtree), intent(inout) :: this
    integer(int32), intent(in) :: a, b
    type(func), intent(in) :: f
    integer(int32) :: l, r, i
    l = a + this%p - 1
    r = b + this%p
    if (l >= r) return
    do i = this%logp, 1, -1
      if (lshift(rshift(l, i), i) /= l) call push(this, rshift(l, i))
      if (lshift(rshift(r, i), i) /= r) call push(this, rshift(r - 1, i))
    end do
    block
      integer(int32) :: l2, r2
      l2 = l
      r2 = r
      do while (l < r)
        if (btest(l, 0)) then
          call all_apply(this, l, f) 
          l = l + 1
        end if
        if (btest(r, 0)) then
          r = r - 1
          call all_apply(this, r, f) 
        end if
        l = rshift(l, 1)
        r = rshift(r, 1)
      end do
      l = l2
      r = r2
    end block
    do i = 1, this%logp
      if (lshift(rshift(l, i), i) /= l) call update(this, rshift(l, i))
      if (lshift(rshift(r, i), i) /= r) call update(this, rshift(r - 1, i))
    end do
  end
  integer(int32) function max_right(this, l, f) result(res)
    class(lazy_segtree), intent(inout) :: this
    integer(int32), intent(in) :: l
    type(monoid) :: x, y
    integer(int32) :: i
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
    do i = this%logp, 1, -1
      call push(this, rshift(res, i))
    end do
    do
      do while (.not.btest(res, 0))
        res = rshift(res, 1)
      end do
      if (.not.f(this%op(x, this%array(res)))) then
        do while (res < this%p)
          call push(this, res)
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
    class(lazy_segtree), intent(inout) :: this
    integer(int32), intent(in) :: r
    type(monoid) :: x, y
    integer(int32) :: i
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
    do i = this%logp, 1, -1
      call push(this, rshift(res - 1, i))
    end do
    do
      res = res - 1
      do while (res > 1 .and. btest(res, 0))
        res = rshift(res, 1)
      end do
      if (.not.f(this%op(this%array(res), x))) then
        do while (res < this%p)
          call push(this, res)
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
end module module_lazy_segtree
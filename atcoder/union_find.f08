module mod_union_find
  use, intrinsic :: iso_fortran_env
  implicit none
  type union_find
    integer(int32) :: n = 0
    integer(int32), allocatable :: par(:), rnk(:), cnt(:)
  contains
    procedure :: find => find
    procedure :: same => same
    procedure :: unite => unite
    procedure :: size => size_of
  end type
  interface union_find
    module procedure :: newuf
  end interface
contains
  type(union_find) function newuf(n) result(res)
    integer(int32), intent(in) :: n
    integer(int32) :: i
    res%n = n
    allocate(res%par(n), res%rnk(n), res%cnt(n))
    res%par = [(i, i = 1, n)]
    res%rnk = 0
    res%cnt = 1
  end
  recursive integer function find(this, i) result(res)
    class(union_find), intent(inout) :: this
    integer(int32), intent(in) :: i
    if (this%par(i) /= i) this%par(i) = find(this, this%par(i))
    res = this%par(i)
  end
  logical function same(this, i, j) result(res)
    class(union_find), intent(inout) :: this
    integer(int32), intent(in) :: i, j
    res = find(this, i) == find(this, j)
  end
  subroutine unite(this, i, j)
    class(union_find), intent(inout) :: this
    integer(int32), intent(in) :: i, j
    integer(int32) :: x, y
    x = find(this, i)
    y = find(this, j)
    if (x == y) return
    if (this%rnk(x) < this%rnk(y)) then
      this%par(x) = y
      this%cnt(y) = this%cnt(y) + this%cnt(x)
    else
      this%par(y) = x
      this%cnt(x) = this%cnt(x) + this%cnt(y)
      if (this%rnk(x) == this%rnk(y)) this%rnk(x) = this%rnk(x) + 1
    end if
  end
  integer(int32) function size_of(this, i) result(res)
    class(union_find), intent(inout) :: this
    integer(int32), intent(in) :: i
    res = this%cnt(find(this, i))
  end
end module mod_union_find
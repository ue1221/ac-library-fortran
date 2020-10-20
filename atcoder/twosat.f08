module module_twosat
  use, intrinsic :: iso_fortran_env
  use strongly_connected_component
  implicit none
  type twosat
    integer(int32), private :: n
    logical, private, allocatable :: ans(:)
    type(scc_graph), private :: scc
  contains
    procedure :: add_clause => add_clause
    procedure :: satisfiable => satisfiable
    procedure :: answer => answer
  end type
  interface twosat
    module procedure :: new_twosat
  end interface
contains
  type(twosat) function new_twosat(n) result(res)
    integer(int32), intent(in) :: n
    res%n = n
    allocate(res%ans(n))
    res%ans = .false.
    res%scc = new_scc_graph(2 * n)
  end
  subroutine add_clause(this, i, f, j, g)
    class(twosat), intent(inout) :: this
    integer(int32), intent(in) :: i, j
    logical, intent(in) :: f, g
    call add_edge(this%scc, &
      2 * i - 1 + merge(0, 1, f), 2 * j - 1 + merge(1, 0, g))
    call add_edge(this%scc, &
      2 * j - 1 + merge(0, 1, g), 2 * i - 1 + merge(1, 0, f))
  end
  logical function satisfiable(this) result(res)
    class(twosat), intent(inout) :: this
    type(array_list) :: id
    integer(int32) :: i
    res = .false.
    id = scc_ids(this%scc)
    do i = 1, this%n
      if (id%array(2 * i - 1) == id%array(2 * i)) return
      this%ans(i) = id%array(2 * i - 1) < id%array(2 * i)
    end do
    res = .true.
  end
  function answer(this) result(res)
    class(twosat), intent(in) :: this
    logical :: res(this%n)
    res = this%ans
  end
end module
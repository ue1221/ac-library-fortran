module strongly_connected_component
  use, intrinsic :: iso_fortran_env
  implicit none
  type, private :: node
    integer(int32) :: from, to
    type(node), pointer :: prev => null()
    type(node), pointer :: next => null()
  end type
  type, private :: linked_list
    integer(int32) :: size = 0
    type(node), pointer :: head => null()
    type(node), pointer :: tail => null()
  end type
  type array_list
    integer(int32) :: size = 0
    integer, allocatable :: array(:)
  end type
  type vertex_group
    type(array_list), allocatable :: group(:)
  end type
  type csr !compressed sparse row
    integer(int32), allocatable :: start(:)
    integer(int32), allocatable :: elist(:)
  end type
  type scc_graph
    integer(int32), private :: n
    type(linked_list), private :: edges
  contains
    procedure, public :: num_vertices => num_vertices
    procedure, public :: add_edge => add_edge
    procedure, public :: scc_ids => scc_ids
    procedure, public :: scc => scc
  end type
  interface scc_graph
    module procedure :: new_scc_graph
  end interface
  private :: new_node
contains
  function new_node(from, to) result(res)
    integer(int32), intent(in) :: from, to
    type(node), pointer :: res
    allocate(res)
    res%from = from
    res%to = to
  end
  subroutine add_linked(this, from, to)
    type(linked_list), intent(inout) :: this
    integer(int32), intent(in) :: from, to
    type(node), pointer :: now
    now => new_node(from, to)
    if (associated(this%head)) then
      now%prev => this%tail
      this%tail%next => now
    else
      this%head => now
    end if
    this%tail => now
    this%size = this%size + 1
  end
  function pop_linked(this) result(res)
    type(linked_list), intent(inout) :: this
    type(node), pointer :: res, now
    res => this%tail
    now => res%prev
    res%prev => null()
    this%tail => now
    if (associated(now)) then
      now%next => null()
    else
      this%head => null()
    end if
    this%size = this%size - 1
  end
  subroutine add_array(this, val)
    type(array_list), intent(inout) :: this
    integer(int32), intent(in) :: val
    if (.not.allocated(this%array)) allocate(this%array(1))
    if (size(this%array) == this%size) call append_array(this, this%size)
    this%size = this%size + 1
    this%array(this%size) = val
  end
  subroutine append_array(this, size)
    type(array_list), intent(inout) :: this
    integer(int32), intent(in) :: size
    integer(int32) :: array(size)
    array = this%array
    deallocate(this%array)
    allocate(this%array(2 * size))
    this%array(1:size) = array
  end
  integer(int32) function pop_array(this) result(res)
    type(array_list), intent(inout) :: this
    res = -1
    if (this%size == 0 .or. .not.allocated(this%array)) return
    res = this%array(this%size)
    this%size = this%size - 1
  end
  type(csr) function new_csr(n, edges) result(res)
    integer(int32), intent(in) :: n
    type(linked_list), intent(in) :: edges
    integer(int32) :: i, counter(0:n)
    type(node), pointer :: now
    allocate(res%start(0:n), res%elist(edges%size))
    res%start = 0
    now => edges%head
    do while (associated(now))
      res%start(now%from) = res%start(now%from) + 1
      now => now%next
    end do
    do i = 1, n
      res%start(i) = res%start(i) + res%start(i - 1)
    end do
    counter = res%start
    now => edges%head
    do while (associated(now))
      i = now%from - 1
      counter(i) = counter(i) + 1
      res%elist(counter(i)) = now%to
      now => now%next
    end do
  end
  type(scc_graph) function new_scc_graph(n) result(res)
    integer(int32), intent(in) :: n
    res%n = n
  end
  integer(int32) function num_vertices(this) result(res)
    class(scc_graph), intent(inout) :: this
    res = this%n
  end
  subroutine add_edge(this, from, to)
    class(scc_graph), intent(inout) :: this
    integer(int32), intent(in) :: from, to
    call add_linked(this%edges, from, to)
  end
  type(array_list) function scc_ids(this) result(res)
    class(scc_graph), intent(in) :: this
    integer :: now_ord, group_num, v, to, i
    integer, dimension(this%n) :: low, ord, ids
    type(array_list) :: visited
    type(csr) :: g
    now_ord = 0
    group_num = 0
    low = 0
    ord = -1
    ids = 0
    g = new_csr(this%n, this%edges)
    allocate(visited%array(this%n))
    do i = 1, this%n
      if  (ord(i) /= -1) cycle
      visited%size = 0
      call scc_dfs(i, this%n, now_ord, group_num, g, visited, low, ord, ids)
    end do
    allocate(res%array(0:this%n))
    res%array(0) = group_num
    do i = 1, this%n
      res%array(i) = group_num - ids(i)
    end do
  end
  recursive subroutine scc_dfs(v, n, now_ord, group_num, g, visited, low, ord, ids)
    integer(int32), intent(in) :: v, n
    integer(int32), intent(inout) :: now_ord, group_num, low(:), ord(:), ids(:)
    type(csr), intent(in) :: g
    type(array_list), intent(inout) :: visited
    integer(int32) :: i, to, u
    low(v) = now_ord
    ord(v) = now_ord
    now_ord = now_ord + 1
    call add_array(visited, v)
    do i = g%start(v - 1) + 1, g%start(v)
      to = g%elist(i)
      if (ord(to) == -1) then
        call scc_dfs(to, n, now_ord, group_num, g, visited, low, ord, ids)
        low(v) = min(low(v), low(to))
      else
        low(v) = min(low(v), ord(to))
      end if
    end do
    if (low(v) == ord(v)) then
      do
        u = pop_array(visited)
        ord(u) = n
        ids(u) = group_num
        if (u == v) exit
      end do
      group_num = group_num + 1
    end if
  end
  type(vertex_group) function scc(this) result(res)
    class(scc_graph), intent(in) :: this
    type(array_list) :: ids
    integer(int32) :: group_num, i, x
    integer(int32), allocatable :: counts(:)
    ids = scc_ids(this)
    group_num = ids%array(0)
    allocate(counts(group_num))
    counts = 0
    do i = 1, this%n
      x = ids%array(i)
      counts(x) = counts(x) + 1
    end do
    allocate(res%group(group_num))
    do i = 1, group_num
      allocate(res%group(i)%array(counts(i)))
    end do
    do i = 1, this%n
      call add_array(res%group(ids%array(i)), i)
    end do
  end
end module
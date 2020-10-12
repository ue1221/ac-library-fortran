module module_maxflow
  use, intrinsic :: iso_fortran_env
  implicit none
  integer(int32), parameter :: intkind = int32
  integer(intkind), private, parameter :: infty = lshift(1_intkind, 8 * intkind - 2)
  type, private :: edge_
    integer(int32) :: to, rev
    integer(intkind) :: cap
  end type
  type, private :: edge_list_
    integer(int32) :: size = 0
    type(edge_), allocatable :: array(:)
  end type
  type, private :: pair_
    integer(int32) :: first, second
  end type
  type, private :: pair_list_
    integer(int32) :: size = 0
    type(pair_), allocatable :: array(:)
  end type
  type, private :: queue_
    integer(int32) :: head = 0
    integer(int32) :: tail = 0
    integer(int32), allocatable :: array(:)
  end type
  type edge
    integer(int32) :: from, to
    integer(intkind) :: cap, flow
  end type
  type mf_graph
    integer(int32), private :: n
    type(pair_list_), private :: pos
    type(edge_list_), private, allocatable :: g(:)
  contains
    procedure :: add_edge => add_edge
    procedure :: get_edge => get_edge
    procedure :: edges => edges
    procedure :: change_edge => change_edge
    procedure, private :: flow0 => flow0
    procedure, private :: flow1 => flow1
    generic :: flow => flow0, flow1
  end type
contains
  type(edge_) function new_edge_(to, rev, cap) result(res)
    integer(int32), intent(in) :: to, rev
    integer(intkind), intent(in) :: cap
    res%to = to
    res%rev = rev
    res%cap = cap
  end
  subroutine add_edge_(this, to, rev, cap)
    type(edge_list_), intent(inout) :: this
    integer(int32), intent(in) :: to, rev
    integer(intkind), intent(in) :: cap
    if (.not.allocated(this%array)) allocate(this%array(1))
    if (size(this%array) == this%size) call append_edge_(this, this%size)
    this%size = this%size + 1
    this%array(this%size) = new_edge_(to, rev, cap)
  end
  subroutine append_edge_(this, size)
    type(edge_list_), intent(inout) :: this
    integer(int32), intent(in) :: size
    type(edge_) :: array(size)
    array = this%array
    deallocate(this%array)
    allocate(this%array(2 * size))
    this%array(1:size) = array
  end
  type(pair_) function new_pair_(first, second) result(res)
    integer(int32), intent(in) :: first, second
    res%first = first
    res%second = second
  end
  subroutine add_pair_(this, first, second)
    type(pair_list_), intent(inout) :: this
    integer(int32), intent(in) :: first, second
    if (.not.allocated(this%array)) allocate(this%array(1))
    if (size(this%array) == this%size) call append_pair_(this, this%size)
    this%size = this%size + 1
    this%array(this%size) = new_pair_(first, second)
  end
  subroutine append_pair_(this, size)
    type(pair_list_), intent(inout) :: this
    integer(int32), intent(in) :: size
    type(pair_) :: array(size)
    array = this%array
    deallocate(this%array)
    allocate(this%array(2 * size))
    this%array(1:size) = array
  end
  subroutine add_queue_(this, val)
    type(queue_), intent(inout) :: this
    integer(int32), intent(in) :: val
    this%tail = this%tail + 1
    this%array(this%tail) = val
  end
  integer(int32) function pop_queue_(this) result(res)
    type(queue_), intent(inout) :: this
    this%head = this%head + 1
    res = this%array(this%head)
  end
  subroutine reset_queue_(this)
    type(queue_), intent(inout) :: this
    this%head = 0
    this%tail = 0
  end
  type(mf_graph) function new_mf_graph(n) result(res)
    integer(int32), intent(in) :: n
    res%n = n
    allocate(res%g(n))
  end
  subroutine add_edge(this, from, to, cap)
    class(mf_graph), intent(inout) :: this
    integer(int32), intent(in) :: from, to
    integer(intkind), intent(in) :: cap
    integer(int32) :: from_id, to_id
    from_id = this%g(from)%size + 1
    to_id = this%g(to)%size + 1
    call add_pair_(this%pos, from, from_id)
    if (from == to) to_id = to_id + 1
    call add_edge_(this%g(from), to, to_id, cap)
    call add_edge_(this%g(to), from, from_id, 0_intkind)
  end
  type(edge) function get_edge(this, i) result(res)
    class(mf_graph), intent(in) :: this
    integer(int32), intent(in) :: i
    type(edge_) :: e, re
    e = this%g(this%pos%array(i)%first)%array(this%pos%array(i)%second)
    re = this%g(e%to)%array(e%rev)
    res%from = this%pos%array(i)%first
    res%to = e%to
    res%cap = e%cap + re%cap
    res%flow = re%cap
  end
  function edges(this) result(res)
    class(mf_graph), intent(in) :: this
    type(edge) :: res(this%pos%size)
    integer(int32) :: i
    do i = 1, this%pos%size
      res(i) = get_edge(this, i)
    end do
  end
  subroutine change_edge(this, i, new_cap, new_flow)
    class(mf_graph), intent(inout) :: this
    integer(int32), intent(in) :: i
    integer(intkind), intent(in) :: new_cap, new_flow
    type(edge_) :: e
    e = this%g(this%pos%array(i)%first)%array(this%pos%array(i)%second)
    this%g(this%pos%array(i)%first)%array(this%pos%array(i)%second)%cap = &
      new_cap - new_flow
    this%g(e%to)%array(e%rev)%cap = new_flow
  end
  integer(intkind) function flow0(this, s, t) result(res)
    class(mf_graph), intent(inout) :: this
    integer(int32), intent(in) :: s, t
    res = flow1(this, s, t, infty)
  end
  integer(intkind) function flow1(this, s, t, flow_limit) result(res)
    class(mf_graph), intent(inout) :: this
    integer(int32), intent(in) :: s, t
    integer(intkind), intent(in) :: flow_limit
    integer(int32), dimension(this%n) :: level, iter
    integer(intkind) :: flow
    type(queue_) :: queue
    allocate(queue%array(this%n))
    res = 0
    do while (res < flow_limit)
      call flow_bfs(this%g, s, t, level, queue)
      if (level(t) == -1) exit
      iter = 1
      do while (res < flow_limit)
        flow = flow_dfs(this%g, t, s, level, iter, flow_limit - res)
        if (flow == 0) exit
        res = res + flow
      end do
    end do
  end
  subroutine flow_bfs(g, s, t, level, queue)
    type(edge_list_), intent(inout) :: g(:)
    integer(int32), intent(in) :: s, t
    integer(int32), intent(inout) :: level(:)
    type(queue_), intent(inout) :: queue
    integer(int32) :: v, i
    type(edge_) :: e
    level = -1
    level(s) = 0
    call reset_queue_(queue)
    call add_queue_(queue, s)
    do while (queue%tail - queue%head > 0)
      v = pop_queue_(queue)
      do i = 1, g(v)%size
        e = g(v)%array(i)
        if (e%cap == 0 .or. level(e%to) >= 0) cycle
        level(e%to) = level(v) + 1
        if (e%to == t) return
        call add_queue_(queue, e%to)
      end do
    end do
  end
  recursive integer(intkind) function flow_dfs(g, v, s, level, iter, up) result(res)
    type(edge_list_), intent(inout) :: g(:)
    integer(int32), intent(in) :: v, s, level(:)
    integer(int32), intent(inout) :: iter(:)
    integer(intkind), intent(in) :: up
    integer(int32) :: level_v, i
    integer(intkind) :: d
    type(edge_) :: e
    res = up
    if (v == s) return
    res = 0
    level_v = level(v)
    do i = iter(v), g(v)%size
      e = g(v)%array(i)
      if (level_v <= level(e%to) .or. g(e%to)%array(e%rev)%cap == 0) cycle
      d = flow_dfs(g, e%to, s, level, iter, min(up - res, g(e%to)%array(e%rev)%cap))
      if (d <= 0) cycle
      g(v)%array(i)%cap = g(v)%array(i)%cap + d
      g(e%to)%array(e%rev)%cap = g(e%to)%array(e%rev)%cap - d
      res = res + d
      if (res == up) exit
      iter(v) = iter(v) + 1
    end do
  end
  function min_cut(this, s) result(res)
    class(mf_graph), intent(in) :: this
    integer(int32), intent(in) :: s
    logical :: res(this%n)
    integer(int32) :: p, i
    type(edge_) :: e
    type(queue_) :: queue
    res = .false.
    allocate(queue%array(this%n))
    call add_queue_(queue, s)
    do while (queue%tail - queue%head > 0)
      p = pop_queue_(queue)
      res(p) = .true.
      do i = 1, this%g(p)%size
        e = this%g(p)%array(i)
        if (e%cap > 0 .and. .not.res(e%to)) then
          res(e%to) = .true.
          call add_queue_(queue, e%to)
        end if
      end do
    end do
  end
end module
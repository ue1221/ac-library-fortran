module module_mincostflow
  use, intrinsic :: iso_fortran_env
  implicit none
  integer(int32), parameter :: capkind = int32
  integer(int32), parameter :: costkind = int64
  integer(capkind), private, parameter :: infty_cap = lshift(1_capkind, 8 * capkind - 2)
  integer(costkind), private, parameter :: infty_cost = lshift(1_costkind, 8 * costkind - 2)
  type, private :: edge_
    integer(int32) :: to, rev
    integer(capkind) :: cap
    integer(costkind) :: cost
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
  type, private :: node
    integer(costkind) :: key
    integer(int32) :: to
    integer(int32) :: rank = 1
    type(node), pointer :: left => null(), right => null()
  end type
  type, private :: leftist_heap
    integer :: size = 0
    type(node), pointer :: root => null()
  end type
  type edge
    integer(int32) :: from, to
    integer(capkind) :: cap, flow
    integer(costkind) :: cost
  end type
  type pair
    integer(capkind) :: flow
    integer(costkind) :: cost
  end type
  type pair_list
    integer(int32) :: size = 0
    type(pair), allocatable :: array(:)
  end type
  type mcf_graph
    integer(int32), private :: n
    type(pair_list_), private :: pos
    type(edge_list_), private, allocatable :: g(:)
  contains
    procedure :: add_edge => add_edge
    procedure :: get_edge => get_edge
    procedure :: edges => edges
    procedure, private :: flow0 => flow0
    procedure, private :: flow1 => flow1
    generic :: flow => flow0, flow1
    procedure, private :: slope0 => slope0
    procedure, private :: slope1 => slope1
    generic :: slope => slope0, slope1
  end type
contains
  type(edge_) function new_edge_(to, rev, cap, cost) result(res)
    integer(int32), intent(in) :: to, rev
    integer(capkind), intent(in) :: cap
    integer(costkind), intent(in) :: cost
    res%to = to
    res%rev = rev
    res%cap = cap
    res%cost = cost
  end
  subroutine add_edge_(this, to, rev, cap, cost)
    type(edge_list_), intent(inout) :: this
    integer(int32), intent(in) :: to, rev
    integer(capkind), intent(in) :: cap
    integer(costkind), intent(in) :: cost
    if (.not.allocated(this%array)) allocate(this%array(1))
    if (size(this%array) == this%size) call append_edge_(this, this%size)
    this%size = this%size + 1
    this%array(this%size) = new_edge_(to, rev, cap, cost)
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
  function new_node(key, to) result(res)
    integer(costkind), intent(in) :: key
    integer(int32), intent(in) :: to
    type(node), pointer :: res
    allocate(res)
    res%key = key
    res%to = to
  end
  recursive function meld_node(a, b) result(res)
    type(node), pointer, intent(in) :: a, b
    type(node), pointer :: res, tmp
    logical :: f
    if (.not.associated(a)) then
      res => b
      return
    end if
    if (.not.associated(b)) then
      res => a
      return
    end if
    if (a%key < b%key) then
      res => a
      res%right => meld_node(res%right, b)
    else
      res => b
      res%right => meld_node(res%right, a)
    end if
    f = .not.associated(res%left)
    if (.not.f) f = res%left%rank < res%right%rank
    if (f) then
      tmp => res%left
      res%left => res%right
      res%right => tmp
    end if
    res%rank = merge(res%right%rank, 0, associated(res%right)) + 1
  end
  subroutine add_node(this, key, to)
    type(leftist_heap), intent(inout) :: this
    integer(costkind), intent(in) :: key
    integer(int32), intent(in) :: to
    this%root => meld_node(this%root, new_node(key, to))
    this%size = this%size + 1
  end
  integer(int32) function pop_node(this) result(res)
    type(leftist_heap), intent(inout) :: this
    type(node), pointer :: left, right
    res = this%root%to
    left => this%root%left
    right => this%root%right
    deallocate(this%root)
    this%root => meld_node(left, right)
    this%size = this%size - 1
  end
  type(pair) function new_pair(flow, cost) result(res)
    integer(capkind), intent(in) :: flow
    integer(costkind), intent(in) :: cost
    res%flow = flow
    res%cost = cost
  end
  subroutine add_pair(this, flow, cost)
    type(pair_list), intent(inout) :: this
    integer(capkind), intent(in) :: flow
    integer(costkind), intent(in) :: cost
    if (.not.allocated(this%array)) allocate(this%array(1))
    if (size(this%array) == this%size) call append_pair(this, this%size)
    this%size = this%size + 1
    this%array(this%size) = new_pair(flow, cost)
  end
  subroutine append_pair(this, size)
    type(pair_list), intent(inout) :: this
    integer(int32), intent(in) :: size
    type(pair) :: array(size)
    array = this%array
    deallocate(this%array)
    allocate(this%array(2 * size))
    this%array(1:size) = array
  end
  type(mcf_graph) function new_mcf_graph(n) result(res)
    integer(int32), intent(in) :: n
    res%n = n
    allocate(res%g(n))
  end
  subroutine add_edge(this, from, to, cap, cost)
    class(mcf_graph), intent(inout) :: this
    integer(int32), intent(in) :: from, to
    integer(capkind), intent(in) :: cap
    integer(costkind), intent(in) :: cost
    integer(int32) :: from_id, to_id
    from_id = this%g(from)%size + 1
    to_id = this%g(to)%size + 1
    call add_pair_(this%pos, from, from_id)
    if (from == to) to_id = to_id + 1
    call add_edge_(this%g(from), to, to_id, cap, cost)
    call add_edge_(this%g(to), from, from_id, 0_capkind, -cost)
  end
  type(edge) function get_edge(this, i) result(res)
    class(mcf_graph), intent(in) :: this
    integer(int32), intent(in) :: i
    type(edge_) :: e, re
    e = this%g(this%pos%array(i)%first)%array(this%pos%array(i)%second)
    re = this%g(e%to)%array(e%rev)
    res%from = this%pos%array(i)%first
    res%to = e%to
    res%cap = e%cap + re%cap
    res%flow = re%cap
    res%cost = e%cost
  end
  function edges(this) result(res)
    class(mcf_graph), intent(in) :: this
    type(edge) :: res(this%pos%size)
    integer(int32) :: i
    do i = 1, this%pos%size
      res(i) = get_edge(this, i)
    end do
  end
  type(pair) function flow0(this, s, t) result(res)
    class(mcf_graph), intent(inout) :: this
    integer(int32), intent(in) :: s, t
    res = flow1(this, s, t, infty_cap)
  end
  type(pair) function flow1(this, s, t, flow_limit) result(res)
    class(mcf_graph), intent(inout) :: this
    integer(int32), intent(in) :: s, t
    integer(capkind), intent(in) :: flow_limit
    type(pair_list) :: tmp
    tmp = slope1(this, s, t, flow_limit)
    res = tmp%array(tmp%size)
  end
  function slope0(this, s, t) result(res)
    class(mcf_graph), intent(inout) :: this
    integer(int32), intent(in) :: s, t
    type(pair_list) :: res
    res = slope1(this, s, t, infty_cap)
  end
  function slope1(this, s, t, flow_limit) result(res)
    class(mcf_graph), intent(inout) :: this
    integer(int32), intent(in) :: s, t
    integer(capkind), intent(in) :: flow_limit
    type(pair_list) :: res
    integer(costkind), dimension(this%n) :: dual, dist
    integer(int32), dimension(this%n) :: pv, pe
    logical, dimension(this%n) :: vis
    integer(capkind) :: flow, c
    integer(costkind) :: cost, prev_cost_per_flow, d
    integer(int32) :: v
    type(edge_) :: e
    dual = 0
    flow = 0
    cost = 0
    prev_cost_per_flow = -1
    call add_pair(res, flow, cost)
    do while (flow < flow_limit)
      if (.not.dual_ref(this%g, s, t, dual, dist, pv, pe, vis)) exit
      c = flow_limit - flow
      v = t
      do while (v /= s)
        c = min(c, this%g(pv(v))%array(pe(v))%cap)
        v = pv(v)
      end do
      v = t
      do while (v /= s)
        this%g(pv(v))%array(pe(v))%cap = this%g(pv(v))%array(pe(v))%cap - c
        e = this%g(pv(v))%array(pe(v))
        this%g(v)%array(e%rev)%cap = this%g(v)%array(e%rev)%cap + c
        v = pv(v)
      end do
      d = -dual(s)
      flow = flow + c
      cost = cost + int(c, costkind) * d
      if (prev_cost_per_flow == d) then
        if (res%size > 0) res%size = res%size - 1
      end if
      call add_pair(res, flow, cost)
      prev_cost_per_flow = d
    end do
  end
  logical function dual_ref(g, s, t, dual, dist, pv, pe, vis) result(res)
    type(edge_list_), intent(inout) :: g(:)
    integer(int32), intent(in) :: s, t
    integer(costkind), intent(inout) :: dual(:), dist(:)
    integer(int32), intent(out) :: pv(:), pe(:)
    logical, intent(inout) :: vis(:)
    type(leftist_heap) :: pq
    integer(int32) :: v, i
    integer(costkind) :: cost
    type(edge_) :: e
    res = .false.
    dist = infty_cost
    pv = -1
    pe = -1
    vis = .false.
    dist(s) = 0
    call add_node(pq, 0_costkind, s)
    do while (pq%size > 0)
      v = pop_node(pq)
      if (vis(v)) cycle
      vis(v) = .true.
      if (v == t) then
        do while (pq%size > 0)
          v = pop_node(pq)
        end do
        exit
      end if
      do i = 1, g(v)%size
        e = g(v)%array(i)
        if (vis(e%to) .or. e%cap == 0) cycle
        cost = e%cost - dual(e%to) + dual(v)
        if (dist(e%to) > dist(v) + cost) then
          dist(e%to) = dist(v) + cost
          pv(e%to) = v
          pe(e%to) = i
          call add_node(pq, dist(e%to), e%to)
        end if
      end do
    end do
    if (.not.vis(t)) return
    do v = 1, size(g)
      if (.not.vis(v)) cycle
      dual(v) = dual(v) - (dist(t) - dist(v))
    end do
    res = .true.
  end
end module
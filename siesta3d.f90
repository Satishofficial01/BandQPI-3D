!=======================================================================
!  siesta3d.f90   (v4)  --  3D band structure of 2D materials from SIESTA
!
!  usage:  siesta3d  system.fdf  [--fast | --fine | --budget=SECONDS] [--noshow]
!
!  Everything is taken from the fdf (SystemLabel, LatticeConstant,
!  BandLinesScale, %block BandLines, %include) and from <SystemLabel>.HSX,
!  the self-consistent Hamiltonian H(R) and overlap S(R) that SIESTA writes
!  (format v2 of recent SIESTA, or the original format of older SIESTA / HONPAS; "SaveHS true" if the file is missing).
!
!  PHYSICS
!   Bloch sums  H(k) = sum_R H(R) e^{ik.R},  S(k) = sum_R S(R) e^{ik.R}
!   generalised eigenproblem  H(k) c = E S(k) c   (LAPACK ZHEGV / ZHEGVX)
!   -> exact DFT eigenvalues at any k, the same operation SIESTA performs
!      for BandLines (Soler et al. 2002). No fitting, no interpolation.
!
!   Mirror symmetry (flat sheets, plane = xy): orbitals with (-1)^(l+|m|)=-1
!   (p_z, d_xz, d_yz, ...) are odd under z -> -z and do not mix with the
!   even ones. Diagonalising the odd block alone gives the pi / pi* bands
!   exactly (graphene, h-BN, C/BN hybrids, ...), i.e. the bands of the
!   tight-binding picture (Castro Neto et al. 2009).
!
!  GRIDS  (sizes chosen automatically from a time budget)
!   * periodic N x N grid of the BZ  -> DOS (linear triangles), electron
!     counts, band-edge search, JDOS, constant-energy contours
!   * whole BZ on a grid shaped like the Wigner-Seitz polygon (exact edges,
!     corners such as K are grid points)
!   * extended rectangular window (as VASPKIT task 231 / the MATLAB plot)
!     -> KX.grd KY.grd BAND_HOMO.grd BAND_LUMO.grd (VASPKIT format)
!   * zoom around the smallest direct gap
!
!  VALIDATION (written to validation.txt, printed at the end)
!   bands vs SIESTA <label>.bands, Hermiticity of H(k), positivity of S(k),
!   periodicity E(k+G)=E(k), time reversal E(-k)=E(k), electron count,
!   band edges vs path; the plotting script adds surface-vs-path and
!   DOS-vs-SIESTA checks.
!
!  Build:  gfortran -O3 -march=native -fopenmp siesta3d.f90 -lopenblas -o siesta3d.x
!          (or -llapack -lblas)
!=======================================================================
module consts
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi       = 3.14159265358979323846_dp
  real(dp), parameter :: Ry2eV    = 13.605693122994_dp
  real(dp), parameter :: Bohr2A   = 0.529177210903_dp
  real(dp), parameter :: hbar_eVs = 6.582119569e-16_dp
  real(dp), parameter :: hb2m     = 7.619964231_dp       ! hbar^2/m_e (eV Ang^2)
end module consts

!=======================================================================
module fdf_lite
  use consts
  implicit none
  integer, parameter :: maxl = 20000
  character(len=512) :: lines(maxl)
  integer :: nlines = 0
contains
  recursive subroutine fdf_load(fname)
    character(len=*), intent(in) :: fname
    integer :: u, ios, ic
    character(len=512) :: l
    open(newunit=u, file=fname, status='old', iostat=ios)
    if (ios /= 0) then
       write(*,'(a)') ' ERROR: cannot open fdf file '//trim(fname); stop 1
    end if
    do
       read(u,'(a)',iostat=ios) l
       if (ios /= 0) exit
       ic = index(l, '#'); if (ic > 0) l = l(1:ic-1)
       ic = index(l, '!'); if (ic > 0) l = l(1:ic-1)
       l = adjustl(l)
       if (len_trim(l) == 0) cycle
       if (norm(word(l,1)) == '%include') then
          call fdf_load(join(dirname(fname), word(l,2))); cycle
       end if
       if (nlines < maxl) then
          nlines = nlines + 1; lines(nlines) = l
       end if
    end do
    close(u)
  end subroutine fdf_load

  function norm(s) result(r)          ! fdf labels: case and . _ - ignored
    character(len=*), intent(in) :: s
    character(len=len(s)) :: r
    integer :: i, j, c
    r = ' '; j = 0
    do i = 1, len_trim(s)
       if (s(i:i)=='.' .or. s(i:i)=='_' .or. s(i:i)=='-') cycle
       j = j + 1; c = iachar(s(i:i))
       if (c >= 65 .and. c <= 90) then; r(j:j) = achar(c+32); else; r(j:j) = s(i:i); end if
    end do
  end function norm

  function word(l, n) result(w)
    character(len=*), intent(in) :: l
    integer, intent(in) :: n
    character(len=256) :: w
    integer :: i, k, st, le
    w = ' '; k = 0; i = 1; le = len_trim(l)
    do while (i <= le)
       do while (i <= le .and. (l(i:i)==' ' .or. l(i:i)==char(9)))
          i = i + 1
       end do
       if (i > le) exit
       st = i
       do while (i <= le .and. l(i:i)/=' ' .and. l(i:i)/=char(9))
          i = i + 1
       end do
       k = k + 1
       if (k == n) then; w = l(st:i-1); return; end if
    end do
  end function word

  function fdf_string(label, def) result(v)
    character(len=*), intent(in) :: label, def
    character(len=256) :: v
    integer :: i
    v = def
    do i = 1, nlines
       if (norm(word(lines(i),1)) == norm(label)) then; v = word(lines(i),2); return; end if
    end do
  end function fdf_string

  real(dp) function fdf_length(label, def_bohr) result(v)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: def_bohr
    integer :: i, ios
    character(len=256) :: un, w2
    v = def_bohr
    do i = 1, nlines
       if (norm(word(lines(i),1)) == norm(label)) then
          w2 = word(lines(i),2); read(w2, *, iostat=ios) v
          un = norm(word(lines(i),3))
          if (un == 'ang' .or. un == 'angstrom') v = v/Bohr2A
          if (un == 'nm') v = 10.0_dp*v/Bohr2A
          return
       end if
    end do
  end function fdf_length

  subroutine fdf_block(label, i1, i2)
    character(len=*), intent(in) :: label
    integer, intent(out) :: i1, i2
    integer :: i
    i1 = 0; i2 = -1
    do i = 1, nlines
       if (norm(word(lines(i),1)) == '%block' .and. norm(word(lines(i),2)) == norm(label)) then
          i1 = i + 1
          do i2 = i1, nlines
             if (norm(word(lines(i2),1)) == '%endblock') exit
          end do
          i2 = i2 - 1; return
       end if
    end do
  end subroutine fdf_block

  function dirname(f) result(d)
    character(len=*), intent(in) :: f
    character(len=512) :: d
    integer :: i
    i = index(f, '/', back=.true.)
    if (i == 0) then; d = '.'; else; d = f(1:i-1); end if
  end function dirname

  function join(d, f) result(p)
    character(len=*), intent(in) :: d, f
    character(len=1024) :: p
    if (f(1:1) == '/') then; p = f; else; p = trim(d)//'/'//trim(f); end if
  end function join
end module fdf_lite

!=======================================================================
module hsx
  use consts
  implicit none
  integer :: na_u, no_u, nspin, nspecies, nsc(3), n_s, kcell(3,3)
  ! spinor (non-collinear nspin=4 / spin-orbit nspin=8) Hamiltonians: the
  ! basis is (orbital, spin) with 2*no_u functions; nsraw = nspin in the file
  logical :: soc = .false.
  logical :: hsx_legacy = .false.     ! original HSX format (SIESTA 4.0 / HONPAS) read
  integer :: nsraw = 1, nbas = 0
  character(len=20), allocatable :: spec_label(:)
  real(dp) :: ucell(3,3), Ef, qtot, temp, rcell(3,3)
  integer, allocatable :: isc_off(:,:), isa(:), numh(:), listhptr(:), listh(:)
  integer, allocatable :: col(:), scell(:), orb_l(:), orb_m(:), orb_atom(:), orb_n(:)
  real(dp), allocatable :: xa(:,:), H(:,:), S(:), Rvec(:,:)
contains
  subroutine read_hsx(fname)
    character(len=*), intent(in) :: fname
    integer :: iu, version, is_dp, is, io, nnz, ispin, ios, ia, iorb, i, im
    integer, allocatable :: nos(:), lasto(:)
    type sp_t
       integer, allocatable :: l(:), n(:)
    end type sp_t
    type(sp_t), allocatable :: sp(:)
    integer, allocatable :: nq(:), zq(:)
    real(dp), allocatable :: zval(:)
    character(len=20), allocatable :: label(:)
    real(4), allocatable :: b4(:)
    real(dp) :: kd(3)
    open(newunit=iu, file=fname, form='unformatted', status='old', iostat=ios)
    if (ios /= 0) then
       write(*,'(a)') ' ERROR: cannot open '//trim(fname)
       write(*,'(a)') '        add  "SaveHS true"  to the fdf and rerun SIESTA.'
       stop 1
    end if
    read(iu, iostat=ios) version
    if (ios /= 0 .or. version /= 2) then
       ! not the versioned format v2 of recent SIESTA: try the original format written by
       ! SIESTA 4.0 / 4.1 (iohs / m_hsx) and by HONPAS (hybrid functionals)
       close(iu)
       call read_hsx_legacy(fname)
       return
    end if
    read(iu) is_dp
    read(iu) na_u, no_u, nspin, nspecies, nsc
    nsraw = nspin
    if (nspin /= 1 .and. nspin /= 2 .and. nspin /= 4 .and. nspin /= 8) then
       write(*,'(a,i0)') ' ERROR: unknown spin setting in HSX, nspin = ', nspin; stop 1
    end if
    soc = nspin >= 4
    nbas = merge(2*no_u, no_u, soc)
    n_s = product(nsc)
    read(iu) ucell, Ef, qtot, temp
    allocate(isc_off(3,n_s), xa(3,na_u), isa(na_u), lasto(na_u))
    read(iu) isc_off, xa, isa, lasto
    allocate(label(nspecies), zval(nspecies), nos(nspecies), sp(nspecies))
    read(iu) (label(is), zval(is), nos(is), is=1,nspecies)
    allocate(spec_label(nspecies)); spec_label = label
    do is = 1, nspecies
       allocate(nq(nos(is)), sp(is)%l(nos(is)), sp(is)%n(nos(is)), zq(nos(is)))
       read(iu) (nq(io), sp(is)%l(io), zq(io), io=1,nos(is))
       sp(is)%n = nq
       deallocate(nq, zq)
    end do
    read(iu) kcell, kd
    allocate(numh(no_u), listhptr(no_u))
    read(iu) numh
    nnz = sum(numh)
    listhptr(1) = 0
    do io = 2, no_u
       listhptr(io) = listhptr(io-1) + numh(io-1)
    end do
    allocate(listh(nnz), H(nnz,nspin), S(nnz))
    do io = 1, no_u
       read(iu) listh(listhptr(io)+1:listhptr(io)+numh(io))
    end do
    if (is_dp == 1) then
       do ispin = 1, nspin
          do io = 1, no_u
             read(iu) H(listhptr(io)+1:listhptr(io)+numh(io), ispin)
          end do
       end do
       do io = 1, no_u
          read(iu) S(listhptr(io)+1:listhptr(io)+numh(io))
       end do
    else
       allocate(b4(maxval(numh)))
       do ispin = 1, nspin
          do io = 1, no_u
             read(iu) b4(1:numh(io))
             H(listhptr(io)+1:listhptr(io)+numh(io), ispin) = b4(1:numh(io))
          end do
       end do
       do io = 1, no_u
          read(iu) b4(1:numh(io))
          S(listhptr(io)+1:listhptr(io)+numh(io)) = b4(1:numh(io))
       end do
    end if
    close(iu)
    allocate(col(nnz), scell(nnz), Rvec(3,n_s))
    col = mod(listh-1, no_u) + 1
    scell = (listh-1)/no_u + 1
    do is = 1, n_s
       Rvec(:,is) = matmul(ucell, real(isc_off(:,is),dp))
    end do
    ! orbital -> (atom, l, m); SIESTA real harmonics, m = -l..l per shell
    allocate(orb_l(no_u), orb_m(no_u), orb_atom(no_u), orb_n(no_u))
    orb_n = 0
    iorb = 0
    do ia = 1, na_u
       is = isa(ia); i = 1
       do while (i <= nos(is))
          do im = -sp(is)%l(i), sp(is)%l(i)
             iorb = iorb + 1
             if (iorb <= no_u) then
                orb_atom(iorb) = ia; orb_l(iorb) = sp(is)%l(i); orb_m(iorb) = im; orb_n(iorb) = sp(is)%n(i)
             end if
          end do
          i = i + 2*sp(is)%l(i) + 1
       end do
    end do
    if (iorb /= no_u) then
       write(*,'(a)') ' WARNING: orbital (l,m) bookkeeping mismatch; mirror analysis disabled'
       orb_l = -1
    end if
    ! orbital -> atom directly from lasto (always reliable)
    do ia = 1, na_u
       if (ia == 1) then
          orb_atom(1:lasto(1)) = 1
       else
          orb_atom(lasto(ia-1)+1:lasto(ia)) = ia
       end if
    end do
    call recip(ucell, rcell)
  end subroutine read_hsx

  !--- original HSX format (SIESTA 4.0, SIESTA 4.1 m_hsx, HONPAS):
  !    no_u, no_s, nspin, maxnh / gamma / indxuo(no_s) / numh / listh (supercell columns) /
  !    H (Ry, single) per spin / S / qtot, temp / xij (Bohr, single) / species, orbitals, atoms.
  !    It has no cell, no atomic positions, no Fermi level and no supercell indices: the cell and
  !    the positions are read from <label>.XV (written by every SIESTA / HONPAS run, Bohr), the
  !    lattice vector of every matrix element is R = xij - (x_j - x_i), and E_F comes from
  !    <label>.EIG (or .bands).
  subroutine read_hsx_legacy(fname)
    character(len=*), intent(in) :: fname
    integer :: iu, no_s, maxnh, nnz, io, ispin, ios, ia, is, i, j, k, jo, nR, iv(3), uxv, nl, ib
    integer :: ia0, ja0, ish, l0, cnt
    logical :: gamma, ex
    integer, allocatable :: indxuo(:), nofis(:), iaorb(:), iphorb(:), Rint(:,:), lofio(:,:), cfg(:,:), zt(:,:)
    real(4), allocatable :: b4(:), x4(:,:)
    real(dp), allocatable :: xij(:,:), zv(:)
    character(len=20), allocatable :: lab(:)
    character(len=1024) :: base, f2
    real(dp) :: r(3), fr(3), ainv(3,3), dmax
    integer, parameter :: maxR = 4000
    open(newunit=iu, file=fname, form='unformatted', status='old')
    read(iu, iostat=ios) no_u, no_s, nspin, maxnh
    if (ios /= 0 .or. no_u <= 0 .or. no_s < no_u .or. nspin < 1) then
       write(*,'(a)') ' ERROR: '//trim(fname)//' is neither the HSX format v2 (recent SIESTA) nor the'
       write(*,'(a)') '        original format (older SIESTA / HONPAS). Rerun with  SaveHS true.'; stop 1
    end if
    if (nspin /= 1 .and. nspin /= 2 .and. nspin /= 4 .and. nspin /= 8) then
       write(*,'(a,i0)') ' ERROR: unknown spin setting in HSX, nspin = ', nspin; stop 1
    end if
    nsraw = nspin; soc = nspin >= 4; nbas = merge(2*no_u, no_u, soc)
    read(iu) gamma
    allocate(indxuo(no_s))
    if (.not. gamma) then
       read(iu) indxuo
    else
       indxuo = [(i, i = 1, no_s)]
    end if
    allocate(numh(no_u), listhptr(no_u))
    read(iu) numh
    nnz = sum(numh)
    listhptr(1) = 0
    do io = 2, no_u
       listhptr(io) = listhptr(io-1) + numh(io-1)
    end do
    allocate(listh(nnz), H(nnz,nspin), S(nnz), b4(maxval(numh)), x4(3,maxval(numh)), xij(3,nnz))
    do io = 1, no_u
       read(iu) listh(listhptr(io)+1:listhptr(io)+numh(io))
    end do
    do ispin = 1, nspin
       do io = 1, no_u
          read(iu) b4(1:numh(io)); H(listhptr(io)+1:listhptr(io)+numh(io), ispin) = b4(1:numh(io))
       end do
    end do
    do io = 1, no_u
       read(iu) b4(1:numh(io)); S(listhptr(io)+1:listhptr(io)+numh(io)) = b4(1:numh(io))
    end do
    read(iu) qtot, temp
    xij = 0.0_dp
    if (.not. gamma) then
       do io = 1, no_u
          read(iu) x4(:,1:numh(io)); xij(:,listhptr(io)+1:listhptr(io)+numh(io)) = x4(:,1:numh(io))
       end do
    end if
    read(iu) nspecies
    allocate(lab(nspecies), zv(nspecies), nofis(nspecies))
    read(iu) (lab(is), zv(is), nofis(is), is=1,nspecies)
    allocate(spec_label(nspecies)); spec_label = lab
    allocate(lofio(maxval(nofis),nspecies), cfg(maxval(nofis),nspecies), zt(maxval(nofis),nspecies))
    do is = 1, nspecies
       do io = 1, nofis(is)
          read(iu) cfg(io,is), lofio(io,is), zt(io,is)
       end do
    end do
    read(iu) na_u
    allocate(isa(na_u), iaorb(no_u), iphorb(no_u))
    read(iu) isa
    read(iu) (iaorb(io), iphorb(io), io=1,no_u)
    close(iu)
    ! ---- cell and positions from <label>.XV
    base = fname(1:len_trim(fname)-4)
    f2 = trim(base)//'.XV'
    inquire(file=f2, exist=ex)
    if (.not. ex) then
       write(*,'(a)') ' ERROR: this HSX file (older SIESTA / HONPAS format) has no cell; '//trim(f2)//' is needed.'
       stop 1
    end if
    open(newunit=uxv, file=f2, status='old')
    do i = 1, 3
       read(uxv,*) ucell(:,i)
    end do
    read(uxv,*) nl
    if (nl /= na_u) then
       write(*,'(a)') ' ERROR: '//trim(f2)//' has a different number of atoms than the HSX file'; stop 1
    end if
    allocate(xa(3,na_u))
    do ia = 1, na_u
       read(uxv,*) i, j, xa(:,ia)
    end do
    close(uxv)
    ! ---- E_F (eV in .EIG / .bands)
    Ef = 0.0_dp
    f2 = trim(base)//'.EIG'; inquire(file=f2, exist=ex)
    if (.not. ex) then; f2 = trim(base)//'.bands'; inquire(file=f2, exist=ex); end if
    if (ex) then
       open(newunit=uxv, file=f2, status='old'); read(uxv,*) Ef; close(uxv)
       Ef = Ef/Ry2eV
    else
       write(*,'(a)') ' WARNING: no .EIG or .bands next to the HSX: E_F set from the electron count'
    end if
    ! ---- lattice vector of every element: R = xij - (x_jo - x_io), integer in units of the cell
    call recip(ucell, rcell)
    ainv = transpose(rcell)/(2.0_dp*pi)
    allocate(Rint(3,maxR), col(nnz), scell(nnz))
    nR = 1; Rint(:,1) = 0; dmax = 0.0_dp
    do io = 1, no_u
       ia0 = iaorb(io)
       do i = listhptr(io)+1, listhptr(io)+numh(io)
          jo = indxuo(listh(i)); ja0 = iaorb(jo); col(i) = jo
          if (gamma) then; scell(i) = 1; cycle; end if
          r = xij(:,i) - (xa(:,ja0) - xa(:,ia0))
          fr = matmul(ainv, r); iv = nint(fr)
          dmax = max(dmax, maxval(abs(matmul(ucell, real(iv,dp)) - r)))
          k = 0
          do j = 1, nR
             if (all(Rint(:,j) == iv)) then; k = j; exit; end if
          end do
          if (k == 0) then
             nR = nR + 1
             if (nR > maxR) then; write(*,'(a)') ' ERROR: too many lattice vectors in HSX'; stop 1; end if
             Rint(:,nR) = iv; k = nR
          end if
          scell(i) = k
       end do
    end do
    if (dmax > 0.05_dp) then
       write(*,'(a,f8.4,a)') ' ERROR: HSX and '//trim(base)//'.XV do not match (lattice vectors off by ', dmax, ' Bohr)'
       stop 1
    end if
    n_s = nR
    allocate(isc_off(3,n_s), Rvec(3,n_s))
    isc_off = Rint(:,1:n_s)
    do i = 1, 3
       nsc(i) = 2*maxval(abs(isc_off(i,:))) + 1
    end do
    do i = 1, n_s
       Rvec(:,i) = matmul(ucell, real(isc_off(:,i),dp))
    end do
    kcell = 0
    ! listh now holds unit-cell columns + supercell index, as for v2 files
    do i = 1, nnz
       listh(i) = col(i) + (scell(i)-1)*no_u
    end do
    ! ---- orbitals: atom, l, m
    allocate(orb_l(no_u), orb_m(no_u), orb_atom(no_u), orb_n(no_u))
    do io = 1, no_u
       orb_atom(io) = iaorb(io)
       is = isa(iaorb(io)); ib = iphorb(io)
       orb_l(io) = lofio(ib,is); orb_n(io) = cfg(ib,is)
       ! m: position inside its (2l+1) shell; shells are consecutive runs of equal l, n and zeta
       ish = ib; cnt = 0; l0 = lofio(ib,is)
       do while (ish > 1)
          if (lofio(ish-1,is) /= l0 .or. cfg(ish-1,is) /= cfg(ib,is) .or. zt(ish-1,is) /= zt(ib,is)) exit
          ish = ish - 1
       end do
       orb_m(io) = mod(ib - ish, 2*l0 + 1) - l0
    end do
    hsx_legacy = .true.
  end subroutine read_hsx_legacy

  subroutine recip(a, b)
    real(dp), intent(in) :: a(3,3)
    real(dp), intent(out) :: b(3,3)
    b(:,1) = cross(a(:,2), a(:,3)); b(:,2) = cross(a(:,3), a(:,1)); b(:,3) = cross(a(:,1), a(:,2))
    b = 2.0_dp*pi*b/dot_product(a(:,1), b(:,1))
  end subroutine recip

  pure function cross(x, y) result(z)
    real(dp), intent(in) :: x(3), y(3)
    real(dp) :: z(3)
    z = [x(2)*y(3)-x(3)*y(2), x(3)*y(1)-x(1)*y(3), x(1)*y(2)-x(2)*y(1)]
  end function cross

  !--- H(k), S(k) (upper triangle) of the orbitals with map(io) > 0 --------
  subroutine build_hk(k, ss, map, nb, Hk, Sk)
    real(dp), intent(in) :: k(3)
    integer, intent(in) :: ss, map(no_u), nb
    complex(dp), intent(out) :: Hk(nb,nb), Sk(nb,nb)
    complex(dp) :: ph(n_s)
    integer :: io, ind, a, b, j
    do j = 1, n_s
       ph(j) = exp(cmplx(0.0_dp, dot_product(k, Rvec(:,j)), dp))
    end do
    Hk = (0.0_dp,0.0_dp); Sk = (0.0_dp,0.0_dp)
    if (soc) then
       call build_hk_spinor(ph, Hk, Sk); return
    end if
    do io = 1, no_u
       a = map(io); if (a == 0) cycle
       do ind = listhptr(io)+1, listhptr(io)+numh(io)
          b = map(col(ind)); if (b < a) cycle      ! also skips b == 0
          Hk(a,b) = Hk(a,b) + H(ind,ss)*ph(scell(ind))
          Sk(a,b) = Sk(a,b) + S(ind)*ph(scell(ind))
       end do
    end do
    do a = 1, nb
       Hk(a,a) = cmplx(real(Hk(a,a)), 0.0_dp, dp); Sk(a,a) = cmplx(real(Sk(a,a)), 0.0_dp, dp)
    end do
  end subroutine build_hk

  !--- 2x2 spin box of file element ind, as it enters H(k) at (jo, io):
  !    nspin=8: components [1..8]; uu = 1 + i5, dd = 2 + i6, ud = 3 - i4, du = 7 + i8
  !    nspin=4: uu = 1, dd = 2, ud = 3 - i4, du = 3 + i4
  !    (checked element by element against sisl and against SIESTA's .bands
  !     for a SIESTA 5.4.2 spin-orbit run: |dH| < 1e-12 Ry)
  pure function spin_box(ind) result(M)
    integer, intent(in) :: ind
    complex(dp) :: M(2,2)
    if (nsraw == 8) then
       M(1,1) = cmplx(H(ind,1), H(ind,5), dp); M(1,2) = cmplx(H(ind,3), -H(ind,4), dp)
       M(2,1) = cmplx(H(ind,7), H(ind,8), dp); M(2,2) = cmplx(H(ind,2), H(ind,6), dp)
    else
       M(1,1) = H(ind,1); M(2,2) = H(ind,2)
       M(1,2) = cmplx(H(ind,3), -H(ind,4), dp); M(2,1) = cmplx(H(ind,3), H(ind,4), dp)
    end if
  end function spin_box

  !--- spinor H(k), S(k) (full matrices, 2*no_u). SIESTA stores the sparse
  !    matrix transposed in the orbital index: element (io, jo + R) of the
  !    file is <jo|H|io + R>^T, i.e. it enters H(k) at (jo, io) with phase
  !    exp(-i k.R) (spin box not transposed), as in sisl.
  subroutine build_hk_spinor(ph, Hk, Sk)
    complex(dp), intent(in) :: ph(n_s)
    complex(dp), intent(inout) :: Hk(2*no_u,2*no_u), Sk(2*no_u,2*no_u)
    complex(dp) :: M(2,2), p
    integer :: io, jo, ind
    do io = 1, no_u
       do ind = listhptr(io)+1, listhptr(io)+numh(io)
          jo = col(ind); p = conjg(ph(scell(ind)))
          M = spin_box(ind)*p
          Hk(2*jo-1:2*jo, 2*io-1:2*io) = Hk(2*jo-1:2*jo, 2*io-1:2*io) + M
          Sk(2*jo-1, 2*io-1) = Sk(2*jo-1, 2*io-1) + S(ind)*p
          Sk(2*jo,   2*io)   = Sk(2*jo,   2*io)   + S(ind)*p
       end do
    end do
  end subroutine build_hk_spinor

  !--- eigenvalues il..iu (eV, absolute) of the block selected by map ------
  subroutine eigk(k, ss, map, nb, il, iu, e)
    real(dp), intent(in) :: k(3)
    integer, intent(in) :: ss, map(no_u), nb, il, iu
    real(dp), intent(out) :: e(iu-il+1)
    complex(dp), allocatable :: Hk(:,:), Sk(:,:), work(:), z(:,:)
    real(dp), allocatable :: rwork(:), w(:)
    integer, allocatable :: iwork(:), ifail(:)
    integer :: info, m
    allocate(Hk(nb,nb), Sk(nb,nb), w(nb), rwork(7*nb), iwork(5*nb), ifail(nb), work(3*nb))
    call build_hk(k, ss, map, nb, Hk, Sk)
    if (il == 1 .and. iu == nb) then
       call zhegv(1, 'N', 'U', nb, Hk, nb, Sk, nb, w, work, 3*nb, rwork, info)
    else
       allocate(z(1,1))
       call zhegvx(1, 'N', 'I', 'U', nb, Hk, nb, Sk, nb, 0.0_dp, 0.0_dp, il, iu, &
                   0.0_dp, m, w, z, 1, work, 3*nb, rwork, iwork, ifail, info)
    end if
    if (info /= 0) then
       write(*,*) 'LAPACK info =', info, '  k (1/Bohr) =', k
       stop ' diagonalisation failed (overlap S(k) not positive definite?)'
    end if
    e = w(1:iu-il+1)*Ry2eV
  end subroutine eigk

  !--- validation helpers: full H(k) both triangles, no symmetrisation -----
  subroutine herm_check(k, ss, dev, smin)
    real(dp), intent(in) :: k(3)
    integer, intent(in) :: ss
    real(dp), intent(out) :: dev, smin
    complex(dp), allocatable :: Hk(:,:), Sk(:,:), work(:)
    real(dp), allocatable :: w(:), rwork(:)
    complex(dp) :: ph(n_s)
    integer :: io, ind, j, info
    allocate(Hk(nbas,nbas), Sk(nbas,nbas), w(nbas), rwork(3*nbas), work(3*nbas))
    do j = 1, n_s
       ph(j) = exp(cmplx(0.0_dp, dot_product(k, Rvec(:,j)), dp))
    end do
    Hk = 0; Sk = 0
    if (soc) then
       call build_hk_spinor(ph, Hk, Sk)
       dev = maxval(abs(Hk - conjg(transpose(Hk))))*Ry2eV
       call zheev('N', 'U', nbas, Sk, nbas, w, work, 3*nbas, rwork, info)
       smin = w(1); return
    end if
    do io = 1, no_u
       do ind = listhptr(io)+1, listhptr(io)+numh(io)
          Hk(io,col(ind)) = Hk(io,col(ind)) + H(ind,ss)*ph(scell(ind))
          Sk(io,col(ind)) = Sk(io,col(ind)) + S(ind)*ph(scell(ind))
       end do
    end do
    dev = maxval(abs(Hk - conjg(transpose(Hk))))*Ry2eV
    call zheev('N', 'U', no_u, Sk, no_u, w, work, 3*no_u, rwork, info)
    smin = w(1)
  end subroutine herm_check
end module hsx

!=======================================================================
module fftmod
  use consts
  implicit none
contains
  subroutine fft1d(x, isign)
    complex(dp), intent(inout) :: x(0:)
    integer, intent(in) :: isign
    integer :: n, i, j, m, mmax, istep
    complex(dp) :: w, wp, t
    n = size(x); j = 0
    do i = 0, n-2
       if (i < j) then; t = x(j); x(j) = x(i); x(i) = t; end if
       m = n/2
       do while (m >= 1 .and. j >= m)
          j = j - m; m = m/2
       end do
       j = j + m
    end do
    mmax = 1
    do while (mmax < n)
       istep = 2*mmax
       wp = exp(cmplx(0.0_dp, isign*pi/mmax, dp)); w = (1.0_dp,0.0_dp)
       do m = 0, mmax-1
          do i = m, n-1, istep
             j = i + mmax; t = w*x(j); x(j) = x(i) - t; x(i) = x(i) + t
          end do
          w = w*wp
       end do
       mmax = istep
    end do
  end subroutine fft1d
  subroutine fft2d(a, isign)
    complex(dp), intent(inout) :: a(0:,0:)
    integer, intent(in) :: isign
    integer :: i
    complex(dp), allocatable :: c(:)
    allocate(c(0:size(a,1)-1))
    do i = 0, size(a,2)-1
       call fft1d(a(:,i), isign)
    end do
    do i = 0, size(a,1)-1
       c = a(i,:); call fft1d(c, isign); a(i,:) = c
    end do
  end subroutine fft2d
end module fftmod

!=======================================================================
program siesta3d
  ! v6: 1D / 2D / 3D periodic systems, exact linear-simplex integration,
  !     own Fermi level, exact constant-energy contours, contour JDOS and
  !     Bloch-overlap (wave-function) weighted QPI in the extended zone,
  !     QPI energy sweep for E(q) line cuts, grid-convergence checks
  use consts
  use fdf_lite
  use hsx
  !$ use omp_lib
  implicit none
  character(len=1024) :: fdfname, hsxname, outdir, arg, exedir, cmd, fdfdir
  character(len=2048) :: plotargs = ''
  character(len=256)  :: slabel, bscale
  logical :: show, ex, user_qpi
  real(dp) :: budget, alat, t0, tdiag, twall0, qpi_user(3)
  integer :: nthreads, iarg, i, u, uv, uv2, nqpi_user, npreq = 0, ncores = 1, est
  integer :: nqpi_tmp = 0
  real(dp) :: ewtmp(16)
  character(len=4096) :: relaunch
  character(len=64) :: envv
  logical :: colour = .false., npuser = .false.
  ! --- geometry --------------------------------------------------------
  integer :: per(3), dimk, ip(3)
  real(dp) :: Bp(3,3), B1(3), B2(3), eu(3), ev(3), en(3), b1p(2), b2p(2)
  ! --- orbital sectors -------------------------------------------------
  integer, allocatable :: map_all(:), map_odd(:)
  integer :: n_odd
  logical :: planar
  ! --- bands -----------------------------------------------------------
  integer :: blo, bhi, nbw, olo, ohi, nbo, nocc(2), nocc_o(2), ib1, ib2, ob1, ob2
  real(dp) :: Nel(2), Nel_o(2), Ef_eV, EfOwn, ewin, dosEF(2), kT
  real(dp), allocatable :: Eper(:,:,:,:), Eoper(:,:,:,:), E3(:,:,:,:,:), E1(:,:,:)
  real(dp), allocatable :: egrid(:), dos(:,:), dos_o(:,:)
  integer :: Np, Na, Nb, N3, N1, m_bz, m_ext, nkz
  ! --- band edges (k stored as plane-frame coordinates u,v,n in 1/Ang) -
  real(dp) :: kvbm(3,2), kcbm(3,2), evbm(2), ecbm(2), kdg(3,2), edg(2), gap
  real(dp) :: mvb(3,2), mcb(3,2), vF(2), kzc(3), kdgo(3,2), edgo(2), vFo(2)
  integer :: s_v, s_c, s_d
  logical :: metal(2), dirac(2), diraco(2), spin_degenerate
  ! --- BZ polygon (plane) -----------------------------------------------
  real(dp) :: V(2,24), bz_ext(2)
  integer :: nv
  ! --- path ------------------------------------------------------------
  integer :: npath_pts, nnodes
  real(dp), allocatable :: kpath(:,:), xpath(:), Epath(:,:,:)
  real(dp) :: knode(3,60), xnode(60)
  character(len=16) :: lnode(60)
  integer :: prog_last = -10
  ! --- QPI -------------------------------------------------------------
  integer :: nqe
  real(dp) :: qpiE(3)
  ! --- v6: spin-degenerate shortcut, eigenvectors for Bloch-weighted QPI
  integer :: nsp_eff = 1, qlo = 1, qhi = 1, nqb = 0, nsw = 0
  real(dp) :: occf = 2.0_dp                         ! electrons per band state (1: spinors or nspin=2)
  ! --- v7: projections (Mulliken weights) on groups = species, or l if one species
  integer :: ngrp = 1
  integer, allocatable :: grp_of(:)
  character(len=24) :: grp_name(96)
  real(4), allocatable :: Wper(:,:,:,:,:), W3(:,:,:,:,:,:), W1(:,:,:,:), Wpath(:,:,:,:)
  real(dp), allocatable :: pdos(:,:,:)
  logical :: trs = .true.                           ! E(-k) = E(k) holds (false: SOC + magnetism)
  logical :: me_on = .false.
  integer :: zs = 1          ! Bloch vectors kept on every zs-th grid point (memory); nearest stored vertex used
  complex(4), allocatable :: Zg(:,:,:,:,:)          ! (nbas, nqb, nsp_eff, 0:Np-1, 0:Np-1), R gauge
  real(dp), allocatable :: tau_uv(:,:), Esw(:)
  integer, allocatable :: oa1(:), oa2(:)
  real(dp) :: mz_soc = 0.0_dp, mz_soc_dg = -1.0_dp
  integer :: mz_soc_n = 0
  ! spin-orbit: band edges of the states with spin-up (<sz> >= 0) and spin-down (<sz> < 0) character
  ! (1 = up, 2 = down): VBM, CBM, smallest direct gap, <sz> at VBM / CBM; k in 1/Ang
  real(dp) :: sce(3,2) = 0.0_dp, scsz(2,2) = 0.0_dp, sck(3,3,2) = 0.0_dp
  logical :: sc_ok(2) = .false.
  real(dp) :: Qwin(2), cec_dE = 0.0_dp, cec_dk = 0.0_dp, dos_conv = -1.0_dp, ef_conv = 0.0_dp, dos_order = 1.0_dp, dos_diff = -1.0_dp, NelOwn(2) = 0.0_dp, spin_split_max = 0.0_dp

  !------------------------------------------------------------ arguments
  twall0 = wtime()
  fdfname = ''; show = .true.; budget = 240.0_dp; user_qpi = .false.; nqpi_user = 0
  iarg = 0
  do while (iarg < command_argument_count())
     iarg = iarg + 1
     call get_command_argument(iarg, arg)
     if (trim(arg) == '-np' .or. trim(arg) == '--np' .or. trim(arg) == '-n' .or. trim(arg) == '--cores') then
        iarg = iarg + 1; call get_command_argument(iarg, arg); read(arg, *, iostat=est) npreq
        if (est /= 0 .or. npreq < 1) call bad_arg('-np needs a number of cores >= 1, not "'//trim(arg)//'"')
        cycle
     else if (arg(1:4) == '-np=' .or. arg(1:5) == '--np=' .or. arg(1:8) == '--cores=') then
        read(arg(index(arg, '=')+1:), *, iostat=est) npreq
        if (est /= 0 .or. npreq < 1) call bad_arg('-np needs a number of cores >= 1, not "'//trim(arg)//'"')
        cycle
     else if (arg(1:3) == '-np' .and. len_trim(arg) > 3) then
        read(arg(4:), *, iostat=est) npreq
        if (est /= 0 .or. npreq < 1) call bad_arg('-np needs a number of cores >= 1, not "'//trim(arg)//'"')
        cycle
     end if
     if (trim(arg) == '--version' .or. trim(arg) == '-V') then
        write(*,'(a)') 'BandQPI-3D v1.0  ·  Satish Kumar'; stop
     else if (trim(arg) == '--help' .or. trim(arg) == '-h') then
        fdfname = ''; exit
     else if (trim(arg) == '--noshow') then; show = .false.
     else if (trim(arg) == '--fast') then; budget = 60.0_dp
     else if (trim(arg) == '--fine') then; budget = 900.0_dp
     else if (arg(1:9) == '--budget=') then
        read(arg(10:), *, iostat=est) budget
        if (est /= 0 .or. budget <= 0.0_dp) call bad_arg('--budget needs a time in seconds, not "'//trim(arg(10:))//'"')
     else if (arg(1:6) == '--qpi=') then
        user_qpi = .true.; qpi_user = huge(1.0_dp)
        call parse_list(arg(7:), qpi_user, nqpi_user)
     else if (arg(1:2) /= '--') then; fdfname = arg
     else if (arg(1:7) == '--ewin=') then
        ! checked here, before the long computation (the plotter checks it again)
        ewtmp = huge(1.0_dp); call parse_list(arg(8:), ewtmp, nqpi_tmp)
        if (nqpi_tmp /= 2) call bad_arg('--ewin needs two energies E1 < E2 (eV), e.g. --ewin=-3,3')
        if (ewtmp(1) >= ewtmp(2)) call bad_arg('--ewin needs two energies E1 < E2 (eV), e.g. --ewin=-3,3')
        plotargs = trim(plotargs)//' '//trim(arg)
     else if (arg(1:8) == '--smear=' .or. arg(1:13) == '--bandcolors=' .or. &
              arg(1:9) == '--kunits=' .or. arg(1:9) == '--gaptol=') then
        plotargs = trim(plotargs)//' '//trim(arg)    ! figure options, passed to the plotter (checked there)
     else
        call bad_arg('unknown option '//trim(arg))
     end if
  end do
  if (len_trim(fdfname) == 0) then
     write(*,'(a)') ' BandQPI-3D v1.0'
     write(*,'(a)') ' usage: siesta3d.x system.fdf [-np N] [--fast|--fine|--budget=S]'
     write(*,'(a)') '          [--qpi=E1,E2,E3] [--ewin=E1,E2] [--smear=eV]'
     write(*,'(a)') '          [--bandcolors=character|index] [--kunits=bohr|ang]'
     call exit(2)
  end if
  ! folder of the executable (also when started through PATH or a symbolic link): /proc/<pid>/exe
  exedir = ''
  write(envv,'(a,i0,a)') '/proc/', getpid(), '/exe'
  call exe_path(trim(envv), exedir)
  if (len_trim(exedir) == 0) then
     call get_command_argument(0, arg); exedir = dirname(arg)
  end if

  !------------------------------------------------------------ cores (-np N)
  ! OpenMP threads run over k-points; LAPACK/BLAS inside each thread must be serial. Both are fixed
  ! through the environment, which OpenBLAS reads at start-up only -> relaunch this program once.
  ncores = 1
  !$ ncores = omp_get_num_procs()
  call get_environment_variable('BANDQPI_CHILD', envv, status=est)
  if (est /= 0) then
     npuser = (npreq > 0)
     if (npreq <= 0) npreq = ncores
     npreq = min(npreq, ncores)
     write(relaunch,'(a,i0,a)') 'env BANDQPI_CHILD=1 BANDQPI_NCORES=', ncores, ' OMP_NUM_THREADS='
     write(relaunch,'(a,i0,a)') trim(relaunch), npreq, ' OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 ' // &
          'GOTO_NUM_THREADS=1 BLIS_NUM_THREADS=1 OMP_DYNAMIC=false'
     write(relaunch,'(a,i0)') trim(relaunch)//' BANDQPI_NP=', npreq
     if (npuser) relaunch = trim(relaunch)//' BANDQPI_NPUSER=1'
     call get_command_argument(0, arg)
     relaunch = trim(relaunch)//' "'//trim(arg)//'"'
     do iarg = 1, command_argument_count()
        call get_command_argument(iarg, arg)
        relaunch = trim(relaunch)//' "'//trim(arg)//'"'
     end do
     call execute_command_line(trim(relaunch), exitstat=est)
     call exit(est)
  end if
  call get_environment_variable('BANDQPI_NCORES', envv, status=est)
  if (est == 0) read(envv, *, iostat=est) ncores
  call get_environment_variable('BANDQPI_NP', envv, status=est)
  if (est == 0) read(envv, *, iostat=est) npreq
  call get_environment_variable('BANDQPI_NPUSER', envv, status=est)
  npuser = (est == 0)
  colour = isatty(6)

  !------------------------------------------------------------ input
  call fdf_load(trim(fdfname))
  fdfdir = dirname(fdfname)
  slabel = fdf_string('SystemLabel', 'siesta')
  hsxname = join(fdfdir, trim(slabel)//'.HSX')
  call read_hsx(trim(hsxname))
  if (soc) nspin = 1                                ! one spinor channel of dimension 2*no_u
  occf = merge(1.0_dp, 2.0_dp, soc .or. nspin == 2)
  alat = fdf_length('LatticeConstant', 1.0_dp)
  Ef_eV = Ef*Ry2eV; EfOwn = Ef_eV; kT = temp*Ry2eV
  nthreads = 1
  !$ nthreads = omp_get_max_threads()
  outdir = join(fdfdir, trim(slabel)//'_3dbands')
  call execute_command_line('mkdir -p "'//trim(outdir)//'"')
  open(newunit=uv, file=trim(outdir)//'/validation.txt', status='replace')
  write(uv,'(a)') '# status | check | value | note'

  call banner()
  write(*,'(a)') '  ── input ──────────────────────────────────────────────────────────────────'
  write(*,'(a,a)')            '  fdf / HSX                : ', trim(fdfname)//'  /  '//trim(hsxname)
  if (hsx_legacy) write(*,'(a)') '  HSX format               : original (older SIESTA / HONPAS); cell from .XV'
  write(*,'(a,i0,a,i0,a,i0)') '  atoms / orbitals / nspin : ', na_u, ' / ', no_u, ' / ', nspin
  write(*,'(a,3i3)')          '  auxiliary supercell      : ', nsc
  write(*,'(a,f11.5,a,f7.4,a)') '  Fermi energy (SIESTA)    : ', Ef_eV, ' eV   (kT = ', kT, ' eV)'
  write(*,'(a,f10.3)')        '  valence electrons        : ', qtot
  if (soc) write(*,'(a,a,i0,a)') '  spin treatment           : ', merge('spin-orbit (nspin 8)   ', 'non-collinear (nspin 4)', nsraw == 8), &
                                 nbas, ' spinor basis functions'
  write(*,'(a,i0,a,i0,a)', advance='no') '  cores used               : ', nthreads, ' of ', ncores, ' '
  if (npuser .or. nthreads < ncores) then
     write(*,'(a,i0,a)') '  (-np ', nthreads, ')'
  else
     write(*,'(a)') '  (all; limit with  -np N)'
  end if

  call setup_geometry()
  call setup_groups()
  if (dimk == 0) then                              ! 0D: molecule / cluster
     if (user_qpi) write(*,'(a)') '  note: --qpi is ignored for a molecule (no k, no QPI)'
     call molecule_mode()
     close(uv)
     call run_plotter()
     stop
  end if
  call setup_sectors()
  call time_diag()
  call check_time_reversal()
  ewin = 6.0_dp
  call band_window(map_all, nbas, blo, bhi); nbw = bhi - blo + 1
  if (planar) then
     call band_window(map_odd, n_odd, olo, ohi); nbo = ohi - olo + 1
  end if

  call make_path()
  call compute_path()
  call spin_check()
  nsp_eff = merge(1, nspin, spin_degenerate)
  call check_against_siesta_bands()
  call siesta_kgrid_eigs()
  call choose_grids()
  if (dimk >= 2) call setup_qpi_vectors()
  if (dimk >= 2) then
     call ws_polygon()
     call periodic_grid_2d()
  end if
  if (dimk == 3) call periodic_grid_3d()
  if (dimk == 1) call periodic_grid_1d()
  call dos_and_counts()
  call refine_dos()
  if (soc) call soc_moment()
  call band_edges()
  if (soc .and. .not. trs) call soc_spin_edges()
  if (dimk >= 2) then
     call grid_bz()
     call grid_ext()
     call grid_zoom()
     call contours_and_qpi()
  end if
  call validation_checks()
  call write_summary()
  close(uv)
  call run_plotter()

contains

  !==================================================================
  ! spin-orbit, magnetic: S_z is not conserved, but each spinor state has a definite spin CHARACTER
  ! <sz> = sum of the up-spinor minus the down-spinor Mulliken weights.  States with <sz> >= 0 are
  ! "spin up", the others "spin down".  Band edges and gaps of each character are taken over the full
  ! periodic k-grid (occupied = band index <= nocc), so the spin-up and spin-down gaps are separate.
  subroutine soc_spin_edges()
    integer :: a, b, c, j, n1, n2, n3, ngh, ch, jg
    real(dp) :: f(3), e, sz, tot, vb(2), cb(2), dgk
    real(4) :: wl(ngrp)
    real(dp) :: ev, ec
    ngh = ngrp/2
    sce(1,:) = -huge(1.0_dp); sce(2,:) = huge(1.0_dp); sce(3,:) = huge(1.0_dp)
    n1 = 1; n2 = 1; n3 = 1
    select case (dimk)
    case (1); n1 = N1
    case (2); n1 = Na; n2 = Nb
    case (3); n1 = N3; n2 = N3; n3 = N3
    case default; return
    end select
    do c = 0, n3-1
    do b = 0, n2-1
    do a = 0, n1-1
       f = [real(a,dp)/n1, real(b,dp)/n2, real(c,dp)/n3]
       vb = -huge(1.0_dp); cb = huge(1.0_dp)
       do j = 1, nbw
          jg = blo + j - 1
          select case (dimk)
          case (1); e = E1(j,1,a); wl = W1(:,j,1,a)
          case (2); e = Eper(j,1,a,b); wl = Wper(:,j,1,a,b)
          case default; e = E3(j,1,a,b,c); wl = W3(:,j,1,a,b,c)
          end select
          tot = sum(real(wl,dp)); if (tot <= 0.0_dp) cycle
          sz = (sum(real(wl(1:ngh),dp)) - sum(real(wl(ngh+1:ngrp),dp)))/tot
          ch = merge(1, 2, sz >= 0.0_dp)
          if (jg <= nocc(1)) then
             vb(ch) = max(vb(ch), e)
             if (e > sce(1,ch)) then
                sce(1,ch) = e; scsz(1,ch) = sz; sck(:,1,ch) = fold3(kf(f(1)*Bp(:,1) + f(2)*Bp(:,2) + f(3)*Bp(:,3)))
             end if
          else
             cb(ch) = min(cb(ch), e)
             if (e < sce(2,ch)) then
                sce(2,ch) = e; scsz(2,ch) = sz; sck(:,2,ch) = fold3(kf(f(1)*Bp(:,1) + f(2)*Bp(:,2) + f(3)*Bp(:,3)))
             end if
          end if
       end do
       do ch = 1, 2
          if (vb(ch) > -huge(1.0_dp) .and. cb(ch) < huge(1.0_dp)) then
             dgk = cb(ch) - vb(ch)
             if (dgk < sce(3,ch)) then
                sce(3,ch) = dgk; sck(:,3,ch) = fold3(kf(f(1)*Bp(:,1) + f(2)*Bp(:,2) + f(3)*Bp(:,3)))
             end if
          end if
       end do
    end do
    end do
    end do
    do ch = 1, 2
       sc_ok(ch) = sce(1,ch) > -huge(1.0_dp) .and. sce(2,ch) < huge(1.0_dp)
       if (.not. sc_ok(ch)) cycle
       ev = sce(1,ch) - Ef_eV; ec = sce(2,ch) - Ef_eV
       write(*,'(a,a,f8.4,a,f8.4,a,f7.4,a)') merge('  spin-up states (SOC)     ', '  spin-down states (SOC)   ', ch == 1), &
            ': VBM', ev, '  CBM', ec, '  gap', sce(2,ch) - sce(1,ch), ' eV'
    end do
  end subroutine soc_spin_edges

  !==================================================================
  ! spin-orbit runs: spin moment m_z = N_up - N_dn of ALL occupied states (not only the plotted window),
  ! Mulliken spin-z weights summed over the up / down spinor groups, on two uniform k-grids (n/2 and n);
  ! the difference between the grids N, N/2, N/4 gives the k-grid error (Richardson; limit 1 %)
  subroutine soc_moment()
    integer :: it, n, nk, q, a, b, c, iu, nd, ngh, j
    real(dp) :: mz(2), tt, k(3), acc
    real(dp), allocatable :: e(:)
    real(4), allocatable :: wg(:,:)
    if (trs) then
       mz_soc = 0.0_dp; mz_soc_dg = 0.0_dp; mz_soc_n = 0
       write(*,'(a)') '  spin moment (SOC)        : 0 by time-reversal symmetry'
       return
    end if
    tt = wtime()
    ngh = ngrp/2
    iu = min(nbas, maxval(nocc(1:nspin)) + 8)
    do it = 1, 2
       select case (dimk)
       case (1); n = merge(16, 32, it == 1); nd = 1
       case (2); n = merge(8, 16, it == 1); nd = 2
       case (3); n = merge(4, 8, it == 1); nd = 3
       case default; n = 1; nd = 0
       end select
       nk = n**nd
       acc = 0.0_dp
       !$omp parallel do schedule(dynamic) private(a,b,c,k,e,wg,j) reduction(+:acc)
       do q = 0, nk - 1
          a = mod(q, n); b = mod(q/n, n); c = q/(n*n)
          k = (real(a,dp)/n)*B1
          if (nd >= 2) k = k + (real(b,dp)/n)*B2
          if (nd == 3) k = k + (real(c,dp)/n)*Bp(:,3)
          if (.not. allocated(e)) allocate(e(iu), wg(ngrp, iu))
          call eigk_w(k, 1, 1, iu, e, 1, iu, wg)
          do j = 1, iu
             if (e(j) <= EfOwn) acc = acc + sum(real(wg(1:ngh, j), dp)) - sum(real(wg(ngh+1:ngrp, j), dp))
          end do
       end do
       !$omp end parallel do
       mz(it) = acc/nk
       mz_soc_n = n
    end do
    mz_soc = mz(2); mz_soc_dg = abs(mz(2) - mz(1))
    write(*,'(a,f9.4,a,f7.4,a,i0,a,i0,a,f6.1,a)') '  spin moment m_z (SOC)    : ', mz_soc, ' muB  (n vs n/2: ', &
         mz_soc_dg, '; ', n, '^', nd, ' k; ', wtime()-tt, ' s)'
  end subroutine soc_moment

  !==================================================================
  subroutine print_list(head, items)
    ! "head item item ..." wrapped at 78 columns, continuation lines indented under the first item
    character(len=*), intent(in) :: head, items(:)
    character(len=1024) :: line
    integer :: i
    line = head
    do i = 1, size(items)
       if (len_trim(line) + 1 + len_trim(items(i)) > 78 .and. len_trim(line) > len(head)) then
          write(*,'(a)') trim(line); line = repeat(' ', len(head))
       end if
       if (len_trim(line) < len(head)) then
          line = line(1:len(head))//trim(items(i))
       else
          line = trim(line)//' '//trim(items(i))
       end if
    end do
    write(*,'(a)') trim(line)
  end subroutine print_list

  !==================================================================
  ! precision target: the estimated DOS error of the k-grid N (Richardson from N, N/2, N/4) below 1 %.
  ! If not, the periodic grid is refined by the factor that the observed error law asks for (at least x1.26,
  ! at most x2 per direction; 1D x2 or more) and everything that depends on it is recomputed, within a wall-time cap.
  subroutine refine_dos()
    real(dp), parameter :: tol = 0.01_dp
    real(dp) :: per_k, est, tlim, mem, f
    integer :: it, nN1, nNa, nNb, nN3
    if (dos_conv < 0.0_dp .or. dos_conv < tol) return
    per_k = tdiag*nsp_eff/real(nthreads,dp)
    tlim = max(4.0_dp*budget, 900.0_dp)
    do it = 1, 5
       if (dos_conv < tol) exit
       write(*,'(a,f6.2,a)') '  DOS k-grid error         : ', 100*dos_conv, ' %  > 1 %  ->  refining the k-grid'
       ! the error falls ~ h^p with the observed order p (2 for smooth DOS, less near 1D van Hove peaks)
       f = min(2.0_dp, max(1.26_dp, 1.1_dp*(dos_conv/tol)**(1.0_dp/max(dos_order, 0.5_dp))))
       select case (dimk)
       case (1)
          nN1 = 4*nint(max(2.0_dp, f)*N1/4); est = per_k*nN1/2
       case (2)
          nNa = 4*nint(f*Na/4); nNb = 4*nint(f*Nb/4); est = per_k*real(nNa,dp)*nNb/2
          if (me_on) then
             if (pick_zs(nNa, nNb) == 0) then
                write(*,'(a)') '  DOS refinement           : stopped (memory of the Bloch-QPI vectors); use --fine'
                exit
             end if
          end if
       case (3)
          nN3 = 4*nint(f*N3/4); if (nN3 <= N3) nN3 = N3 + 4
          est = per_k*real(nN3,dp)**3/2
       case default
          exit
       end select
       if (wtime() - twall0 + est > tlim) then
          write(*,'(a,f7.0,a)') '  DOS refinement           : stopped (needs ~', est, ' s more); use --budget or --fine'
          exit
       end if
       deallocate(egrid, dos, pdos)
       if (allocated(dos_o)) deallocate(dos_o)
       select case (dimk)
       case (1)
          deallocate(E1, W1); N1 = nN1
          write(*,'(a,i0)') '  1D grid (refined)        : ', N1
          call periodic_grid_1d()
       case (2)
          deallocate(Eper, Wper); if (allocated(Eoper)) deallocate(Eoper)
          Na = nNa; Nb = nNb
          if (me_on) then
             zs = pick_zs(Na, Nb)
             deallocate(Zg); allocate(Zg(nbas, nqb, nsp_eff, 0:Na/zs-1, 0:Nb/zs-1))
          end if
          write(*,'(a,i0,a,i0)') '  plane grid (refined)     : ', Na, ' x ', Nb
          call periodic_grid_2d()
       case (3)
          deallocate(E3, W3); N3 = nN3
          write(*,'(a,i0,a)') '  3D grid (refined)        : ', N3, '^3'
          call periodic_grid_3d()
       end select
       call dos_and_counts()
    end do
    if (dos_conv < tol) write(*,'(a,f6.2,a)') '  DOS k-grid error         : ', 100*dos_conv, ' %  (target < 1 %)'
  end subroutine refine_dos

  !==================================================================
  subroutine banner()
    ! one box, 79 columns wide: version in the top border; left the emblem (probability density of the
    ! hydrogen-like 3d orbital (n,l,m) = (3,2,1), inferno colours, two pixels per character cell),
    ! right the name, purpose and author. Without colour the emblem is drawn with shade characters.
    integer, parameter :: EMB(18,19) = reshape([ &
         16,16,232,232,232,232,232,232,232,16,232,232,232,232,232,232,232,16,16,232,232,232,233,233,233,233,232, &
         16,232,233,233,233,233,232,232,232,232,232,232,233,234,17,17,17,233,16,233,17,17,17,234,233,232,232,232, &
         232,233,235,53,53,53,53,17,16,17,53,53,53,53,235,233,232,232,233,234,53,53,89,125,89,53,16,53,89,125,89, &
         53,53,234,233,232,233,17,53,89,131,202,202,125,16,125,202,202,131,89,53,17,233,232,233,17,53,125,202,221, &
         229,208,16,208,229,221,202,125,53,17,233,232,233,17,53,89,202,229,229,229,16,229,229,229,202,89,53,17, &
         233,232,232,233,17,53,125,208,229,229,16,229,229,208,125,53,17,233,232,16,16,16,16,16,16,16,16,16,16,16, &
         16,16,16,16,16,16,16,232,232,233,17,53,125,208,229,229,16,229,229,208,125,53,17,233,232,232,233,17,53, &
         89,202,229,229,229,16,229,229,229,202,89,53,17,233,232,233,17,53,125,202,221,229,208,16,208,229,221,202, &
         125,53,17,233,232,233,17,53,89,131,202,202,125,16,125,202,202,131,89,53,17,233,232,233,234,53,53,89,125, &
         89,53,16,53,89,125,89,53,53,234,233,232,232,233,235,53,53,53,53,17,16,17,53,53,53,53,235,233,232,232,232, &
         232,233,234,17,17,17,233,16,233,17,17,17,234,233,232,232,16,232,232,232,233,233,233,233,232,16,232,233, &
         233,233,233,232,232,232,16,16,232,232,232,232,232,232,232,16,232,232,232,232,232,232,232,16], [18, 19])
    integer, parameter :: EMS(9,19) = reshape([ &
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1,0,1,1,1,0,0,1,1,1,0,1,1,1,1,0,1,2,3,1,2, &
         3,1,0,0,1,2,4,2,4,3,1,0,0,1,2,4,3,4,4,1,0,0,0,1,4,3,4,2,1,0,0,0,0,0,0,0,0,0,0,0,0,1,4,3,4,2,1,0,0,1,2, &
         4,3,4,4,1,0,0,1,2,4,2,4,3,1,0,0,1,2,3,1,2,3,1,0,0,1,1,1,0,1,1,1,1,0,1,1,1,0,1,1,1,0,0,0,0,0,0,0,0,1,0, &
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], [9, 19])
    character(len=60) :: txt(9)
    character(len=8) :: col(9)
    logical :: bold(9)
    character(len=*), parameter :: ver = ' BandQPI-3D · version 1.0 '
    character(len=*), parameter :: shd(0:4) = [character(len=3) :: ' ', '░', '▒', '▓', '█']
    integer, parameter :: W = 77
    integer :: n, r, c, pr
    character(len=8) :: a1, a2
    txt = ''; col = ''; bold = .false.
    txt(2) = 'BandQPI-3D';                                   col(2) = '214'; bold(2) = .true.
    txt(3) = 'version 1.0';                                  col(3) = '250'
    txt(5) = 'bands · DOS · quasiparticle interference';     col(5) = '117'
    txt(6) = 'SIESTA / HONPAS Hamiltonians · 0D to 3D';      col(6) = '117'
    txt(8) = 'Satish Kumar · Computational Condensed Matter Physics'; bold(8) = .true.
    txt(9) = '97satish02@gmail.com | ORCID: 0009-0000-2018-1564'
    n = (W - dwidth(ver))/2
    write(*,'(a)') ''
    call bedge('╭', '╮', n, ver, W - dwidth(ver) - n)
    do r = 1, 9
       if (colour) write(*,'(a)', advance='no') achar(27)//'[38;5;67m'
       write(*,'(a)', advance='no') ' │ '
       do c = 1, 19
          if (colour) then
             write(a1,'(i0)') EMB(2*r-1, c); write(a2,'(i0)') EMB(2*r, c)
             write(*,'(a)', advance='no') achar(27)//'[38;5;'//trim(a1)//';48;5;'//trim(a2)//'m▀'
          else
             write(*,'(a)', advance='no') trim(shd(EMS(r, c)))
             if (EMS(r, c) == 0) write(*,'(a)', advance='no') ' '
          end if
       end do
       if (colour) write(*,'(a)', advance='no') achar(27)//'[0m'
       write(*,'(a)', advance='no') '  '
       if (colour .and. len_trim(col(r)) > 0) then
          if (bold(r)) then
             write(*,'(a)', advance='no') achar(27)//'[1;38;5;'//trim(col(r))//'m'
          else
             write(*,'(a)', advance='no') achar(27)//'[38;5;'//trim(col(r))//'m'
          end if
       else if (colour .and. bold(r)) then
          write(*,'(a)', advance='no') achar(27)//'[1m'
       end if
       write(*,'(a)', advance='no') trim(txt(r))
       pr = W - 1 - 19 - 2 - dwidth(trim(txt(r)))
       if (colour) write(*,'(a)', advance='no') achar(27)//'[0;38;5;67m'
       write(*,'(a)', advance='no') repeat(' ', max(pr, 0))//'│'
       if (colour) write(*,'(a)', advance='no') achar(27)//'[0m'
       write(*,*)
    end do
    call bedge('╰', '╯', W, '', 0)
    write(*,'(a)') ''
  end subroutine banner

  integer function dwidth(t)                 ! display width of a UTF-8 string (continuation bytes skipped)
    character(len=*), intent(in) :: t
    integer :: i
    dwidth = 0
    do i = 1, len(t)
       if (iachar(t(i:i)) < 128 .or. iachar(t(i:i)) >= 192) dwidth = dwidth + 1
    end do
  end function dwidth

  subroutine bedge(cl, cr, n1, mid, n2)
    character(len=*), intent(in) :: cl, cr, mid
    integer, intent(in) :: n1, n2
    integer :: i
    if (colour) write(*,'(a)', advance='no') achar(27)//'[38;5;67m'
    write(*,'(a)', advance='no') ' '//cl
    do i = 1, n1; write(*,'(a)', advance='no') '─'; end do
    if (len(mid) > 0) then
       if (colour) write(*,'(a)', advance='no') achar(27)//'[1;38;5;214m'
       write(*,'(a)', advance='no') mid
       if (colour) write(*,'(a)', advance='no') achar(27)//'[0;38;5;67m'
    end if
    do i = 1, n2; write(*,'(a)', advance='no') '─'; end do
    write(*,'(a)', advance='no') cr
    if (colour) write(*,'(a)', advance='no') achar(27)//'[0m'
    write(*,*)
  end subroutine bedge

  subroutine bline(txt, col, bold)
    ! one box row: txt centred in the 77 inner columns
    character(len=*), intent(in) :: txt, col
    logical, intent(in) :: bold
    integer :: pl, pr
    pl = (77 - dwidth(txt))/2; pr = 77 - dwidth(txt) - pl
    if (colour) write(*,'(a)', advance='no') achar(27)//'[38;5;67m'
    write(*,'(a)', advance='no') ' │'
    if (colour) write(*,'(a)', advance='no') achar(27)//'[0m'
    write(*,'(a)', advance='no') repeat(' ', pl)
    if (colour .and. len(col) > 0 .and. bold) write(*,'(a)', advance='no') achar(27)//'[1;38;5;'//col//'m'
    if (colour .and. len(col) > 0 .and. .not. bold) write(*,'(a)', advance='no') achar(27)//'[38;5;'//col//'m'
    if (colour .and. len(col) == 0 .and. bold) write(*,'(a)', advance='no') achar(27)//'[1m'
    write(*,'(a)', advance='no') txt
    if (colour) write(*,'(a)', advance='no') achar(27)//'[0;38;5;67m'
    write(*,'(a)', advance='no') repeat(' ', pr)//'│'
    if (colour) write(*,'(a)', advance='no') achar(27)//'[0m'
    write(*,*)
  end subroutine bline

  !==================================================================
  function katxt(kv) result(t)
    ! "at k = (...) 1/Ang"; for a wire the single coordinate along the axis and its fraction of |b|
    real(dp), intent(in) :: kv(3)
    character(len=120) :: t
    real(dp) :: bl
    if (dimk == 1) then
       bl = norm2(B1) / Bohr2A
       write(t,'(a,f8.4,a,f5.3,a)') '  at k =', kv(1), ' 1/Ang  (', abs(kv(1)) / max(bl, 1.0e-12_dp), ' b)'
    else
       write(t,'(a,3f8.4,a)') '  k = (', kv, ') 1/Ang'
    end if
  end function katxt
  !==================================================================
  real(dp) function wtime()
    integer(8) :: c, r
    call system_clock(c, r); wtime = real(c,dp)/real(r,dp)
  end function wtime

  subroutine parse_list(s, x, n)
    character(len=*), intent(in) :: s
    real(dp), intent(inout) :: x(3)
    integer, intent(out) :: n
    integer :: a, b, ios
    n = 0; a = 1
    do while (a <= len_trim(s) .and. n < 3)
       b = index(s(a:), ',')
       if (b == 0) then; b = len_trim(s) + 1; else; b = a + b - 1; end if
       n = n + 1; read(s(a:b-1), *, iostat=ios) x(n)
       if (ios /= 0) n = n - 1
       a = b + 1
    end do
  end subroutine parse_list

  ! plane-frame coordinates (u,v,n) in 1/Ang  ->  Cartesian k in 1/Bohr
  function kc(k) result(kk)
    real(dp), intent(in) :: k(3)
    real(dp) :: kk(3)
    kk = (k(1)*eu + k(2)*ev + k(3)*en)*Bohr2A
  end function kc
  function kf(kk) result(k)                    ! Cartesian 1/Bohr -> frame 1/Ang
    real(dp), intent(in) :: kk(3)
    real(dp) :: k(3)
    k = [dot_product(kk,eu), dot_product(kk,ev), dot_product(kk,en)]/Bohr2A
  end function kf

  subroutine vcheck(stat, name, val, note)
    character(len=*), intent(in) :: stat, name, note
    real(dp), intent(in) :: val
    write(uv,'(a,a,a,a,es11.3,a,a)') trim(stat), ' | ', trim(name), ' | ', val, ' | ', trim(note)
    flush(uv)
  end subroutine vcheck

  subroutine progress(i, n, txt)
    integer, intent(in) :: i, n
    character(len=*), intent(in) :: txt
    integer :: p
    p = int(100.0_dp*i/max(n,1))
    if (p >= prog_last + 10 .and. p < 100) then
       prog_last = (p/10)*10
       write(*,'(a,a,a,a,i3,a,a)', advance='no') achar(13), '  ', pad25(txt), ': ', prog_last, ' %  ', &
            repeat('#', prog_last/10)//repeat('.', 10 - prog_last/10)
       flush(6)
    end if
  end subroutine progress
  function pad25(t) result(r)
    character(len=*), intent(in) :: t
    character(len=25) :: r
    r = adjustl(t)
  end function pad25
  subroutine progress_done(txt, tt)
    character(len=*), intent(in) :: txt
    real(dp), intent(in) :: tt
    write(*,'(a,a,a,a,f8.1,a)') achar(13), '  ', pad25(txt), ': ', wtime()-tt, ' s                 '
    prog_last = -10
  end subroutine progress_done

  !==================================================================
  subroutine setup_geometry()
    integer :: n
    real(dp) :: z(3)
    per = 0
    do i = 1, 3
       if (nsc(i) > 1 .or. kcell(i,i) > 1) per(i) = 1
    end do
    dimk = sum(per)
    if (dimk == 0) then
       write(*,'(a)') '  periodic dimensions      : 0D   (molecule / cluster: levels at Gamma)'
       B1 = 0.0_dp; B2 = 0.0_dp; b1p = 0.0_dp; b2p = 0.0_dp
       return
    end if
    n = 0; ip = 0
    do i = 1, 3
       if (per(i) == 1) then; n = n + 1; ip(n) = i; end if
    end do
    Bp = 0.0_dp
    do i = 1, dimk
       Bp(:,i) = rcell(:,ip(i))
    end do
    B1 = Bp(:,1); B2 = Bp(:,2)
    z = [0.0_dp, 0.0_dp, 1.0_dp]
    if (dimk >= 2) then
       en = cross(B1, B2); en = en/norm2(en)
       if (dot_product(en, z) < 0.0_dp) en = -en
       eu = [1.0_dp,0.0_dp,0.0_dp] - en(1)*en
       if (norm2(eu) < 1e-6_dp) eu = [0.0_dp,1.0_dp,0.0_dp] - en(2)*en
       eu = eu/norm2(eu); ev = cross(en, eu)
    else
       eu = B1/norm2(B1)
       en = z - dot_product(z,eu)*eu
       if (norm2(en) < 1e-6_dp) en = [1.0_dp,0.0_dp,0.0_dp] - eu(1)*eu
       en = en/norm2(en); ev = cross(en, eu)
    end if
    b1p = [dot_product(B1,eu), dot_product(B1,ev)]/Bohr2A
    b2p = [dot_product(B2,eu), dot_product(B2,ev)]/Bohr2A
    write(*,'(a,i0,a)', advance='no') '  periodic dimensions      : ', dimk, 'D'
    if (dimk == 3) write(*,'(a)', advance='no') '   (bulk crystal: 3D DOS and gap)'
    if (dimk == 1) write(*,'(a)', advance='no') '   (wire / chain)'
    write(*,*)
  end subroutine setup_geometry

  ! mirror (z -> -z) sectors for a flat sheet in the xy plane
  subroutine setup_sectors()
    integer :: io, ind
    real(dp) :: zspread, cmax
    logical :: oa, ob
    allocate(map_all(no_u), map_odd(no_u))
    do io = 1, no_u
       map_all(io) = io
    end do
    map_odd = 0; n_odd = 0; planar = .false.
    if (soc) then
       write(*,'(a)') '  mirror analysis          : spinor bands (spin-orbit) -> bands by index'
       return
    end if
    if (dimk /= 2 .or. any(orb_l < 0)) return
    if (abs(abs(en(3)) - 1.0_dp) > 1e-6_dp) then
       write(*,'(a)') '  mirror analysis          : sheet not in the xy plane -> bands by index'
       return
    end if
    zspread = maxval(xa(3,:)) - minval(xa(3,:))
    if (zspread > 2e-3_dp) then
       write(*,'(a,f7.3,a)') '  mirror analysis          : not flat (dz = ', zspread*Bohr2A, &
            ' Ang) -> bands by index'
       return
    end if
    cmax = 0.0_dp
    do io = 1, no_u
       oa = mod(orb_l(io) + abs(orb_m(io)), 2) == 1
       do ind = listhptr(io)+1, listhptr(io)+numh(io)
          ob = mod(orb_l(col(ind)) + abs(orb_m(col(ind))), 2) == 1
          if (oa .neqv. ob) cmax = max(cmax, maxval(abs(H(ind,:))), abs(S(ind)))
       end do
    end do
    do io = 1, no_u
       if (mod(orb_l(io) + abs(orb_m(io)), 2) == 1) then
          n_odd = n_odd + 1; map_odd(io) = n_odd
       end if
    end do
    if (cmax < 1e-6_dp .and. n_odd > 0) planar = .true.
    write(*,'(a,i0,a,es9.2,a)') '  mirror-odd (pi) orbitals : ', n_odd, '   (even-odd coupling ', cmax*Ry2eV, ' eV)'
    if (planar) call vcheck('PASS', 'mirror symmetry: even-odd coupling (eV)', cmax*Ry2eV, &
                            'pi and sigma bands separate exactly')
  end subroutine setup_sectors

  ! projection groups: one group per atomic shell of the SIESTA basis, (species, n, l) -> "Mo-4d", "S-3p",
  ! named as SIESTA lists them in <label>.ORB_INDX (polarisation shells keep SIESTA's n, e.g. C-3d).
  ! The plotter sums shells to elements or to s/p/d/f exactly. More than 48 shells: shells of equal
  ! (species, l) merged ("Mo-d"); still more: species only. Spin-orbit: every group split into its
  ! spinor components _up / _dn.
  subroutine setup_groups()
    integer :: b, io, is, mode, k, n, key(no_u)
    integer, parameter :: maxg = 48
    character(len=1), parameter :: lname(0:3) = ['s','p','d','f']
    integer :: gs(maxg*2), gn(maxg*2), gl(maxg*2), perm(maxg*2), inv(maxg*2)
    allocate(grp_of(nbas))
    if (.not. all(orb_l >= 0)) then
       ngrp = min(nspecies, maxg)
       do is = 1, ngrp
          grp_name(is) = adjustl(spec_label(is))
       end do
       do b = 1, nbas
          io = merge((b+1)/2, b, soc)
          grp_of(b) = min(isa(orb_atom(io)), maxg)
       end do
    else
       do mode = 1, 3
          ngrp = 0
          do io = 1, no_u
             is = isa(orb_atom(io))
             n = merge(orb_n(io), 0, mode == 1)
             k = 0
             do b = 1, ngrp
                if (gs(b) == is .and. (mode == 3 .or. (gl(b) == min(orb_l(io), 3) .and. gn(b) == n))) then
                   k = b; exit
                end if
             end do
             if (k == 0) then
                if (ngrp == maxg*2) exit
                ngrp = ngrp + 1; k = ngrp
                gs(k) = is; gl(k) = min(orb_l(io), 3); gn(k) = n
             end if
             key(io) = k
          end do
          if (ngrp <= maxg) exit
       end do
       ! order: species as in the fdf, then n, then l (Cr-3s Cr-3p Cr-3d Cr-4s ...): insertion sort + remap
       do b = 1, ngrp
          perm(b) = b
       end do
       do b = 2, ngrp
          k = perm(b); n = b - 1
          do while (n >= 1)
             if (gs(perm(n))*10000 + gn(perm(n))*10 + gl(perm(n)) <= gs(k)*10000 + gn(k)*10 + gl(k)) exit
             perm(n+1) = perm(n); n = n - 1
          end do
          perm(n+1) = k
       end do
       do b = 1, ngrp
          inv(perm(b)) = b
       end do
       gs(1:ngrp) = gs(perm(1:ngrp)); gn(1:ngrp) = gn(perm(1:ngrp)); gl(1:ngrp) = gl(perm(1:ngrp))
       do io = 1, no_u
          key(io) = inv(key(io))
       end do
       do b = 1, ngrp
          if (mode == 3) then
             grp_name(b) = trim(adjustl(spec_label(gs(b))))
          else if (gn(b) > 0) then
             write(grp_name(b), '(a,a,i0,a)') trim(adjustl(spec_label(gs(b)))), '-', gn(b), lname(gl(b))
          else
             grp_name(b) = trim(adjustl(spec_label(gs(b))))//'-'//lname(gl(b))
          end if
       end do
       do b = 1, nbas
          io = merge((b+1)/2, b, soc)
          grp_of(b) = key(io)
       end do
    end if
    if (soc) then
       ! spinors: resolve every group into its spin-up and spin-down components
       ! (basis 2*io-1 = up, 2*io = down), as in SIESTA's spin-orbit PDOS
       do is = 1, ngrp
          grp_name(ngrp+is) = trim(grp_name(is))//'_dn'
          grp_name(is) = trim(grp_name(is))//'_up'
       end do
       do b = 1, nbas
          if (mod(b, 2) == 0) grp_of(b) = grp_of(b) + ngrp
       end do
       ngrp = 2*ngrp
    end if
    call print_list('  projections              : ', grp_name(1:ngrp))
  end subroutine setup_groups

  ! eigenvalues il..iu, Mulliken group weights of bands wl..wu
  !   w_g(n) = sum_{b in g} Re[ c_b^* (S(k) c)_b ],  sum_g w_g = 1
  ! and (optionally) the eigenvectors il..iu (R gauge, c^+ S c = 1)
  subroutine eigk_w(k, ss, il, iu, e, wl, wu, wg, z)
    real(dp), intent(in) :: k(3)
    integer, intent(in) :: ss, il, iu, wl, wu
    real(dp), intent(out) :: e(iu-il+1)
    real(4), intent(out) :: wg(ngrp, wl:wu)
    complex(dp), intent(out), optional :: z(nbas, iu-il+1)
    complex(dp), allocatable :: Hk(:,:), Sk(:,:), S0(:,:), work(:), zz(:,:), Sz(:,:)
    real(dp), allocatable :: rwork(:), w(:)
    integer, allocatable :: iwork(:), ifail(:)
    integer :: info, m, nb, lw, n, b
    real(dp) :: acc(ngrp)
    nb = nbas; lw = 34*nb
    allocate(Hk(nb,nb), Sk(nb,nb), S0(nb,nb), w(nb), rwork(7*nb), iwork(5*nb), ifail(nb), work(lw), &
             zz(nb, iu-il+1), Sz(nb, wu-wl+1))
    call build_hk(k, ss, map_all, nb, Hk, Sk)
    if (soc) then
       S0 = Sk
    else                                        ! build_hk fills the upper triangle only
       S0 = Sk
       do b = 1, nb
          do n = b+1, nb
             S0(n,b) = conjg(Sk(b,n))
          end do
       end do
    end if
    call zhegvx(1, 'V', 'I', 'U', nb, Hk, nb, Sk, nb, 0.0_dp, 0.0_dp, il, iu, &
                0.0_dp, m, w, zz, nb, work, lw, rwork, iwork, ifail, info)
    if (info /= 0) then
       write(*,*) 'LAPACK info =', info, '  k (1/Bohr) =', k
       stop ' diagonalisation with eigenvectors failed'
    end if
    e = w(1:iu-il+1)*Ry2eV
    call zgemm('N', 'N', nb, wu-wl+1, nb, (1.0_dp,0.0_dp), S0, nb, zz(:,wl-il+1:wu-il+1), nb, &
               (0.0_dp,0.0_dp), Sz, nb)
    do n = wl, wu
       acc = 0.0_dp
       do b = 1, nb
          acc(grp_of(b)) = acc(grp_of(b)) + real(conjg(zz(b,n-il+1))*Sz(b,n-wl+1))
       end do
       wg(:,n) = real(acc, 4)
    end do
    if (present(z)) z = zz
  end subroutine eigk_w

  ! with spin-orbit coupling AND a magnetic moment, E(-k) /= E(k): then all
  ! grids are computed in full instead of using the half-grid shortcut
  subroutine check_time_reversal()
    real(dp) :: e1(nbas), e2(nbas), kr(3), d
    integer :: t
    real(dp), parameter :: rr(3,2) = reshape([0.1734_dp,0.2917_dp,0.1123_dp, 0.3812_dp,0.0791_dp,0.2333_dp], [3,2])
    trs = .true.
    if (.not. soc) return
    d = 0.0_dp
    do t = 1, 2
       kr = rr(1,t)*Bp(:,1) + rr(2,t)*Bp(:,2) + rr(3,t)*Bp(:,3)
       call eigk(kr, 1, map_all, nbas, 1, nbas, e1)
       call eigk(-kr, 1, map_all, nbas, 1, nbas, e2)
       d = max(d, maxval(abs(e1 - e2)))
    end do
    trs = d < 1e-5_dp
    if (trs) then
       write(*,'(a)') '  time reversal            : E(-k) = E(k) (Kramers pairs) -> half grids'
    else
       write(*,'(a,es9.2,a)') '  time reversal            : broken, max|E(k)-E(-k)| =', d, ' eV -> full grids'
    end if
  end subroutine check_time_reversal

  subroutine time_diag()
    real(dp) :: et(nbas)
    integer :: nrep
    call eigk(0.137_dp*B1 + 0.071_dp*Bp(:,2) + 0.03_dp*Bp(:,3), 1, map_all, nbas, 1, nbas, et)
    t0 = wtime(); nrep = 0
    do while (nrep < 3 .or. (wtime() - t0 < 0.3_dp .and. nrep < 2000))
       call eigk((0.137_dp+0.01_dp*nrep)*B1 + 0.071_dp*Bp(:,2), 1, map_all, nbas, 1, nbas, et)
       nrep = nrep + 1
    end do
    tdiag = (wtime() - t0)/nrep
    write(*,'(/a)') '  ── calculation ────────────────────────────────────────────────────────────'
    write(*,'(a,f9.4,a)') '  time per diagonalisation : ', tdiag, ' s'
  end subroutine time_diag

  !==================================================================
  subroutine band_window(map, nb, lo, hi)
    integer, intent(in) :: map(nbas), nb
    integer, intent(out) :: lo, hi
    real(dp) :: emn(nb), emx(nb), e(nb), f(3)
    integer :: a, b, c, ss, n6, nb2, nc
    emn = huge(1.0_dp); emx = -huge(1.0_dp); n6 = 6
    nb2 = merge(n6-1, 0, dimk >= 2); nc = merge(3, 0, dimk == 3)
    do ss = 1, nspin
       do a = 0, n6-1
          do b = 0, nb2
             do c = 0, nc
                f = [real(a,dp)/n6, real(b,dp)/n6, real(c,dp)/4]
                call eigk(f(1)*Bp(:,1) + f(2)*Bp(:,2) + f(3)*Bp(:,3), ss, map, nb, 1, nb, e)
                emn = min(emn, e); emx = max(emx, e)
             end do
          end do
       end do
    end do
    lo = nb; hi = 1
    do a = 1, nb
       if (emx(a) >= Ef_eV - ewin .and. emn(a) <= Ef_eV + ewin) then
          lo = min(lo, a); hi = max(hi, a)
       end if
    end do
    if (hi < lo) then; lo = 1; hi = min(nb, 4); end if
    lo = max(1, lo - 1); hi = min(nb, hi + 1)
  end subroutine band_window

  !==================================================================
  subroutine make_path()
    integer :: i1, i2, l, np, iseg, ios2, nseg, ntot
    real(dp) :: f(3), lenp, dl
    character(len=256) :: lab
    integer, allocatable :: nps(:)
    nnodes = 0
    call fdf_block('BandLines', i1, i2)
    bscale = norm(fdf_string('BandLinesScale', 'pi/a'))
    if (i1 > 0 .and. i2 >= i1) then
       do l = i1, i2
          read(lines(l), *, iostat=ios2) np, f
          if (ios2 /= 0 .or. nnodes >= 60) cycle
          nnodes = nnodes + 1
          lab = word(lines(l), 5)
          if (norm(lab) == '\gamma' .or. norm(lab) == 'gamma' .or. norm(lab) == 'g') lab = 'G'
          if (bscale == 'reciprocallatticevectors') then
             knode(:,nnodes) = f(1)*rcell(:,1) + f(2)*rcell(:,2) + f(3)*rcell(:,3)
          else
             knode(:,nnodes) = f*pi/alat
          end if
          if (len_trim(lab) == 0) then                 ! unnamed node: Gamma, or X / Y / Z at a zone face
             lab = '.'
             if (maxval(abs(knode(:,nnodes))) < 1.0e-8_dp) then
                lab = 'G'
             else if (bscale == 'reciprocallatticevectors' .and. count(abs(f) > 1.0e-8_dp) == 1 .and. &
                      abs(maxval(abs(f)) - 0.5_dp) < 1.0e-6_dp) then
                lab = 'XYZ'(maxloc(abs(f), 1):maxloc(abs(f), 1))
             end if
          end if
          lnode(nnodes) = lab
          xnode(nnodes) = real(np, dp)
       end do
       write(*,'(a,i0,a)') '  band path                : fdf BandLines (', nnodes, ' nodes)'
    end if
    if (nnodes < 2) call auto_path()
    nseg = nnodes - 1
    allocate(nps(nseg))
    lenp = 0.0_dp
    do iseg = 1, nseg
       lenp = lenp + norm2(knode(:,iseg+1) - knode(:,iseg))
    end do
    do iseg = 1, nseg
       if (xnode(iseg+1) >= 1.0_dp) then
          nps(iseg) = nint(xnode(iseg+1))
       else
          nps(iseg) = max(10, nint(240.0_dp*norm2(knode(:,iseg+1)-knode(:,iseg))/max(lenp,1e-12_dp)))
       end if
    end do
    ntot = 1 + sum(nps)
    allocate(kpath(3,ntot), xpath(ntot))
    kpath(:,1) = knode(:,1); xpath(1) = 0.0_dp; npath_pts = 1; xnode(1) = 0.0_dp
    do iseg = 1, nseg
       do l = 1, nps(iseg)
          npath_pts = npath_pts + 1
          kpath(:,npath_pts) = knode(:,iseg) + (knode(:,iseg+1)-knode(:,iseg))*real(l,dp)/nps(iseg)
          dl = norm2(kpath(:,npath_pts) - kpath(:,npath_pts-1))/Bohr2A
          xpath(npath_pts) = xpath(npath_pts-1) + dl
       end do
       xnode(iseg+1) = xpath(npath_pts)
    end do
  end subroutine make_path

  subroutine auto_path()
    real(dp) :: ang, r, B3(3)
    if (dimk == 1) then
       nnodes = 3
       knode(:,1) = -0.5_dp*B1; knode(:,2) = 0.0_dp; knode(:,3) = 0.5_dp*B1
       lnode(1:3) = ['X  ', 'G  ', 'X  ']
    else
       ang = acos(dot_product(B1,B2)/(norm2(B1)*norm2(B2)))*180.0_dp/pi
       r = norm2(B1)/norm2(B2)
       if (abs(r-1.0_dp) < 1e-3_dp .and. (abs(ang-60.0_dp) < 0.5_dp .or. abs(ang-120.0_dp) < 0.5_dp)) then
          nnodes = 4
          knode(:,1) = 0.0_dp
          if (abs(ang-120.0_dp) < 0.5_dp) then
             knode(:,2) = (B1 + 2.0_dp*B2)/3.0_dp
          else
             knode(:,2) = (2.0_dp*B1 + B2)/3.0_dp
          end if
          knode(:,3) = 0.5_dp*B1; knode(:,4) = 0.0_dp
          lnode(1:4) = ['G  ', 'K  ', 'M  ', 'G  ']
       else
          nnodes = 5
          knode(:,1) = 0.0_dp; knode(:,2) = 0.5_dp*B1
          knode(:,3) = 0.5_dp*(B1+B2); knode(:,4) = 0.5_dp*B2; knode(:,5) = 0.0_dp
          lnode(1:5) = ['G  ', 'X  ', 'S  ', 'Y  ', 'G  ']
       end if
       if (dimk == 3) then                      ! add the third direction
          B3 = Bp(:,3)
          nnodes = nnodes + 1
          knode(:,nnodes) = 0.5_dp*B3; lnode(nnodes) = 'Z'
       end if
    end if
    xnode(1:nnodes) = 0.0_dp
    write(*,'(a)') '  band path                : automatic (no BandLines in fdf)'
  end subroutine auto_path

  subroutine compute_path()
    integer :: p, ss
    real(dp), allocatable :: Eo(:,:,:)
    real(dp) :: tt, kk(3)
    tt = wtime()
    allocate(Epath(nbas, nspin, npath_pts), Wpath(ngrp, blo:bhi, nspin, npath_pts))
    if (planar) allocate(Eo(n_odd, nspin, npath_pts))
    !$omp parallel do schedule(dynamic) private(ss)
    do p = 1, npath_pts
       do ss = 1, nspin
          call eigk_w(kpath(:,p), ss, 1, nbas, Epath(:,ss,p), blo, bhi, Wpath(:,:,ss,p))
          if (planar) call eigk(kpath(:,p), ss, map_odd, n_odd, 1, n_odd, Eo(:,ss,p))
       end do
    end do
    !$omp end parallel do
    open(newunit=u, file=trim(outdir)//'/path.dat', status='replace')
    write(u,'(a,i0,a,i0,a,i0)') '# nspin ', nspin, ' nbands ', nbas, ' nodd ', merge(n_odd, 0, planar)
    do p = 1, nnodes
       write(u,'(a,f12.6,1x,a)') '# node ', xnode(p), trim(lnode(p))
    end do
    write(u,'(a)') '# x(1/Ang) ku kv kn  E(all bands, spin1)[spin2]  [E(odd/pi, spin1)[spin2]]  eV absolute'
    do p = 1, npath_pts
       kk = kf(kpath(:,p))
       write(u,'(4f12.6)', advance='no') xpath(p), kk
       write(u,'(*(f13.6))', advance='no') Epath(:,:,p)
       if (planar) write(u,'(*(f13.6))', advance='no') Eo(:,:,p)
       write(u,*)
    end do
    close(u)
    open(newunit=u, file=trim(outdir)//'/path_proj.dat', status='replace')
    write(u,'(a,i0,a,i0,a,i0,a,i0)') '# ngrp ', ngrp, ' band_lo ', blo, ' band_hi ', bhi, ' nspin ', nspin
    write(u,'(a,*(1x,a))') '# groups', (trim(grp_name(p)), p=1,ngrp)
    write(u,'(a)') '# Mulliken weights per k: (group fastest, band, spin)'
    do p = 1, npath_pts
       write(u,'(*(1x,f7.4))') Wpath(:,:,:,p)
    end do
    close(u)
    write(*,'(a,i0,a,f7.1,a)') '  band path                : ', npath_pts, ' k-points  (', wtime()-tt, ' s)'
  end subroutine compute_path

  ! eigenvalues on SIESTA's own Monkhorst-Pack grid (PDOS.kgrid.MonkhorstPack,
  ! else kgrid.MonkhorstPack): lets the plotter rebuild SIESTA's .DOS exactly
  subroutine siesta_kgrid_eigs()
    integer :: i1, i2, i, j, n(3), mm(3,3), ios, u, a, b, c, ss, nk
    real(dp) :: dsp(3), bc(3,3), f(3), k(3), w(3)
    real(dp), allocatable :: e(:)
    character(len=32) :: which
    which = 'PDOS.kgrid.MonkhorstPack'
    call fdf_block(which, i1, i2)
    if (i2 < i1+2) then
       which = 'kgrid.MonkhorstPack'; call fdf_block(which, i1, i2)
    end if
    if (i2 < i1+2) return
    do i = 1, 3
       read(lines(i1+i-1), *, iostat=ios) mm(i,:), w(1)
       if (ios /= 0) then
          read(lines(i1+i-1), *, iostat=ios) mm(i,:); w(1) = 0.0_dp
          if (ios /= 0) return
       end if
       dsp(i) = w(1)
    end do
    do i = 1, 3
       do j = 1, 3
          if (i /= j .and. mm(i,j) /= 0) then
             write(*,'(a)') '  SIESTA k-grid            : non-diagonal grid -> .DOS check by smearing'
             return
          end if
       end do
       n(i) = max(1, mm(i,i))
    end do
    call recip(ucell, bc)
    nk = n(1)*n(2)*n(3)
    allocate(e(nbas))
    open(newunit=u, file=trim(outdir)//'/siesta_kgrid_eig.dat', status='replace')
    write(u,'(a,3i5,a,i0,a,i0,a)') '# grid ', n, '  from '//trim(which)//'  nk ', nk, ' nspin ', nspin, &
         '  E (eV, absolute), one line per k and spin'
    do c = 0, n(3)-1
       do b = 0, n(2)-1
          do a = 0, n(1)-1
             f = ([a, b, c] + dsp)/real(n, dp)
             f = f - nint(f)
             k = f(1)*bc(:,1) + f(2)*bc(:,2) + f(3)*bc(:,3)
             do ss = 1, nspin
                call eigk(k, ss, map_all, nbas, 1, nbas, e)
                write(u,'(i2,*(f12.5))') ss, e
             end do
          end do
       end do
    end do
    close(u)
    write(*,'(a,3i4,a)') '  SIESTA k-grid (.DOS)     : ', n, '  eigenvalues saved for the .DOS check'
  end subroutine siesta_kgrid_eigs

  ! compare with SIESTA's own <label>.bands; points are matched by their
  ! distance along the BandLines path, so the point count need not agree
  subroutine check_against_siesta_bands()
    character(len=1024) :: bf
    integer :: ub, ios, nbb, nss, nkk, ik, iseg, nb2
    real(dp) :: efb, rd(2), d, x, seglen, cum, lenB
    real(dp), allocatable :: eb(:), xs(:), es(:,:), emine(:)
    real(dp) :: k(3)
    bf = join(fdfdir, trim(slabel)//'.bands')
    inquire(file=trim(bf), exist=ex)
    if (.not. ex) then
       call vcheck('SKIP', 'bands vs SIESTA .bands', 0.0_dp, 'no '//trim(slabel)//'.bands next to the fdf'); return
    end if
    open(newunit=ub, file=trim(bf), status='old')
    read(ub,*,iostat=ios) efb
    read(ub,*,iostat=ios) rd
    read(ub,*,iostat=ios) rd
    read(ub,*,iostat=ios) nbb, nss, nkk
    if (ios /= 0 .or. nss /= nspin .or. nbb > nbas .or. nkk < 2) then
       call vcheck('SKIP', 'bands vs SIESTA .bands', 0.0_dp, 'header does not match this HSX (other run?)')
       close(ub); return
    end if
    allocate(eb(nbb*nss), xs(nkk), es(nbb*nss, nkk))
    do ik = 1, nkk
       read(ub,*,iostat=ios) xs(ik), eb
       if (ios /= 0) exit
       es(:,ik) = eb
    end do
    close(ub)
    if (ios /= 0) then
       call vcheck('SKIP', 'bands vs SIESTA .bands', 0.0_dp, 'could not parse the .bands file'); return
    end if
    lenB = 0.0_dp                           ! path length in 1/Bohr, as SIESTA writes it
    do iseg = 1, nnodes-1
       lenB = lenB + norm2(knode(:,iseg+1) - knode(:,iseg))
    end do
    if (abs(xs(nkk) - lenB) > 1e-3_dp*max(lenB,1e-6_dp) + 1e-5_dp) then
       call vcheck('SKIP', 'bands vs SIESTA .bands', xs(nkk)-lenB, &
            'path length differs: BandLines in the fdf changed after the SIESTA run')
       return
    end if
    nb2 = nbb
    allocate(emine(nbas*nspin))
    d = 0.0_dp
    do ik = 1, nkk
       x = xs(ik); cum = 0.0_dp; k = knode(:,nnodes)
       do iseg = 1, nnodes-1
          seglen = norm2(knode(:,iseg+1) - knode(:,iseg))
          if (x <= cum + seglen + 1e-9_dp .and. seglen > 1e-12_dp) then
             k = knode(:,iseg) + (knode(:,iseg+1) - knode(:,iseg))*min(1.0_dp, (x - cum)/seglen)
             exit
          end if
          cum = cum + seglen
       end do
       block
         real(dp) :: e(nbas)
         integer :: ss
         do ss = 1, nspin
            call eigk(k, ss, map_all, nbas, 1, nbas, e)
            ! physically relevant window E_F +- 10 eV (states far above are
            ! sensitive to the conditioning of S(k) and irrelevant here)
            d = max(d, maxval(abs(es((ss-1)*nbb+1:ss*nbb, ik) - e(1:nbb)), &
                              mask=abs(es((ss-1)*nbb+1:ss*nbb, ik) - efb) < 10.0_dp))
         end do
       end block
    end do
    if (d < 2e-3_dp) then
       call vcheck('PASS', 'bands vs SIESTA .bands, E_F+-10 eV: max |dE| (eV)', d, 'HSX reproduces SIESTA (print precision)')
    else if (d < 2e-2_dp) then
       call vcheck('WARN', 'bands vs SIESTA .bands, E_F+-10 eV: max |dE| (eV)', d, 'small deviation: HSX from another SCF step?')
    else
       call vcheck('FAIL', 'bands vs SIESTA .bands, E_F+-10 eV: max |dE| (eV)', d, 'HSX and .bands do not belong to the same run')
    end if
    write(*,'(a,es9.2,a,i0,a)') '  check vs SIESTA .bands   : max |dE| = ', d, ' eV  (', nkk, ' points)'
  end subroutine check_against_siesta_bands

  !==================================================================
  subroutine choose_grids()
    real(dp) :: per_k, avail, x
    integer :: c
    per_k = tdiag*nsp_eff/real(nthreads,dp)
    if (me_on) per_k = per_k*1.15_dp
    if (planar) per_k = per_k*(1.0_dp + (real(n_odd,dp)/nbas)**3 + 0.05_dp)
    avail = max(budget - (wtime() - twall0) - 900*tdiag/nthreads, 0.3_dp*budget)/per_k
    Np = 0; N3 = 0; N1 = 0
    if (dimk == 1) then
       ! 1D: no contours / QPI sweeps to reserve time for; time reversal halves the k-points diagonalised
       avail = max(budget - (wtime() - twall0) - 80*tdiag/nthreads, 0.3_dp*budget)/per_k
       N1 = max(200, min(4000, int(1.6_dp*avail))); N1 = 2*(N1/2)          ! even: DOS check on N/2
    else
       ! periodic plane grid (DOS, E_F, contours, QPI): the accuracy-critical one
       ! equal k spacing along b1 and b2 (anisotropic cells), even sizes
       Np = int(sqrt(2.0_dp*merge(0.25_dp, 0.40_dp, dimk == 3)*avail))
       Np = max(32, min(320, Np))
       x = sqrt(norm2(b1p)/norm2(b2p))
       Na = 2*max(8, min(400, nint(Np*x/2.0_dp))); Nb = 2*max(8, min(400, nint(Np/x/2.0_dp)))
       if (dimk == 3) then
          N3 = 8
          do while (0.5_dp*real(N3+4,dp)**3 <= 0.35_dp*avail .and. N3 < 64)
             N3 = N3 + 4
          end do
       end if
       m_bz  = max(41, min(201, int(sqrt(2.0_dp*0.15_dp*avail))))
       m_ext = max(41, min(201, int(sqrt(2.0_dp*0.15_dp*avail))))
       nkz   = max(31, min(151, int(sqrt(0.10_dp*avail))))
       if (mod(m_ext,2) == 0) m_ext = m_ext + 1
       if (mod(nkz,2) == 0) nkz = nkz + 1
    end if
    write(*,'(a,f7.0,a)') '  time budget              : ', budget, ' s'
    if (dimk == 1) write(*,'(a,i0)') '  1D grid                  : ', N1
    if (dimk >= 2) write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0)') '  grids: plane-periodic ', Na, ' x ', Nb, '   whole-BZ ~', m_bz, &
         '^2   extended ', m_ext, '^2   zoom ', nkz
    if (dimk == 3) write(*,'(a,i0,a)') '  3D grid (DOS, gap)       : ', N3, '^3'
  end subroutine choose_grids

  !==================================================================
  logical function in_bz(kxy)
    real(dp), intent(in) :: kxy(2)
    real(dp) :: g(2)
    integer :: n1, n2
    in_bz = .true.
    do n1 = -2, 2
       do n2 = -2, 2
          if (n1 == 0 .and. n2 == 0) cycle
          g = n1*b1p + n2*b2p
          if (dot_product(kxy,g) > 0.5_dp*dot_product(g,g)*(1+1e-9_dp)) then
             in_bz = .false.; return
          end if
       end do
    end do
  end function in_bz

  ! fold Cartesian-frame k (u,v,n) into the first BZ of the full lattice
  function fold3(k) result(r)
    real(dp), intent(in) :: k(3)
    real(dp) :: r(3), kk(3), best, g(3)
    integer :: n1, n2, n3, m3
    r = k; best = norm2(k); m3 = merge(2, 0, dimk == 3)
    do n1 = -2, 2
       do n2 = merge(-2, 0, dimk >= 2), merge(2, 0, dimk >= 2)
          do n3 = -m3, m3
             g = kf(n1*Bp(:,1) + n2*Bp(:,2) + n3*Bp(:,3))
             kk = k + g
             if (norm2(kk) < best - 1e-9_dp) then; best = norm2(kk); r = kk; end if
          end do
       end do
    end do
  end function fold3

  subroutine ws_polygon()
    real(dp) :: G(2,24), A(2,2), r(2), det, ang(24), tmp(2), x
    integer :: ng, i1, i2, n1, n2, q
    logical :: dup
    ng = 0
    do n1 = -2, 2
       do n2 = -2, 2
          if (n1 == 0 .and. n2 == 0) cycle
          ng = ng + 1; G(:,ng) = n1*b1p + n2*b2p
       end do
    end do
    nv = 0
    do i1 = 1, ng
       do i2 = i1+1, ng
          A(1,:) = G(:,i1); A(2,:) = G(:,i2)
          det = A(1,1)*A(2,2) - A(1,2)*A(2,1)
          if (abs(det) < 1e-10_dp) cycle
          r(1) = ( A(2,2)*0.5_dp*dot_product(G(:,i1),G(:,i1)) - A(1,2)*0.5_dp*dot_product(G(:,i2),G(:,i2)))/det
          r(2) = (-A(2,1)*0.5_dp*dot_product(G(:,i1),G(:,i1)) + A(1,1)*0.5_dp*dot_product(G(:,i2),G(:,i2)))/det
          if (.not. in_bz(r*(1.0_dp - 1e-9_dp))) cycle
          dup = .false.
          do q = 1, nv
             if (norm2(V(:,q) - r) < 1e-7_dp) dup = .true.
          end do
          if (.not. dup .and. nv < 24) then
             nv = nv + 1; V(:,nv) = r; ang(nv) = atan2(r(2), r(1))
          end if
       end do
    end do
    do i1 = 1, nv-1
       do i2 = i1+1, nv
          if (ang(i2) < ang(i1)) then
             x = ang(i1); ang(i1) = ang(i2); ang(i2) = x
             tmp = V(:,i1); V(:,i1) = V(:,i2); V(:,i2) = tmp
          end if
       end do
    end do
    bz_ext = [maxval(abs(V(1,1:nv))), maxval(abs(V(2,1:nv)))]
  end subroutine ws_polygon

  !==================================================================
  ! periodic grid of the b1-b2 plane through Gamma (2D materials: the BZ)
  subroutine periodic_grid_2d()
    integer :: npts, q, a, b, ss, pa, pb, j
    integer, allocatable :: li(:), lj(:)
    real(dp) :: tt, k(3)
    complex(dp), allocatable :: zv(:,:)
    allocate(Eper(nbw, nspin, 0:Na-1, 0:Nb-1), li(Na*Nb), lj(Na*Nb), Wper(ngrp, nbw, nspin, 0:Na-1, 0:Nb-1))
    if (planar) allocate(Eoper(nbo, nspin, 0:Na-1, 0:Nb-1))
    npts = 0
    do b = 0, Nb-1
       do a = 0, Na-1
          pa = mod(Na-a, Na); pb = mod(Nb-b, Nb)
          if (.not. trs .or. b*Na + a <= pb*Na + pa) then
             npts = npts + 1; li(npts) = a; lj(npts) = b
          end if
       end do
    end do
    tt = wtime()
    !$omp parallel do schedule(dynamic) private(a,b,ss,k,zv,j)
    do q = 1, npts
       a = li(q); b = lj(q)
       k = (real(a,dp)/Na)*B1 + (real(b,dp)/Nb)*B2
       do ss = 1, nsp_eff
          if (.not. allocated(zv)) allocate(zv(nbas, nbw))
          call eigk_w(k, ss, blo, bhi, Eper(:,ss,a,b), blo, bhi, Wper(:,:,ss,a,b), zv)
          if (me_on .and. mod(a,zs) == 0 .and. mod(b,zs) == 0) then
             do j = 1, nqb                        ! Euclidean-normalised, R gauge
                Zg(:,j,ss,a/zs,b/zs) = cmplx(zv(:,qlo-blo+j)/sqrt(sum(abs(zv(:,qlo-blo+j))**2)), kind=4)
             end do
          end if
          if (planar) call eigk(k, ss, map_odd, n_odd, olo, ohi, Eoper(:,ss,a,b))
       end do
       !$ if (omp_get_thread_num() == 0) &
       call progress(q, npts, 'plane grid (periodic)    ')
    end do
    !$omp end parallel do
    call progress_done('plane grid (periodic)    ', tt)
    do b = 0, Nb-1
       do a = 0, Na-1
          pa = mod(Na-a, Na); pb = mod(Nb-b, Nb)
          if (trs .and. b*Na + a > pb*Na + pa) then
             Eper(:,1:nsp_eff,a,b) = Eper(:,1:nsp_eff,pa,pb)
             Wper(:,:,1:nsp_eff,a,b) = Wper(:,:,1:nsp_eff,pa,pb)
             if (planar) Eoper(:,1:nsp_eff,a,b) = Eoper(:,1:nsp_eff,pa,pb)
             if (me_on .and. mod(a,zs) == 0 .and. mod(b,zs) == 0) &
                Zg(:,:,:,a/zs,b/zs) = conjg(Zg(:,:,:,pa/zs,pb/zs))     ! time reversal, R gauge
          end if
       end do
    end do
    if (nsp_eff < nspin) then
       Eper(:,2,:,:) = Eper(:,1,:,:); Wper(:,:,2,:,:) = Wper(:,:,1,:,:)
       if (planar) Eoper(:,2,:,:) = Eoper(:,1,:,:)
    end if
    open(newunit=u, file=trim(outdir)//'/bz_grid.dat', status='replace')
    write(u,'(a,i0,a,i0,a,i0,a,i0,a,i0)') '# Na ', Na, ' Nb ', Nb, ' nspin ', nspin, ' band_lo ', blo, ' band_hi ', bhi
    write(u,'(a,4f14.8)') '# B1 B2 in-plane (1/Ang): ', b1p, b2p
    write(u,'(a)') '# i j ku kv  E(band_lo..hi, spin1)[spin2]  (eV - E_F)'
    do b = 0, Nb-1
       do a = 0, Na-1
          write(u,'(2i5,2f11.6)', advance='no') a, b, (real(a,dp)/Na)*b1p + (real(b,dp)/Nb)*b2p
          write(u,'(*(f10.4))') Eper(:,:,a,b) - Ef_eV
       end do
    end do
    close(u)
    if (soc) call write_bz_spin()
  end subroutine periodic_grid_2d

  !==================================================================
  ! spin-orbit: spin character <sz> of every state of bz_grid.dat (same order: band fastest, then
  ! i, then j), as int8 = nint(127 <sz>); used for the spin-resolved JDOS and checked by the audit
  subroutine write_bz_spin()
    integer :: a, b, j, ngh, us
    integer(1), allocatable :: q(:)
    real(dp) :: t
    ngh = ngrp/2
    allocate(q(nbw))
    open(newunit=us, file=trim(outdir)//'/bz_spin.bin', access='stream', form='unformatted', status='replace')
    do b = 0, Nb-1
       do a = 0, Na-1
          do j = 1, nbw
             t = sum(real(Wper(:,j,1,a,b),dp))
             q(j) = int(nint(127.0_dp*(sum(real(Wper(1:ngh,j,1,a,b),dp)) - sum(real(Wper(ngh+1:ngrp,j,1,a,b),dp))) &
                    /max(t, 1.0e-9_dp)), 1)
          end do
          write(us) q
       end do
    end do
    close(us)
  end subroutine write_bz_spin

  subroutine periodic_grid_3d()
    integer :: npts, q, a, b, c, ss, pa, pb, pc, idx, pidx
    integer, allocatable :: la(:), lb(:), lc(:)
    real(dp) :: tt, k(3)
    allocate(E3(nbw, nspin, 0:N3-1, 0:N3-1, 0:N3-1), la(N3**3), lb(N3**3), lc(N3**3), &
             W3(ngrp, nbw, nspin, 0:N3-1, 0:N3-1, 0:N3-1))
    npts = 0
    do c = 0, N3-1
       do b = 0, N3-1
          do a = 0, N3-1
             pa = mod(N3-a, N3); pb = mod(N3-b, N3); pc = mod(N3-c, N3)
             idx = (c*N3 + b)*N3 + a; pidx = (pc*N3 + pb)*N3 + pa
             if (.not. trs .or. idx <= pidx) then
                npts = npts + 1; la(npts) = a; lb(npts) = b; lc(npts) = c
             end if
          end do
       end do
    end do
    tt = wtime()
    !$omp parallel do schedule(dynamic) private(a,b,c,ss,k)
    do q = 1, npts
       a = la(q); b = lb(q); c = lc(q)
       k = (real(a,dp)/N3)*Bp(:,1) + (real(b,dp)/N3)*Bp(:,2) + (real(c,dp)/N3)*Bp(:,3)
       do ss = 1, nsp_eff
          call eigk_w(k, ss, blo, bhi, E3(:,ss,a,b,c), blo, bhi, W3(:,:,ss,a,b,c))
       end do
       !$ if (omp_get_thread_num() == 0) &
       call progress(q, npts, '3D BZ grid (periodic)    ')
    end do
    !$omp end parallel do
    call progress_done('3D BZ grid (periodic)    ', tt)
    do c = 0, N3-1
       do b = 0, N3-1
          do a = 0, N3-1
             pa = mod(N3-a, N3); pb = mod(N3-b, N3); pc = mod(N3-c, N3)
             idx = (c*N3 + b)*N3 + a; pidx = (pc*N3 + pb)*N3 + pa
             if (trs .and. idx > pidx) then
                E3(:,:,a,b,c) = E3(:,:,pa,pb,pc); W3(:,:,:,a,b,c) = W3(:,:,:,pa,pb,pc)
             end if
          end do
       end do
    end do
    if (nsp_eff < nspin) then
       E3(:,2,:,:,:) = E3(:,1,:,:,:); W3(:,:,2,:,:,:) = W3(:,:,1,:,:,:)
    end if
  end subroutine periodic_grid_3d

  subroutine periodic_grid_1d()
    integer :: a, ss
    real(dp) :: tt
    allocate(E1(nbw, nspin, 0:N1-1), W1(ngrp, nbw, nspin, 0:N1-1))
    tt = wtime()
    !$omp parallel do schedule(dynamic) private(ss)
    do a = 0, merge(N1/2, N1-1, trs)
       do ss = 1, nsp_eff
          call eigk_w((real(a,dp)/N1)*B1, ss, blo, bhi, E1(:,ss,a), blo, bhi, W1(:,:,ss,a))
       end do
    end do
    !$omp end parallel do
    if (trs) then
       do a = N1/2+1, N1-1
          E1(:,:,a) = E1(:,:,N1-a); W1(:,:,:,a) = W1(:,:,:,N1-a)
       end do
    end if
    if (nsp_eff < nspin) then
       E1(:,2,:) = E1(:,1,:); W1(:,:,2,:) = W1(:,:,1,:)
    end if
    call write_grid1d()
    call progress_done('1D grid (periodic)       ', tt)
  end subroutine periodic_grid_1d

  !==================================================================
  ! 0D (molecule / cluster): molecular levels from the generalised problem H c = E S c at Gamma
  ! (the HSX of a 0D run has a single cell, so H(k=0) is the full molecular Hamiltonian)
  subroutine molecule_mode()
    integer :: ss, ie, ne, u1, nlev, k, ihomo(2), nfill, kk, s2, i2
    real(dp), allocatable :: e(:,:), occ(:,:), dd(:,:), pd(:,:,:)
    real(4), allocatable :: wg(:,:,:)
    real(dp) :: dev, smin, sg, homo, lumo, de, emin, x, nel0, tt, mom0, hs, ls
    logical, allocatable :: used(:,:)
    tt = wtime()
    write(*,'(/a)') '  ── calculation ────────────────────────────────────────────────────────────'
    if (.not. allocated(map_all)) then
       allocate(map_all(no_u)); map_all = [(k, k = 1, no_u)]
    end if
    nlev = nbas
    allocate(e(nlev, nspin), wg(ngrp, nlev, nspin), occ(nlev, nspin), used(nlev, nspin))
    do ss = 1, nspin
       call eigk_w([0.0_dp, 0.0_dp, 0.0_dp], ss, 1, nlev, e(:,ss), 1, nlev, wg(:,:,ss))
    end do
    write(*,'(a,i0,a,f8.2,a)') '  molecular levels         : ', nlev*nspin, '   (', wtime()-tt, ' s)'
    ! T = 0 aufbau filling of the (integer) electron count: lowest levels first, both spins together
    occ = 0.0_dp; used = .false.; nel0 = real(nint(qtot), dp); nfill = 0
    x = 0.0_dp
    do while (x < nel0 - 1e-6_dp)
       kk = 0; s2 = 0; homo = huge(1.0_dp)
       do ss = 1, nspin
          do k = 1, nlev
             if (.not. used(k,ss) .and. e(k,ss) < homo) then; homo = e(k,ss); kk = k; s2 = ss; end if
          end do
       end do
       if (kk == 0) exit
       used(kk,s2) = .true.
       occ(kk,s2) = min(occf, nel0 - x); x = x + occ(kk,s2)
    end do
    homo = -huge(1.0_dp); lumo = huge(1.0_dp); ihomo = 0
    do ss = 1, nspin
       do k = 1, nlev
          if (occ(k,ss) > 1e-9_dp .and. e(k,ss) > homo) homo = e(k,ss)
          if (occ(k,ss) < occf - 1e-9_dp .and. e(k,ss) < lumo) lumo = e(k,ss)
          if (occ(k,ss) > 1e-9_dp) ihomo(ss) = k
       end do
    end do
    EfOwn = 0.5_dp*(homo + lumo)
    ! checks
    call herm_check([0.0_dp, 0.0_dp, 0.0_dp], 1, dev, smin)
    call vcheck(merge('PASS', 'FAIL', dev < 1e-6_dp), 'Hermiticity max|H-H^+| (eV)', dev, 'HSX file complete and consistent')
    if (smin > 1e-3_dp) then
       call vcheck('PASS', 'smallest eigenvalue of S (> 1e-6: positive definite)', smin, 'overlap well conditioned')
    else
       call vcheck(merge('PASS', 'FAIL', smin > 1e-6_dp), 'smallest eigenvalue of S (> 1e-6: positive definite)', smin, &
            'diffuse orbitals, nearly dependent but positive definite: H c = E S c is solved stably')
    end if
    call vcheck(merge('PASS', 'FAIL', abs(sum(occ) - qtot) < 1e-3_dp), 'electrons in occupied levels - qtot', &
                sum(occ) - qtot, 'aufbau filling at T = 0 reproduces the SIESTA charge')
    call vcheck(merge('PASS', 'INFO', lumo - homo > 1e-3_dp), 'HOMO-LUMO gap (eV)', lumo - homo, &
                'closed shell if > 0; INFO = open shell / degenerate HOMO')
    ! levels
    open(newunit=u1, file=trim(outdir)//'/levels.dat', status='replace')
    write(u1,'(a,i0,a,i0,a,i0)') '# ngrp ', ngrp, ' nspin ', nspin, ' nlev ', nlev
    write(u1,'(a,*(1x,a))') '# groups', (trim(grp_name(k)), k=1,ngrp)
    write(u1,'(a)') '# spin level E-E_F(SIESTA) (eV)  occupation  Mulliken weights (groups)'
    do ss = 1, nspin
       do k = 1, nlev
          write(u1,'(i3,i6,f12.5,f8.3,*(1x,f9.6))') ss, k, e(k,ss) - Ef_eV, occ(k,ss), wg(:,k,ss)
       end do
    end do
    close(u1)
    ! DOS / PDOS of discrete levels: Gaussian of sigma = 50 meV (display convention, stated in the files)
    sg = 0.05_dp; de = 0.002_dp
    emin = min(minval(e) - Ef_eV, -ewin) - 0.5_dp
    ne = int((max(maxval(e) - Ef_eV, ewin) + 0.5_dp - emin)/de) + 1
    ne = min(ne, int(2*15.0_dp/de) + 1)
    emin = max(emin, -15.0_dp)
    allocate(dd(ne, nspin), pd(ne, ngrp, nspin)); dd = 0.0_dp; pd = 0.0_dp
    do ss = 1, nspin
       do k = 1, nlev
          do ie = 1, ne
             x = emin + (ie-1)*de - (e(k,ss) - Ef_eV)
             if (abs(x) > 6*sg) cycle
             dd(ie,ss) = dd(ie,ss) + exp(-0.5_dp*(x/sg)**2)/(sg*sqrt(2*pi))
             pd(ie,:,ss) = pd(ie,:,ss) + real(wg(:,k,ss), dp)*exp(-0.5_dp*(x/sg)**2)/(sg*sqrt(2*pi))   ! not clipped: shells add up to the DOS
          end do
       end do
    end do
    open(newunit=u1, file=trim(outdir)//'/dos.dat', status='replace')
    write(u1,'(a)') '# E-E_F(SIESTA) (eV)  DOS spin1 [spin2]  states/eV per spin channel; molecular levels, Gaussian 0.05 eV'
    do ie = 1, ne
       write(u1,'(f10.4,*(es14.6))') emin + (ie-1)*de, dd(ie,:)
    end do
    close(u1)
    open(newunit=u1, file=trim(outdir)//'/pdos.dat', status='replace')
    write(u1,'(a,i0,a,i0)') '# ngrp ', ngrp, ' nspin ', nspin
    write(u1,'(a,*(1x,a))') '# groups', (trim(grp_name(k)), k=1,ngrp)
    write(u1,'(a)') '# E-E_F(SIESTA)  PDOS(group fastest, spin) states/eV, Mulliken, Gaussian 0.05 eV'
    do ie = 1, ne
       write(u1,'(f10.4,*(es14.6))') emin + (ie-1)*de, ((pd(ie,k,ss), k=1,ngrp), ss=1,nspin)
    end do
    close(u1)
    ! summary
    open(newunit=u, file=trim(outdir)//'/summary.txt', status='replace')
    write(u,'(a,a)')     'system_label : ', trim(slabel)
    write(u,'(a)')       'engine_version : 1.0'
    write(u,'(a,a)')     'fdf_dir : ', trim(fdfdir)
    write(u,'(a,a)')     'fdf_file : ', trim(fdfname)
    write(u,'(a,i0)')    'n_atoms : ', na_u
    write(u,'(a,i0)')    'n_orbitals : ', no_u
    write(u,'(a,i0)')    'nspin : ', nspin
    write(u,'(a,i0)')    'spin_degenerate : 0'
    write(u,'(a,i0)')    'soc : ', merge(nsraw, 0, soc)
    write(u,'(a)')       'time_reversal : 1'
    write(u,'(a)')       'periodic_dims : 0'
    write(u,'(a,f14.6)') 'E_F_eV : ', Ef_eV
    write(u,'(a,f14.6)') 'E_F_own_minus_siesta_eV : ', EfOwn - Ef_eV
    write(u,'(a,f14.6)') 'kT_siesta_eV : ', kT
    write(u,'(a,f14.6)') 'n_electrons : ', qtot
    write(u,'(a,i0)')    'n_levels : ', nlev
    write(u,'(a,2i6)')   'homo_index : ', ihomo
    write(u,'(a,f14.6)') 'HOMO_eV : ', homo - Ef_eV
    write(u,'(a,f14.6)') 'LUMO_eV : ', lumo - Ef_eV
    write(u,'(a,f14.6)') 'overall_gap_eV : ', max(lumo - homo, 0.0_dp)
    write(u,'(a)')       'gap_type : HOMO-LUMO'
    mom0 = 0.0_dp
    if (nspin == 2) mom0 = sum(occ(:,1)) - sum(occ(:,2))
    if (lumo - homo > 1e-3_dp .and. abs(mom0) < 1e-6_dp) then
       write(u,'(a)') 'spin1_type : closed shell (molecule)'
    else if (lumo - homo > 1e-3_dp) then
       write(u,'(a)') 'spin1_type : open shell (molecule, spin-polarised)'
    else
       write(u,'(a)') 'spin1_type : open shell (molecule, partially filled degenerate level)'
    end if
    if (nspin == 2) write(u,'(a,2f12.5)') 'electrons_up_down : ', sum(occ(:,1)), sum(occ(:,2))
    if (nspin == 2) write(u,'(a,f12.5)') 'magnetic_moment_muB : ', sum(occ(:,1)) - sum(occ(:,2))
    close(u)
    write(*,'(/a)') '  ── physics ────────────────────────────────────────────────────────────────'
    write(*,'(a,f10.4,a)') '  HOMO  E-E_F              : ', homo - Ef_eV, ' eV'
    write(*,'(a,f10.4,a)') '  LUMO  E-E_F              : ', lumo - Ef_eV, ' eV'
    write(*,'(a,f10.4,a)') '  HOMO-LUMO gap            : ', lumo - homo, ' eV'
    if (nspin == 2) then
       do ss = 1, 2
          hs = -huge(1.0_dp); ls = huge(1.0_dp)
          do k = 1, nlev
             if (occ(k,ss) > 1e-9_dp) hs = max(hs, e(k,ss))
             if (occ(k,ss) < 1.0_dp - 1e-9_dp) ls = min(ls, e(k,ss))
          end do
          write(*,'(a,a,f10.4,a,f10.4,a,f8.4,a)') '  ', merge('spin up   HOMO / LUMO    :', 'spin down HOMO / LUMO    :', ss == 1), &
               hs - Ef_eV, ' /', ls - Ef_eV, '  (gap ', ls - hs, ' eV)'
       end do
       write(*,'(a,f10.4,a)') '  magnetic moment          : ', sum(occ(:,1)) - sum(occ(:,2)), ' mu_B'
    end if
    if (lumo - homo <= 1e-3_dp) write(*,'(a)') '     (HOMO = LUMO: partially filled degenerate level, open shell)'
    write(*,'(/a)')       '  ── output ─────────────────────────────────────────────────────────────────'
    write(*,'(a,a)')      '  output folder            : ', trim(outdir)
    write(*,'(a,f8.1,a)') '  total compute time       : ', wtime()-twall0, ' s'
  end subroutine molecule_mode

  ! the periodic 1D grid (exact eigenvalues) for the spectral function and the 1D QPI of the plotter
  subroutine write_grid1d()
    integer :: a, u1
    open(newunit=u1, file=trim(outdir)//'/grid1d.dat', status='replace')
    write(u1,'(a,i0,a,i0,a,i0,a,i0,a,f14.8)') '# N1 ', N1, ' nspin ', nspin, ' band_lo ', blo, ' band_hi ', bhi, &
         ' b1_len_1/Ang ', norm2(B1)/Bohr2A
    write(u1,'(a)') '# a  k(1/Ang along the chain)  E(band_lo..hi, spin1)[spin2]  (eV - E_F)'
    do a = 0, N1-1
       write(u1,'(i6,f12.7,*(f10.4))') a, (real(a,dp)/N1)*norm2(B1)/Bohr2A, E1(:,1:nspin,a) - Ef_eV
    end do
    close(u1)
  end subroutine write_grid1d

  !==================================================================
  ! linear simplex integration: DOS on egrid and electrons below e
  ! 1D intervals, 2D triangles (Lehmann-Taut), 3D tetrahedra (Bloechl)
  subroutine simplex_dos(nb, lo, Ebl2, Ebl3, Ebl1, d, ecount, nel, Wb2, Wb3, Wb1, pd)
    integer, intent(in) :: nb, lo
    real(dp), intent(in), optional :: Ebl2(:,:,0:,0:), Ebl3(:,:,0:,0:,0:), Ebl1(:,:,0:)
    real(4), intent(in), optional :: Wb2(:,:,:,0:,0:), Wb3(:,:,:,0:,0:,0:), Wb1(:,:,:,0:)
    real(dp), intent(out), optional :: d(:,:), pd(:,:,:)
    real(dp), intent(in) :: ecount
    real(dp), intent(out) :: nel(2)
    integer :: a, b, c, it, ib, ss, ie, i1, i2, ne, cnr(3,4,6), co2(2,3,2), m, iv4, mb2, ig
    real(dp) :: w, t(4), x, emin, de, g, wv(96), gg
    logical :: dodos, dopd
    dodos = present(d); dopd = present(pd)
    if (dodos) then
       ne = size(d,1); emin = egrid(1); de = egrid(2) - egrid(1); d = 0.0_dp
    end if
    if (dopd) pd = 0.0_dp
    nel = 0.0_dp
    ! 6 tetrahedra of a cube along the 0-7 diagonal
    cnr(:,1,1)=[0,0,0]; cnr(:,2,1)=[1,0,0]; cnr(:,3,1)=[1,1,0]; cnr(:,4,1)=[1,1,1]
    cnr(:,1,2)=[0,0,0]; cnr(:,2,2)=[1,0,0]; cnr(:,3,2)=[1,0,1]; cnr(:,4,2)=[1,1,1]
    cnr(:,1,3)=[0,0,0]; cnr(:,2,3)=[0,1,0]; cnr(:,3,3)=[1,1,0]; cnr(:,4,3)=[1,1,1]
    cnr(:,1,4)=[0,0,0]; cnr(:,2,4)=[0,1,0]; cnr(:,3,4)=[0,1,1]; cnr(:,4,4)=[1,1,1]
    cnr(:,1,5)=[0,0,0]; cnr(:,2,5)=[0,0,1]; cnr(:,3,5)=[1,0,1]; cnr(:,4,5)=[1,1,1]
    cnr(:,1,6)=[0,0,0]; cnr(:,2,6)=[0,0,1]; cnr(:,3,6)=[0,1,1]; cnr(:,4,6)=[1,1,1]
    co2(:,1,1) = [0,0]; co2(:,2,1) = [1,0]; co2(:,3,1) = [0,1]
    co2(:,1,2) = [1,1]; co2(:,2,2) = [1,0]; co2(:,3,2) = [0,1]
    do ss = 1, nspin
       nel(ss) = real(lo-1, dp)
       if (present(Ebl3)) then
          m = size(Ebl3,3); w = 1.0_dp/(6.0_dp*real(m,dp)**3)
          !$omp parallel do private(b,a,it,ib,t,i1,i2,ie,x,g,iv4,wv,ig,gg) reduction(+:nel) collapse(2)
          do c = 0, m-1
             do b = 0, m-1
                do a = 0, m-1
                   do it = 1, 6
                      do ib = 1, nb
                         wv = 0.0_dp
                         do iv4 = 1, 4
                            t(iv4) = Ebl3(ib,ss,mod(a+cnr(1,iv4,it),m),mod(b+cnr(2,iv4,it),m),mod(c+cnr(3,iv4,it),m))
                            if (dopd) wv(1:ngrp) = wv(1:ngrp) + 0.25_dp* &
                                 Wb3(:,ib,ss,mod(a+cnr(1,iv4,it),m),mod(b+cnr(2,iv4,it),m),mod(c+cnr(3,iv4,it),m))
                         end do
                         call sortn(t, 4)
                         nel(ss) = nel(ss) + w*occ_tet(t, ecount)
                         if (dodos) then
                            ! DOS averaged over each energy bin [E - de/2, E + de/2]: exact difference of the
                            ! integrated count, so every state is counted once even for a simplex narrower than de
                            i1 = max(1, floor((t(1)-emin)/de + 0.5_dp) + 1); i2 = min(ne, floor((t(4)-emin)/de + 0.5_dp) + 1)
                            do ie = i1, i2
                               g = w*(occ_tet(t, egrid(ie) + 0.5_dp*de) - occ_tet(t, egrid(ie) - 0.5_dp*de))/de
                               !$omp atomic
                               d(ie,ss) = d(ie,ss) + g
                               if (dopd) then
                                  do ig = 1, ngrp
                                     gg = g*wv(ig)
                                     !$omp atomic
                                     pd(ie,ig,ss) = pd(ie,ig,ss) + gg
                                  end do
                               end if
                            end do
                         end if
                      end do
                   end do
                end do
             end do
          end do
          !$omp end parallel do
       else if (present(Ebl2)) then
          m = size(Ebl2,3); mb2 = size(Ebl2,4); w = 1.0_dp/(2.0_dp*real(m,dp)*real(mb2,dp))
          do b = 0, mb2-1
             do a = 0, m-1
                do it = 1, 2
                   do ib = 1, nb
                      t(1) = Ebl2(ib,ss,mod(a+co2(1,1,it),m),mod(b+co2(2,1,it),mb2))
                      t(2) = Ebl2(ib,ss,mod(a+co2(1,2,it),m),mod(b+co2(2,2,it),mb2))
                      t(3) = Ebl2(ib,ss,mod(a+co2(1,3,it),m),mod(b+co2(2,3,it),mb2))
                      if (dopd) then
                         wv(1:ngrp) = (Wb2(:,ib,ss,mod(a+co2(1,1,it),m),mod(b+co2(2,1,it),mb2)) + &
                                       Wb2(:,ib,ss,mod(a+co2(1,2,it),m),mod(b+co2(2,2,it),mb2)) + &
                                       Wb2(:,ib,ss,mod(a+co2(1,3,it),m),mod(b+co2(2,3,it),mb2)))/3.0_dp
                      end if
                      call sortn(t, 3)
                      nel(ss) = nel(ss) + w*occ_tri(t, ecount)
                      if (dodos) then
                         ! DOS averaged over each energy bin (exact difference of the integrated count)
                         i1 = max(1, floor((t(1)-emin)/de + 0.5_dp) + 1); i2 = min(ne, floor((t(3)-emin)/de + 0.5_dp) + 1)
                         do ie = i1, i2
                            x = egrid(ie)
                            g = w*(occ_tri(t, x + 0.5_dp*de) - occ_tri(t, x - 0.5_dp*de))/de
                            d(ie,ss) = d(ie,ss) + g
                            if (dopd) pd(ie,1:ngrp,ss) = pd(ie,1:ngrp,ss) + g*wv(1:ngrp)
                         end do
                      end if
                   end do
                end do
             end do
          end do
       else
          m = size(Ebl1,3); w = 1.0_dp/real(m,dp)
          do a = 0, m-1
             do ib = 1, nb
                t(1) = Ebl1(ib,ss,a); t(2) = Ebl1(ib,ss,mod(a+1,m)); call sortn(t, 2)
                if (ecount >= t(2)) then
                   nel(ss) = nel(ss) + w
                else if (ecount > t(1)) then
                   nel(ss) = nel(ss) + w*(ecount - t(1))/(t(2) - t(1))
                end if
                if (dodos) then
                   ! a linear band segment fills [t1, t2] uniformly: its states in each energy bin = overlap /
                   ! width (exact; a segment narrower than a bin -- a flat band -- puts its states in its bin)
                   i1 = max(1, floor((t(1)-emin)/de + 0.5_dp) + 1); i2 = min(ne, floor((t(2)-emin)/de + 0.5_dp) + 1)
                   if (dopd) wv(1:ngrp) = 0.5_dp*(Wb1(:,ib,ss,a) + Wb1(:,ib,ss,mod(a+1,m)))
                   do ie = i1, i2
                      x = egrid(ie)
                      if (t(2) - t(1) > 1e-12_dp) then
                         g = w*max(0.0_dp, min(t(2), x + 0.5_dp*de) - max(t(1), x - 0.5_dp*de))/((t(2) - t(1))*de)
                      else
                         g = merge(w/de, 0.0_dp, t(1) >= x - 0.5_dp*de .and. t(1) < x + 0.5_dp*de)
                      end if
                      d(ie,ss) = d(ie,ss) + g
                      if (dopd) pd(ie,1:ngrp,ss) = pd(ie,1:ngrp,ss) + g*wv(1:ngrp)
                   end do
                end if
             end do
          end do
       end if
    end do
  end subroutine simplex_dos

  real(dp) function occ_tri(t, e)
    real(dp), intent(in) :: t(:), e
    if (e <= t(1)) then
       occ_tri = 0.0_dp
    else if (e >= t(3)) then
       occ_tri = 1.0_dp
    else if (e <= t(2)) then
       occ_tri = (e-t(1))**2/max((t(2)-t(1))*(t(3)-t(1)), 1e-30_dp)
    else
       occ_tri = 1.0_dp - (t(3)-e)**2/max((t(3)-t(1))*(t(3)-t(2)), 1e-30_dp)
    end if
  end function occ_tri

  ! linear tetrahedron, Bloechl et al. PRB 49, 16223 (1994), eqs. (A2)-(A4)
  real(dp) function occ_tet(t, e)
    real(dp), intent(in) :: t(:), e
    real(dp) :: e21, e31, e41, e32, e42, e43, x
    e21 = max(t(2)-t(1),1e-30_dp); e31 = max(t(3)-t(1),1e-30_dp); e41 = max(t(4)-t(1),1e-30_dp)
    e32 = max(t(3)-t(2),1e-30_dp); e42 = max(t(4)-t(2),1e-30_dp); e43 = max(t(4)-t(3),1e-30_dp)
    if (e <= t(1)) then
       occ_tet = 0.0_dp
    else if (e < t(2)) then
       x = e - t(1); occ_tet = x**3/(e21*e31*e41)
    else if (e < t(3)) then
       x = e - t(2)
       occ_tet = (e21**2 + 3*e21*x + 3*x**2 - (e31 + e42)*x**3/(e32*e42))/(e31*e41)
    else if (e < t(4)) then
       x = t(4) - e; occ_tet = 1.0_dp - x**3/(e41*e42*e43)
    else
       occ_tet = 1.0_dp
    end if
  end function occ_tet

  real(dp) function dos_tet(t, e)
    real(dp), intent(in) :: t(:), e
    real(dp) :: e21, e31, e41, e32, e42, e43, x
    e21 = max(t(2)-t(1),1e-30_dp); e31 = max(t(3)-t(1),1e-30_dp); e41 = max(t(4)-t(1),1e-30_dp)
    e32 = max(t(3)-t(2),1e-30_dp); e42 = max(t(4)-t(2),1e-30_dp); e43 = max(t(4)-t(3),1e-30_dp)
    if (e <= t(1) .or. e >= t(4)) then
       dos_tet = 0.0_dp
    else if (e < t(2)) then
       x = e - t(1); dos_tet = 3*x**2/(e21*e31*e41)
    else if (e < t(3)) then
       x = e - t(2)
       dos_tet = (3*e21 + 6*x - 3*(e31 + e42)*x**2/(e32*e42))/(e31*e41)
    else
       x = t(4) - e; dos_tet = 3*x**2/(e41*e42*e43)
    end if
  end function dos_tet

  subroutine sortn(t, n)
    integer, intent(in) :: n
    real(dp), intent(inout) :: t(:)
    integer :: a, b
    real(dp) :: x
    do a = 2, n
       x = t(a); b = a - 1
       do while (b >= 1)
          if (t(b) <= x) exit
          t(b+1) = t(b); b = b - 1
       end do
       t(b+1) = x
    end do
  end subroutine sortn

  real(dp) function count_total(e)
    real(dp), intent(in) :: e
    real(dp) :: nn(2)
    if (dimk == 3) then
       call simplex_dos(nbw, blo, Ebl3=E3, ecount=e, nel=nn)
    else if (dimk == 2) then
       call simplex_dos(nbw, blo, Ebl2=Eper, ecount=e, nel=nn)
    else
       call simplex_dos(nbw, blo, Ebl1=E1, ecount=e, nel=nn)
    end if
    count_total = sum(nn(1:nspin))*occf
  end function count_total

  subroutine dos_and_counts()
    integer :: ne, ie, ss, it
    real(dp) :: emin, de, ea, eb2, em, tt, nh(2), qtarget, d1, d2, pord, nq(2)
    real(dp), allocatable :: dh(:,:), dq(:,:)
    logical :: halfdone, qdone
    tt = wtime()
    emin = Ef_eV - ewin; de = 0.002_dp
    ne = int(2*ewin/de) + 1
    allocate(egrid(ne), dos(ne,nspin), pdos(ne,ngrp,nspin))
    do ie = 1, ne; egrid(ie) = emin + (ie-1)*de; end do
    if (dimk == 3) then
       call simplex_dos(nbw, blo, Ebl3=E3, d=dos, ecount=Ef_eV, nel=Nel, Wb3=W3, pd=pdos)
    else if (dimk == 2) then
       call simplex_dos(nbw, blo, Ebl2=Eper, d=dos, ecount=Ef_eV, nel=Nel, Wb2=Wper, pd=pdos)
    else
       call simplex_dos(nbw, blo, Ebl1=E1, d=dos, ecount=Ef_eV, nel=Nel, Wb1=W1, pd=pdos)
    end if
    ! grid convergence: the same DOS from every second grid point (N/2),
    ! both smoothed to 50 meV, compared within E_F +- 1 eV
    ! (an explicit flag: the previous test "dh(1,1) >= 0" failed when rounding gave -1e-18 at a band edge)
    allocate(dh(ne,nspin)); dh = 0.0_dp; halfdone = .false.
    if (dimk == 3 .and. mod(N3,2) == 0) then
       call simplex_dos(nbw, blo, Ebl3=E3(:,:,0::2,0::2,0::2), d=dh, ecount=Ef_eV, nel=nh); halfdone = .true.
    else if (dimk == 2 .and. mod(Na,2) == 0 .and. mod(Nb,2) == 0) then
       call simplex_dos(nbw, blo, Ebl2=Eper(:,:,0::2,0::2), d=dh, ecount=Ef_eV, nel=nh); halfdone = .true.
    else if (dimk == 1 .and. mod(N1,2) == 0) then
       call simplex_dos(nbw, blo, Ebl1=E1(:,:,0::2), d=dh, ecount=Ef_eV, nel=nh); halfdone = .true.
    end if
    if (halfdone) then
       ! error of the DOS on the N grid by Richardson extrapolation: d1 = |DOS(N) - DOS(N/2)|,
       ! d2 = |DOS(N/2) - DOS(N/4)| (both 50 meV, rel. L1 within E_F +- 1 eV, or +- 3 eV in a gap).
       ! Observed order p = log2(d2/d1), limited to [0.5, 2]; error = d1/(2^p - 1). Without an N/4 grid
       ! p = 1 is assumed (error = d1). A slowly converging DOS (p < 1) gets a LARGER estimate, never smaller.
       d1 = conv_metric(dos, dh, egrid, de, ne)
       qdone = .false.
       if (dimk == 3 .and. mod(N3,4) == 0) then
          allocate(dq(ne,nspin)); dq = 0.0_dp
          call simplex_dos(nbw, blo, Ebl3=E3(:,:,0::4,0::4,0::4), d=dq, ecount=Ef_eV, nel=nq); qdone = .true.
       else if (dimk == 2 .and. mod(Na,4) == 0 .and. mod(Nb,4) == 0) then
          allocate(dq(ne,nspin)); dq = 0.0_dp
          call simplex_dos(nbw, blo, Ebl2=Eper(:,:,0::4,0::4), d=dq, ecount=Ef_eV, nel=nq); qdone = .true.
       else if (dimk == 1 .and. mod(N1,4) == 0) then
          allocate(dq(ne,nspin)); dq = 0.0_dp
          call simplex_dos(nbw, blo, Ebl1=E1(:,:,0::4), d=dq, ecount=Ef_eV, nel=nq); qdone = .true.
       end if
       pord = 1.0_dp
       if (qdone) then
          d2 = conv_metric(dh, dq, egrid, de, ne)
          if (d1 > 1e-12_dp .and. d2 > 1e-12_dp) pord = min(2.0_dp, max(0.5_dp, log(d2/d1)/log(2.0_dp)))
       end if
       dos_order = pord; dos_diff = d1
       dos_conv = d1/(2.0_dp**pord - 1.0_dp)
       write(*,'(a,f6.2,a,f5.2,a,f5.2,a)') '  DOS k-grid error         : ', 100*dos_conv, ' %  (N vs N/2 ', 100*d1, &
            ' %, observed order ', pord, ')'
       ef_conv = sum(nh(1:nspin) - Nel(1:nspin))*occf
    else
       write(*,'(a,3i6)') '  DOS convergence          : not computed (grid not divisible by 2): ', Na, Nb, N3
    end if
    ! own Fermi level (T = 0) from the exact simplex count
    ! neutral cells: the electron count is an integer (SIESTA's qtot carries ~1e-5 noise)
    qtarget = qtot
    if (abs(qtot - nint(qtot)) < 1e-3_dp) qtarget = real(nint(qtot), dp)
    ea = Ef_eV - 3.0_dp; eb2 = Ef_eV + 3.0_dp
    do it = 1, 60
       em = 0.5_dp*(ea + eb2)
       if (count_total(em) < qtarget - 1e-7_dp) then; ea = em; else; eb2 = em; end if
       if (eb2 - ea < 1e-5_dp) exit
    end do
    EfOwn = 0.5_dp*(ea + eb2)
    ! upper end of the energy range holding exactly qtarget electrons (gap or
    ! Dirac point between grid points): the T = 0 Fermi level is its middle
    ea = EfOwn; eb2 = EfOwn + 3.0_dp
    do it = 1, 60
       em = 0.5_dp*(ea + eb2)
       if (count_total(em) <= qtarget + 1e-7_dp) then; ea = em; else; eb2 = em; end if
       if (eb2 - ea < 1e-5_dp) exit
    end do
    EfOwn = 0.5_dp*(EfOwn + 0.5_dp*(ea + eb2))
    ! electrons of each spin at the charge-neutral Fermi level of this k-grid: the spin moment is
    ! N_up - N_dn there. SIESTA's E_F comes from its coarser grid with Fermi-Dirac smearing; at that E_F
    ! the dense-grid count is not neutral when a flat band sits at E_F (half-metals).
    if (dimk == 3) then
       call simplex_dos(nbw, blo, Ebl3=E3, ecount=EfOwn, nel=NelOwn)
    else if (dimk == 2) then
       call simplex_dos(nbw, blo, Ebl2=Eper, ecount=EfOwn, nel=NelOwn)
    else
       call simplex_dos(nbw, blo, Ebl1=E1, ecount=EfOwn, nel=NelOwn)
    end if
    ! for an insulator any energy in the gap counts qtot electrons
    do ss = 1, nspin
       dosEF(ss) = dos(nint((Ef_eV-emin)/de)+1, ss)
    end do
    if (nspin == 1) then
       nocc(1) = nint(qtot/occf)
    else
       nocc(1) = nint(Nel(1)); nocc(2) = nint(qtot) - nocc(1)
       if (abs(Nel(2) - nocc(2)) > 0.5_dp) nocc(2) = nint(Nel(2))
    end if
    nocc(nspin+1:2) = nocc(1)
    ! keep a valence and a conduction band inside 1..nbas (and inside the window)
    do ss = 1, 2
       nocc(ss) = max(blo, min(nocc(ss), bhi - 1))
    end do
    ib1 = max(1, minval(nocc(1:nspin)) - 1); ib2 = min(nbas, maxval(nocc(1:nspin)) + 2)
    if (planar) then
       allocate(dos_o(ne,nspin))
       call simplex_dos(nbo, olo, Ebl2=Eoper, d=dos_o, ecount=Ef_eV, nel=Nel_o)
       do ss = 1, nspin
          nocc_o(ss) = nint(Nel_o(ss))
       end do
       nocc_o(nspin+1:2) = nocc_o(1)
       if (minval(nocc_o(1:nspin)) < 1 .or. maxval(nocc_o(1:nspin)) >= n_odd) then
          planar = .false.
       else
          ob1 = max(1, minval(nocc_o(1:nspin)) - 1); ob2 = min(n_odd, maxval(nocc_o(1:nspin)) + 2)
       end if
    end if
    open(newunit=u, file=trim(outdir)//'/dos.dat', status='replace')
    write(u,'(a,i0,a)') '# E-E_F(SIESTA) (eV)  DOS spin1 [spin2]  [pi-DOS ...]   states/eV/cell per spin; ', &
         dimk, 'D linear simplices'
    do ie = 1, ne
       write(u,'(f10.4)', advance='no') egrid(ie) - Ef_eV
       write(u,'(*(es14.6))', advance='no') dos(ie,:)
       if (allocated(dos_o)) write(u,'(*(es14.6))', advance='no') dos_o(ie,:)
       write(u,*)
    end do
    close(u)
    open(newunit=u, file=trim(outdir)//'/pdos.dat', status='replace')
    write(u,'(a,i0,a,i0)') '# ngrp ', ngrp, ' nspin ', nspin
    write(u,'(a,*(1x,a))') '# groups', (trim(grp_name(ie)), ie=1,ngrp)
    write(u,'(a)') '# E-E_F(SIESTA)  PDOS(group fastest, spin) states/eV/cell per spin channel, Mulliken, linear simplices'
    do ie = 1, ne
       write(u,'(f10.4,*(es12.4))') egrid(ie) - Ef_eV, pdos(ie,:,:)
    end do
    close(u)
    write(*,'(a,f8.1,a)') '  DOS / electron count     : ', wtime()-tt, ' s'
  end subroutine dos_and_counts

  ! relative L1 difference of two DOS on egrid, both smoothed with a 50 meV Gaussian, within E_F +- 1 eV
  ! (E_F +- 3 eV if there are almost no states within +- 1 eV: a gap)
  real(dp) function conv_metric(da, db, eg, de, ne)
    integer, intent(in) :: ne
    real(dp), intent(in) :: da(:,:), db(:,:), eg(:), de
    real(dp), allocatable :: gk(:), s1(:), s2(:)
    real(dp) :: num, den, sg, w
    integer :: ng, j, ss, ie, pass
    sg = 0.05_dp; ng = nint(4*sg/de)
    allocate(gk(-ng:ng), s1(ne), s2(ne))
    do j = -ng, ng
       gk(j) = exp(-0.5_dp*(j*de/sg)**2)
    end do
    gk = gk/sum(gk)
    do pass = 1, 2
       w = merge(1.0_dp, 3.0_dp, pass == 1)
       num = 0.0_dp; den = 0.0_dp
       do ss = 1, nsp_eff
          s1 = 0.0_dp; s2 = 0.0_dp
          do ie = 1+ng, ne-ng
             s1(ie) = sum(gk*da(ie-ng:ie+ng,ss)); s2(ie) = sum(gk*db(ie-ng:ie+ng,ss))
          end do
          do ie = 1, ne
             if (abs(eg(ie) - Ef_eV) <= w) then
                num = num + abs(s1(ie) - s2(ie)); den = den + abs(s1(ie))
             end if
          end do
       end do
       if (den*de >= 0.05_dp) exit
    end do
    conv_metric = num/max(den, 1e-12_dp)
  end function conv_metric

  !==================================================================
  ! band nb (sector sec: 1 all, 2 pi) at frame k; nb < 0: gap -nb,-nb+1
  real(dp) function fobj(k, ss, nb, sec)
    real(dp), intent(in) :: k(3)
    integer, intent(in) :: ss, nb, sec
    real(dp) :: e2(2)
    if (nb > 0) then
       if (sec == 1) then
          call eigk(kc(k), ss, map_all, nbas, nb, nb, e2(1:1))
       else
          call eigk(kc(k), ss, map_odd, n_odd, nb, nb, e2(1:1))
       end if
       fobj = e2(1)
    else
       if (sec == 1) then
          call eigk(kc(k), ss, map_all, nbas, -nb, -nb+1, e2)
       else
          call eigk(kc(k), ss, map_odd, n_odd, -nb, -nb+1, e2)
       end if
       fobj = e2(2) - e2(1)
    end if
  end function fobj

  ! compass search over the periodic directions (dimk-dimensional)
  subroutine refine(k, ss, nb, sec, sgn, val)
    real(dp), intent(inout) :: k(3)
    integer, intent(in) :: ss, nb, sec
    real(dp), intent(in) :: sgn
    real(dp), intent(out) :: val
    real(dp) :: step, dirs(3,26), fv(26), tol, bmin
    integer :: it, m, mb, nd, a, b, c
    nd = 0
    do a = -1, 1
       do b = merge(-1, 0, dimk >= 2), merge(1, 0, dimk >= 2)
          do c = merge(-1, 0, dimk == 3), merge(1, 0, dimk == 3)
             if (a == 0 .and. b == 0 .and. c == 0) cycle
             nd = nd + 1; dirs(:,nd) = [real(a,dp), real(b,dp), real(c,dp)]
             dirs(:,nd) = dirs(:,nd)/norm2(dirs(:,nd))
          end do
       end do
    end do
    bmin = norm2(B1)/Bohr2A
    if (dimk >= 2) bmin = min(bmin, norm2(B2)/Bohr2A)
    if (dimk == 3) bmin = min(bmin, norm2(Bp(:,3))/Bohr2A)
    step = 0.5_dp*bmin/real(max(Na, Nb, N3, N1, 8), dp)
    tol = merge(1e-5_dp, 2e-4_dp, nb < 0)
    val = fobj(k, ss, nb, sec)
    do it = 1, 120
       !$omp parallel do
       do m = 1, nd
          fv(m) = fobj(k + step*dirs(:,m), ss, nb, sec)
       end do
       !$omp end parallel do
       mb = maxloc(sgn*fv(1:nd), 1)
       if (sgn*(fv(mb) - val) > 1e-7_dp) then
          k = k + step*dirs(:,mb); val = fv(mb)
       else
          step = 0.5_dp*step
       end if
       if (step < tol) exit
    end do
  end subroutine refine

  ! effective-mass tensor (dimk x dimk), principal values in m_e
  subroutine masses(k, ss, nb, m)
    real(dp), intent(in) :: k(3)
    integer, intent(in) :: ss, nb
    real(dp), intent(out) :: m(3)
    real(dp) :: h, e0, Hs(3,3), ei(3), wk(20), di(3), dj(3)
    integer :: a, b, info
    h = 0.01_dp; m = 0.0_dp; Hs = 0.0_dp
    e0 = fobj(k, ss, nb, 1)
    do a = 1, dimk
       di = 0.0_dp; di(a) = h
       Hs(a,a) = (fobj(k+di, ss, nb, 1) - 2*e0 + fobj(k-di, ss, nb, 1))/h**2
       do b = a+1, dimk
          dj = 0.0_dp; dj(b) = h
          Hs(a,b) = (fobj(k+di+dj, ss, nb, 1) - fobj(k+di-dj, ss, nb, 1) &
                   - fobj(k-di+dj, ss, nb, 1) + fobj(k-di-dj, ss, nb, 1))/(4*h*h)
          Hs(b,a) = Hs(a,b)
       end do
    end do
    call dsyev('N', 'U', dimk, Hs, 3, ei, wk, 20, info)
    do a = 1, dimk
       m(a) = hb2m/sign(max(abs(ei(a)),1e-8_dp), ei(a))
    end do
  end subroutine masses

  real(dp) function cone_v(k, ss, nb, sec)
    real(dp), intent(in) :: k(3)
    integer, intent(in) :: ss, nb, sec
    real(dp) :: h, a, sm
    integer :: m
    h = 0.005_dp; sm = 0.0_dp
    do m = 0, 11
       a = m*pi/6
       sm = sm + fobj(k + h*[cos(a), sin(a), 0.0_dp], ss, -nb, sec)/(2*h)
    end do
    cone_v = sm/12*1e-10_dp/hbar_eVs
  end function cone_v

  subroutine band_edges()
    integer :: ss, iv, ic, l2(2), l3(3), l1(1)
    real(dp) :: k(3), tt, e2(2), f(3)
    tt = wtime()
    if (bhi - blo < 1) then                   ! a single band: nothing to compare
       metal = .true.; dirac = .false.; gap = 0.0_dp; kzc = 0.0_dp
       evbm = Ef_eV; ecbm = Ef_eV; edg = 0.0_dp; kvbm = 0.0_dp; kcbm = 0.0_dp; kdg = 0.0_dp
       vF = 0.0_dp; mvb = 0.0_dp; mcb = 0.0_dp; s_v = 1; s_c = 1; s_d = 1
       return
    end if
    do ss = 1, nsp_eff
       iv = nocc(ss) - blo + 1; ic = iv + 1
       ! VBM
       f = grid_argext(iv, iv, ss, 1.0_dp)
       k = fold3(kf(f(1)*Bp(:,1) + f(2)*Bp(:,2) + f(3)*Bp(:,3)))
       call refine(k, ss, nocc(ss), 1, 1.0_dp, evbm(ss)); kvbm(:,ss) = k
       ! CBM
       f = grid_argext(ic, ic, ss, -1.0_dp)
       k = fold3(kf(f(1)*Bp(:,1) + f(2)*Bp(:,2) + f(3)*Bp(:,3)))
       call refine(k, ss, nocc(ss)+1, 1, -1.0_dp, ecbm(ss)); kcbm(:,ss) = k
       ! smallest direct gap
       f = grid_argext(ic, iv, ss, 0.0_dp)
       k = fold3(kf(f(1)*Bp(:,1) + f(2)*Bp(:,2) + f(3)*Bp(:,3)))
       call refine(k, ss, -nocc(ss), 1, -1.0_dp, edg(ss)); kdg(:,ss) = k
       ! metal if bands overlap, or if E_F cuts through a band (doped)
       metal(ss) = (ecbm(ss) - evbm(ss) < 0.0_dp) .or. &
                   (ecbm(ss) < Ef_eV - 0.02_dp) .or. (evbm(ss) > Ef_eV + 0.02_dp)
       dirac(ss) = (edg(ss) < 0.02_dp .and. dimk == 2 .and. (ecbm(ss) - evbm(ss) >= 0.0_dp .or. &
                    abs(ecbm(ss) - evbm(ss)) < 0.02_dp))
       vF(ss) = 0.0_dp; mvb(:,ss) = 0.0_dp; mcb(:,ss) = 0.0_dp
       if (dirac(ss)) then
          call eigk(kc(kdg(:,ss)), ss, map_all, nbas, nocc(ss), nocc(ss)+1, e2)
          evbm(ss) = e2(1); ecbm(ss) = e2(2); kvbm(:,ss) = kdg(:,ss); kcbm(:,ss) = kdg(:,ss)
          vF(ss) = cone_v(kdg(:,ss), ss, nocc(ss), 1)
       else if (.not. metal(ss)) then
          call masses(kvbm(:,ss), ss, nocc(ss), mvb(:,ss))
          call masses(kcbm(:,ss), ss, nocc(ss)+1, mcb(:,ss))
       end if
       if (planar) then
          iv = nocc_o(ss) - olo + 1
          l2 = minloc(Eoper(iv+1,ss,:,:) - Eoper(iv,ss,:,:)) - 1
          k = fold3(kf((real(l2(1),dp)/Na)*B1 + (real(l2(2),dp)/Nb)*B2))
          call refine(k, ss, -nocc_o(ss), 2, -1.0_dp, edgo(ss)); kdgo(:,ss) = k
          diraco(ss) = edgo(ss) < 0.02_dp
          vFo(ss) = 0.0_dp
          if (diraco(ss)) vFo(ss) = cone_v(kdgo(:,ss), ss, nocc_o(ss), 2)
       end if
    end do
    if (nsp_eff < nspin) then                 ! identical spin channels
       kvbm(:,2) = kvbm(:,1); kcbm(:,2) = kcbm(:,1); evbm(2) = evbm(1); ecbm(2) = ecbm(1)
       kdg(:,2) = kdg(:,1); edg(2) = edg(1); metal(2) = metal(1); dirac(2) = dirac(1)
       vF(2) = vF(1); mvb(:,2) = mvb(:,1); mcb(:,2) = mcb(:,1)
       if (planar) then
          kdgo(:,2) = kdgo(:,1); edgo(2) = edgo(1); diraco(2) = diraco(1); vFo(2) = vFo(1)
       end if
    end if
    s_v = 1; s_c = 1; s_d = 1
    do ss = 1, nspin
       if (evbm(ss) > evbm(s_v)) s_v = ss
       if (ecbm(ss) < ecbm(s_c)) s_c = ss
       if (edg(ss) < edg(s_d)) s_d = ss
    end do
    gap = ecbm(s_c) - evbm(s_v)
    if (any(metal(1:nspin))) gap = min(gap, 0.0_dp)
    kzc = kdg(:,s_d)
    write(*,'(a,f8.1,a)') '  band edges (refined)     : ', wtime()-tt, ' s'
  end subroutine band_edges

  ! fractional coordinates of the grid extremum of band ia (sgn=1 max,
  ! -1 min) or of the smallest difference ia - ib (sgn = 0)
  function grid_argext(ia, ib, ss, sgn) result(f)
    integer, intent(in) :: ia, ib, ss
    real(dp), intent(in) :: sgn
    real(dp) :: f(3)
    integer :: l1(1), l2(2), l3(3)
    f = 0.0_dp
    if (dimk == 3) then
       if (sgn > 0.5_dp) then
          l3 = maxloc(E3(ia,ss,:,:,:))
       else if (sgn < -0.5_dp) then
          l3 = minloc(E3(ia,ss,:,:,:))
       else
          l3 = minloc(E3(ia,ss,:,:,:) - E3(ib,ss,:,:,:))
       end if
       f = real(l3 - 1, dp)/N3
    else if (dimk == 2) then
       if (sgn > 0.5_dp) then
          l2 = maxloc(Eper(ia,ss,:,:))
       else if (sgn < -0.5_dp) then
          l2 = minloc(Eper(ia,ss,:,:))
       else
          l2 = minloc(Eper(ia,ss,:,:) - Eper(ib,ss,:,:))
       end if
       f(1) = real(l2(1) - 1, dp)/Na; f(2) = real(l2(2) - 1, dp)/Nb
    else
       if (sgn > 0.5_dp) then
          l1 = maxloc(E1(ia,ss,:))
       else if (sgn < -0.5_dp) then
          l1 = minloc(E1(ia,ss,:))
       else
          l1 = minloc(E1(ia,ss,:) - E1(ib,ss,:))
       end if
       f(1) = real(l1(1) - 1, dp)/N1
    end if
  end function grid_argext

  ! are the two spin channels identical (non-magnetic spin-polarised run)?
  ! Two collinear spin channels are physically identical when the largest |E_up - E_dn| of all bands on the
  ! whole path is below 2 meV: a non-magnetic solution of a spin-polarised SCF (e.g. a hybrid functional in
  ! HONPAS) keeps residual splittings of ~0.1-1 meV from the finite SCF tolerance; 2 meV is below k_B T at
  ! 25 K and below the accuracy of DFT eigenvalues, while a real exchange splitting is far larger.
  subroutine spin_check()
    real(dp) :: dmax
    spin_degenerate = .false.
    if (nspin /= 2) return
    dmax = maxval(abs(Epath(:,1,:) - Epath(:,2,:)))
    spin_split_max = dmax
    spin_degenerate = dmax < 2e-3_dp
    if (spin_degenerate) then
       write(*,'(a,f7.3,a)') '  spin channels            : identical within ', 1000*dmax, &
            ' meV (non-magnetic) -> one set of bands'
    else
       write(*,'(a,f10.3,a)') '  spin channels            : different (largest splitting ', 1000*dmax, ' meV)'
    end if
  end subroutine spin_check

  !==================================================================
  ! generic evaluation on a list of plane points (u,v) with offset n
  subroutine eval_points(npt, kp, off, partner, Ei, Eo, txt)
    integer, intent(in) :: npt, partner(npt)
    real(dp), intent(in) :: kp(2,npt), off
    real(dp), intent(out) :: Ei(ib1:ib2, nspin, npt), Eo(:,:,:)
    character(len=*), intent(in) :: txt
    integer :: q, ss, nn, j
    integer, allocatable :: lst(:)
    real(dp) :: tt
    integer, allocatable :: pt(:)
    allocate(lst(npt), pt(npt)); nn = 0
    pt = partner
    if (.not. trs) pt = 0                     ! partners are -k images: invalid without time reversal
    do q = 1, npt
       if (pt(q) == 0) then; nn = nn + 1; lst(nn) = q; end if
    end do
    tt = wtime()
    !$omp parallel do schedule(dynamic) private(ss)
    do j = 1, nn
       do ss = 1, nsp_eff
          call eigk(kc([kp(1,lst(j)), kp(2,lst(j)), off]), ss, map_all, nbas, ib1, ib2, Ei(:,ss,lst(j)))
          if (planar) call eigk(kc([kp(1,lst(j)), kp(2,lst(j)), off]), ss, map_odd, n_odd, ob1, ob2, Eo(:,ss,lst(j)))
       end do
       !$ if (omp_get_thread_num() == 0) &
       call progress(j, nn, txt)
    end do
    !$omp end parallel do
    call progress_done(txt, tt)
    do q = 1, npt
       if (pt(q) /= 0) then
          Ei(:,:,q) = Ei(:,:,pt(q))
          if (planar) Eo(:,:,q) = Eo(:,:,pt(q))
       end if
    end do
    if (nsp_eff < nspin) then
       Ei(:,2,:) = Ei(:,1,:)
       if (planar) Eo(:,2,:) = Eo(:,1,:)
    end if
  end subroutine eval_points

  subroutine write_grid(fname, gname, n1, n2, kp, off, Ei, Eo)
    character(len=*), intent(in) :: fname, gname
    integer, intent(in) :: n1, n2
    real(dp), intent(in) :: kp(2,n1*n2), off, Ei(ib1:ib2,nspin,n1*n2), Eo(:,:,:)
    integer :: q, uu
    open(newunit=uu, file=trim(outdir)//'/'//fname, status='replace')
    write(uu,'(a,a,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') '# grid ', gname, ' n1 ', n1, ' n2 ', n2, &
         ' nspin ', nspin, ' band_lo ', ib1, ' band_hi ', ib2, ' nocc1 ', nocc(1), ' nocc2 ', nocc(2)
    write(uu,'(a,i0,a,i0,a,i0,a,i0,a,i0)') '# odd ', merge(1,0,planar), ' odd_lo ', merge(ob1,0,planar), &
         ' odd_hi ', merge(ob2,0,planar), ' nocc_odd1 ', merge(nocc_o(1),0,planar), ' nocc_odd2 ', merge(nocc_o(2),0,planar)
    write(uu,'(a,*(f11.6))') '# vertices ', V(:,1:nv)
    write(uu,'(a,f12.7)') '# plane_offset ', off
    write(uu,'(a)') '# ku kv inBZ  E(band_lo..hi, spin1)[spin2]  [E(odd, ...)]  (1/Ang, eV - E_F); n1 fastest'
    do q = 1, n1*n2
       write(uu,'(2f12.7,i3)', advance='no') kp(:,q), merge(1,0,in_bz(kp(:,q)))
       write(uu,'(*(f11.5))', advance='no') Ei(:,:,q) - Ef_eV
       if (planar) write(uu,'(*(f11.5))', advance='no') Eo(:,:,q) - Ef_eV
       write(uu,*)
    end do
    close(uu)
  end subroutine write_grid

  subroutine grid_bz()
    integer :: ns, nt, npe, is2, it2, ie, q, half, n
    real(dp) :: x, r(2)
    real(dp), allocatable :: kp(:,:), Ei(:,:,:), Eo(:,:,:)
    integer, allocatable :: partner(:)
    npe = 1
    if (nv == 4) then
       x  = norm2(V(:,2)-V(:,1))/norm2(V(:,4)-V(:,1))
       ns = min(4*m_bz, max(21, nint(m_bz*sqrt(x)))); nt = min(4*m_bz, max(21, nint(m_bz/sqrt(x))))
    else
       ns = max(21, m_bz/2); npe = max(8, (2*m_bz)/nv); nt = nv*npe + 1
    end if
    n = ns*nt
    allocate(kp(2,n), partner(n), Ei(ib1:ib2,nspin,n), Eo(max(1,ob2-ob1+1),nspin,n))
    half = (nt-1)/2
    do it2 = 1, nt
       do is2 = 1, ns
          q = is2 + (it2-1)*ns
          if (nv == 4) then
             kp(:,q) = V(:,1) + (V(:,2)-V(:,1))*real(is2-1,dp)/(ns-1) + (V(:,4)-V(:,1))*real(it2-1,dp)/(nt-1)
             partner(q) = 0
             if (q > (ns+1-is2) + (nt-it2)*ns) partner(q) = (ns+1-is2) + (nt-it2)*ns
          else
             ie = min((it2-1)/npe + 1, nv)
             x  = real(it2-1 - (ie-1)*npe, dp)/npe
             r  = V(:,ie) + x*(V(:,mod(ie,nv)+1) - V(:,ie))
             kp(:,q) = r*real(is2-1,dp)/(ns-1)
             partner(q) = 0
             if (is2 == 1 .and. it2 > 1) partner(q) = 1
             if (is2 > 1 .and. it2 > half + 1 .and. mod(nv,2) == 0) partner(q) = is2 + (it2-1-half)*ns
          end if
       end do
    end do
    call eval_points(n, kp, 0.0_dp, partner, Ei, Eo, 'whole BZ (WS polygon)    ')
    call write_grid('bands3d_bz.dat', 'bz', ns, nt, kp, 0.0_dp, Ei, Eo)
  end subroutine grid_bz

  subroutine grid_ext()
    integer :: nx, ny, ix, iy, q, n
    real(dp) :: hx, hy, ratio
    real(dp), allocatable :: kp(:,:), Ei(:,:,:), Eo(:,:,:)
    integer, allocatable :: partner(:)
    hx = 1.5_dp*bz_ext(1); hy = 1.5_dp*bz_ext(2)
    ratio = hx/hy
    nx = max(31, min(4*m_ext, nint(m_ext*sqrt(ratio)))); ny = max(31, min(4*m_ext, nint(m_ext/sqrt(ratio))))
    if (mod(nx,2) == 0) nx = nx + 1
    if (mod(ny,2) == 0) ny = ny + 1
    n = nx*ny
    allocate(kp(2,n), partner(n), Ei(ib1:ib2,nspin,n), Eo(max(1,ob2-ob1+1),nspin,n))
    do iy = 1, ny
       do ix = 1, nx
          q = ix + (iy-1)*nx
          kp(:,q) = [-hx + 2*hx*(ix-1)/real(nx-1,dp), -hy + 2*hy*(iy-1)/real(ny-1,dp)]
          partner(q) = 0
          if (q > (nx+1-ix) + (ny-iy)*nx) partner(q) = (nx+1-ix) + (ny-iy)*nx
       end do
    end do
    call eval_points(n, kp, 0.0_dp, partner, Ei, Eo, 'extended window (VASPKIT)')
    call write_grid('bands3d_ext.dat', 'ext', nx, ny, kp, 0.0_dp, Ei, Eo)
    call wgrd('KX.grd', reshape(kp(1,:), [nx,ny]))
    call wgrd('KY.grd', reshape(kp(2,:), [nx,ny]))
    call wgrd('BAND_HOMO.grd', reshape(Ei(nocc(s_v),s_v,:) - Ef_eV, [nx,ny]))
    call wgrd('BAND_LUMO.grd', reshape(Ei(nocc(s_c)+1,s_c,:) - Ef_eV, [nx,ny]))
  end subroutine grid_ext

  ! zoom around the smallest direct gap; for bulk the plane is shifted
  ! along the normal so that it passes through that k-point
  subroutine grid_zoom()
    integer :: ix, iy, q, n
    real(dp) :: hx, hy
    real(dp), allocatable :: kp(:,:), Ei(:,:,:), Eo(:,:,:)
    integer, allocatable :: partner(:)
    hx = 0.25_dp*bz_ext(1); hy = 0.25_dp*bz_ext(2)
    n = nkz*nkz
    allocate(kp(2,n), partner(n), Ei(ib1:ib2,nspin,n), Eo(max(1,ob2-ob1+1),nspin,n))
    partner = 0
    do iy = 1, nkz
       do ix = 1, nkz
          q = ix + (iy-1)*nkz
          kp(:,q) = kzc(1:2) + [-hx + 2*hx*(ix-1)/real(nkz-1,dp), -hy + 2*hy*(iy-1)/real(nkz-1,dp)]
       end do
    end do
    call eval_points(n, kp, kzc(3), partner, Ei, Eo, 'zoom on smallest gap     ')
    call write_grid('bands3d_zoom.dat', 'zoom', nkz, nkz, kp, kzc(3), Ei, Eo)
  end subroutine grid_zoom

  subroutine wgrd(fname, A)
    character(len=*), intent(in) :: fname
    real(dp), intent(in) :: A(:,:)
    integer :: uu, ii, jj
    open(newunit=uu, file=trim(outdir)//'/'//fname, status='replace')
    do jj = 1, size(A,2)
       write(uu,'(*(es16.8))') (A(ii,jj), ii=1,size(A,1))
    end do
    close(uu)
  end subroutine wgrd

  !==================================================================
  ! eigenvalues il..iu AND eigenvectors (R gauge, z^+ S z = 1) of the full H(k)
  subroutine eigk_v(k, ss, il, iu, e, z)
    real(dp), intent(in) :: k(3)
    integer, intent(in) :: ss, il, iu
    real(dp), intent(out) :: e(iu-il+1)
    complex(dp), intent(out) :: z(nbas, iu-il+1)
    complex(dp), allocatable :: Hk(:,:), Sk(:,:), work(:)
    real(dp), allocatable :: rwork(:), w(:)
    integer, allocatable :: iwork(:), ifail(:)
    integer :: info, m, nb, lw
    nb = nbas; lw = 34*nb
    allocate(Hk(nb,nb), Sk(nb,nb), w(nb), rwork(7*nb), iwork(5*nb), ifail(nb), work(lw))
    call build_hk(k, ss, map_all, nb, Hk, Sk)
    call zhegvx(1, 'V', 'I', 'U', nb, Hk, nb, Sk, nb, 0.0_dp, 0.0_dp, il, iu, &
                0.0_dp, m, w, z, nb, work, lw, rwork, iwork, ifail, info)
    if (info /= 0) then
       write(*,*) 'LAPACK info =', info, '  k (1/Bohr) =', k
       stop ' diagonalisation with eigenvectors failed'
    end if
    e = w(1:iu-il+1)*Ry2eV
  end subroutine eigk_v

  !==================================================================
  ! automatic QPI energies must lie where there are states: an energy inside a gap between two bands
  ! (e.g. CBM + 0.3 eV above a narrow conduction band) gives no contours and an empty QPI map.  Such an
  ! energy (zero simplex DOS) is moved to the nearest energy with states (away from E_F first), 50 meV
  ! into the band.
  subroutine qpi_energies_in_bands()
    integer :: ie, j0, j, k, dir, ne, jj
    real(dp) :: thr, dt(size(egrid))
    ne = size(egrid)
    dt = sum(dos, 2)
    thr = max(1.0e-6_dp, 1.0e-4_dp*maxval(dt))        ! the simplex DOS is exactly zero in a gap
    do ie = 1, nqe
       j0 = nint((qpiE(ie) - egrid(1))/(egrid(2) - egrid(1))) + 1
       if (j0 < 1 .or. j0 > ne) cycle
       if (dt(j0) > thr) cycle
       jj = 0
       do k = 1, ne
          do dir = 1, 2
             j = j0 + merge(1, -1, (dir == 1) .eqv. (qpiE(ie) >= Ef_eV))*k
             if (j >= 1 .and. j <= ne) then
                if (dt(j) > thr) then; jj = j; exit; end if
             end if
          end do
          if (jj > 0) exit
       end do
       if (jj == 0) cycle
       j = jj + nint(sign(0.05_dp, real(jj - j0, dp))/(egrid(2) - egrid(1)))
       if (j >= 1 .and. j <= ne) then
          if (dt(j) > thr) jj = j
       end if
       write(*,'(a,f8.3,a,f8.3,a)') '  QPI energy               : ', qpiE(ie) - Ef_eV, ' eV is in a gap -> ', &
            egrid(jj) - Ef_eV, ' eV'
       qpiE(ie) = egrid(jj)
    end do
  end subroutine qpi_energies_in_bands

  !==================================================================
  ! stride of the stored Bloch vectors: 1 if the full grid fits in 1.2 GB, else 2, 4 or 8 (must divide
  ! the grid); 0 if even every 8th point does not fit
  integer function pick_zs(n1, n2)
    integer, intent(in) :: n1, n2
    real(dp) :: m
    m = real(n1,dp)*n2*nbas*nqb*nsp_eff*8.0_dp
    pick_zs = 1
    do while (m/pick_zs**2 > 1.2e9_dp)
       pick_zs = 2*pick_zs
       if (pick_zs > 8 .or. mod(n1, pick_zs) /= 0 .or. mod(n2, pick_zs) /= 0) then
          pick_zs = 0; return
       end if
    end do
  end function pick_zs

  !==================================================================
  ! bands kept with eigenvectors (for wave-function weighted QPI), atoms
  ! in the plane frame, and the q window of the extended-zone maps
  subroutine setup_qpi_vectors()
    real(dp) :: ew0, mem, hx, hy, mm
    integer :: nv0, ia, ib
    ew0 = ewin; ewin = 2.0_dp
    call band_window(map_all, nbas, qlo, qhi)
    ewin = ew0
    nv0 = nint(qtot/occf)
    qlo = max(blo, min(qlo, nv0)); qhi = min(bhi, max(qhi, nv0+1))
    nqb = qhi - qlo + 1
    mem = real(Na,dp)*Nb*nbas*nqb*nsp_eff*8.0_dp
    allocate(tau_uv(2,na_u), oa1(na_u), oa2(na_u))
    do ia = 1, na_u
       tau_uv(:,ia) = [dot_product(xa(:,ia), eu), dot_product(xa(:,ia), ev)]*Bohr2A
    end do
    oa1 = nbas + 1; oa2 = 0
    do ib = 1, nbas                              ! basis function -> atom (spinors: 2 per orbital)
       ia = orb_atom(merge((ib+1)/2, ib, soc))
       oa1(ia) = min(oa1(ia), ib); oa2(ia) = max(oa2(ia), ib)
    end do
    zs = pick_zs(Na, Nb)
    me_on = zs > 0 .and. nqb >= 1
    if (me_on) then
       allocate(Zg(nbas, nqb, nsp_eff, 0:Na/zs-1, 0:Nb/zs-1))
       if (zs == 1) then
          write(*,'(a,i0,a,i0,a,f7.1,a)') '  Bloch vectors stored     : bands ', qlo, '..', qhi, '  (', mem/1e6_dp, ' MB)'
       else
          write(*,'(a,i0,a,i0,a,f7.1,a,i0,a)') '  Bloch vectors stored     : bands ', qlo, '..', qhi, '  (', &
               mem/zs**2/1e6_dp, ' MB, every ', zs, '. k)'
       end if
    else
       write(*,'(a)') '  Bloch vectors            : too large -> QPI without wave-function weights'
    end if
  end subroutine setup_qpi_vectors

  !==================================================================
  ! Exact constant-energy contours of the linear interpolant on the
  ! triangles of the periodic plane grid, sampled with points of weight
  ! dl/|grad E|.  Each sample remembers the nearest grid vertex (for its
  ! Bloch vector) and its spin and band.
  subroutine contour_samples(e0, maxpts, npt, px, pw, ps, pbd, pva, pvb, nseg, segs, segsb, iunit, ie)
    real(dp), intent(in) :: e0
    integer, intent(in) :: maxpts, iunit, ie
    integer, intent(out) :: npt, nseg
    real(dp), allocatable, intent(out) :: px(:,:), pw(:), segs(:,:)
    integer, allocatable, intent(out) :: ps(:), pbd(:), pva(:), pvb(:), segsb(:,:)
    real(dp) :: t(3), P(2,3), cpts(2,2), g(2), A2(2,2), det, gradn, seglen, ds, lenall, xm(2)
    integer :: ss, ib, a, b, it, c, ncp, j, m, co2(2,3,2), cb, vb(2)
    co2(:,1,1) = [0,0]; co2(:,2,1) = [1,0]; co2(:,3,1) = [0,1]
    co2(:,1,2) = [1,1]; co2(:,2,2) = [1,0]; co2(:,3,2) = [0,1]
    allocate(segs(7, 20000), segsb(4, 20000)); nseg = 0
    do ss = 1, nsp_eff
       do ib = 1, nbw
          do b = 0, Nb-1
             do a = 0, Na-1
                do it = 1, 2
                   do c = 1, 3
                      t(c) = Eper(ib, ss, mod(a+co2(1,c,it),Na), mod(b+co2(2,c,it),Nb))
                      P(:,c) = (real(a+co2(1,c,it),dp)/Na)*b1p + (real(b+co2(2,c,it),dp)/Nb)*b2p
                   end do
                   if (e0 < minval(t) .or. e0 >= maxval(t)) cycle
                   ncp = 0
                   do c = 1, 3
                      j = mod(c,3) + 1
                      if ((t(c) <= e0 .and. t(j) > e0) .or. (t(j) <= e0 .and. t(c) > e0)) then
                         ncp = ncp + 1
                         if (ncp <= 2) cpts(:,ncp) = P(:,c) + (P(:,j)-P(:,c))*(e0-t(c))/(t(j)-t(c))
                      end if
                   end do
                   if (ncp /= 2) cycle
                   A2(1,:) = P(:,2) - P(:,1); A2(2,:) = P(:,3) - P(:,1)
                   det = A2(1,1)*A2(2,2) - A2(1,2)*A2(2,1)
                   if (abs(det) < 1e-14_dp) cycle
                   g(1) = ( A2(2,2)*(t(2)-t(1)) - A2(1,2)*(t(3)-t(1)))/det
                   g(2) = (-A2(2,1)*(t(2)-t(1)) + A2(1,1)*(t(3)-t(1)))/det
                   gradn = max(norm2(g), 1e-6_dp)
                   ! nearest triangle vertex to the segment midpoint
                   xm = 0.5_dp*(cpts(:,1) + cpts(:,2)); cb = 1
                   do c = 2, 3
                      if (norm2(P(:,c) - xm) < norm2(P(:,cb) - xm)) cb = c
                   end do
                   vb = [a + co2(1,cb,it), b + co2(2,cb,it)]
                   if (nseg == size(segs,2)) call grow(segs, segsb)
                   nseg = nseg + 1
                   segs(:,nseg) = [cpts(:,1), cpts(:,2), 1.0_dp/gradn, g]
                   segsb(:,nseg) = [ss, blo+ib-1, vb(1), vb(2)]
                   if (iunit /= 0) write(iunit,'(3i5,4f12.6,es13.5)') ie, ss, blo+ib-1, cpts(:,1), cpts(:,2), 1.0_dp/gradn
                end do
             end do
          end do
       end do
    end do
    lenall = 0.0_dp
    do j = 1, nseg
       lenall = lenall + norm2(segs(3:4,j) - segs(1:2,j))
    end do
    ds = max(lenall/maxpts, 1e-4_dp*min(norm2(b1p), norm2(b2p)))
    npt = 0
    do j = 1, nseg
       npt = npt + max(1, nint(norm2(segs(3:4,j) - segs(1:2,j))/ds))
    end do
    allocate(px(2,npt), pw(npt), ps(npt), pbd(npt), pva(npt), pvb(npt)); npt = 0
    do j = 1, nseg
       seglen = norm2(segs(3:4,j) - segs(1:2,j))
       m = max(1, nint(seglen/ds))
       do c = 1, m
          npt = npt + 1
          px(:,npt) = segs(1:2,j) + (segs(3:4,j) - segs(1:2,j))*(real(c,dp) - 0.5_dp)/m
          pw(npt) = seglen/m*segs(5,j); ps(npt) = segsb(1,j); pbd(npt) = segsb(2,j)
          pva(npt) = segsb(3,j); pvb(npt) = segsb(4,j)
       end do
    end do
  end subroutine contour_samples

  ! spin texture <sigma> on the contour segments (spinor bands): the spinor of
  ! the nearest grid vertex, sx = 2 Re(u*d), sy = 2 Im(u*d), sz = |u|^2 - |d|^2
  ! summed over orbitals (Euclidean-normalised spinor)
  subroutine spin_texture(ie, nseg, segs, segsb)
    integer, intent(in) :: ie, nseg, segsb(:,:)
    real(dp), intent(in) :: segs(:,:)
    integer, save :: ust = 0
    integer :: j, b
    complex(4) :: zv(nbas), cu, cd
    real(dp) :: sv(3)
    if (ie == 1 .or. ust == 0) then
       open(newunit=ust, file=trim(outdir)//'/spin_texture.dat', status='replace')
       write(ust,'(a)') '# ie  u1 v1 u2 v2  <sx> <sy> <sz>   spin expectation on exact contours (spin-orbit)'
    end if
    do j = 1, nseg
       b = segsb(2,j)
       if (b < qlo .or. b > qhi) cycle
       zv = Zg(:, b-qlo+1, 1, modulo(zs*nint(real(segsb(3,j),dp)/zs), Na)/zs, modulo(zs*nint(real(segsb(4,j),dp)/zs), Nb)/zs)
       sv = 0.0_dp
       do b = 1, nbas/2
          cu = zv(2*b-1); cd = zv(2*b)
          sv(1) = sv(1) + 2*real(conjg(cu)*cd); sv(2) = sv(2) + 2*aimag(conjg(cu)*cd)
          sv(3) = sv(3) + abs(cu)**2 - abs(cd)**2
       end do
       write(ust,'(i4,4f11.6,3f9.4)') ie, segs(1:4,j), sv
    end do
    if (ie == nqe) then
       close(ust); ust = 0
    end if
  end subroutine spin_texture

  ! contour JDOS folded into the periodic q-cell, density per (1/Ang)^2
  subroutine jdos_periodic(npt, px, pw, ps, Nq, qm)
    integer, intent(in) :: npt, ps(npt), Nq
    real(dp), intent(in) :: px(2,npt), pw(npt)
    real(dp), intent(out) :: qm(0:Nq-1, 0:Nq-1)
    real(dp) :: q(2), fu, fv, du, dv, w, Binv(2,2), det, cell
    integer :: i2, j, iu0, iv0, ua, va
    det = b1p(1)*b2p(2) - b1p(2)*b2p(1)
    Binv = reshape([b2p(2), -b1p(2), -b2p(1), b1p(1)], [2,2])/det
    cell = abs(det)/real(Nq,dp)**2
    qm = 0.0_dp
    !$omp parallel do schedule(dynamic) private(j,q,fu,fv,iu0,iv0,du,dv,w,ua,va) reduction(+:qm)
    do i2 = 1, npt
       do j = 1, npt
          if (ps(j) /= ps(i2)) cycle                      ! spin-conserving scattering
          q = px(:,j) - px(:,i2)
          fu = Binv(1,1)*q(1) + Binv(1,2)*q(2); fv = Binv(2,1)*q(1) + Binv(2,2)*q(2)
          fu = (fu - floor(fu))*Nq; fv = (fv - floor(fv))*Nq
          iu0 = int(fu); iv0 = int(fv); du = fu - iu0; dv = fv - iv0
          w = pw(i2)*pw(j)
          ua = mod(iu0, Nq); va = mod(iv0, Nq)
          qm(ua, va) = qm(ua, va) + w*(1-du)*(1-dv)
          qm(mod(ua+1,Nq), va) = qm(mod(ua+1,Nq), va) + w*du*(1-dv)
          qm(ua, mod(va+1,Nq)) = qm(ua, mod(va+1,Nq)) + w*(1-du)*dv
          qm(mod(ua+1,Nq), mod(va+1,Nq)) = qm(mod(ua+1,Nq), mod(va+1,Nq)) + w*du*dv
       end do
    end do
    !$omp end parallel do
    qm = qm/cell/(2.0_dp*pi)**4
  end subroutine jdos_periodic

  ! Bloch-overlap weighted JDOS in the extended zone (no folding):
  !   P(q) = sum_s Int dl1 dl2/(|v1||v2|) |M(k1,k2,q)|^2 delta(q - k2 + k1 - G)
  !   M = sum_A exp(-i G.tau_A) sum_{a in A} conj(y_a(k1)) y_a(k2),
  !   y(k) = z(k) exp(-i k.tau)  (cell-periodic part, atomic gauge)
  ! i.e. the q-Fourier component of psi_k1^* psi_k2 for on-site orbital
  ! densities (form factor 1). Suppresses pseudospin/spin-forbidden
  ! backscattering and gives the Bragg-replica intensities.
  subroutine qpi_weighted(npt, px, pw, ps, pbd, pva, pvb, nx, ny, qm, ok)
    integer, intent(in) :: npt, ps(npt), pbd(npt), pva(npt), pvb(npt), nx, ny
    real(dp), intent(in) :: px(2,npt), pw(npt)
    real(dp), intent(out) :: qm(nx, ny)
    logical, intent(out) :: ok
    complex(4), allocatable :: Y(:,:), eG(:,:)
    real(dp), allocatable :: Gl(:,:)
    complex(4) :: QA(na_u), Mq
    real(dp) :: kv(2), q0(2), q(2), fx, fy, dx, dy, w, rmax, ph, mm
    integer :: p, ia, ng, n1, n2, i2, j, ig, ix, iy, nmax
    allocate(Y(nbas, npt))
    ok = .true.
    do p = 1, npt
       if (pbd(p) < qlo .or. pbd(p) > qhi) then
          Y(:,p) = (0.0, 0.0); ok = .false.; cycle
       end if
       ia = zs*nint(real(pva(p),dp)/zs); ig = zs*nint(real(pvb(p),dp)/zs)     ! nearest stored vertex
       kv = (real(ia,dp)/Na)*b1p + (real(ig,dp)/Nb)*b2p
       Y(:,p) = Zg(:, pbd(p)-qlo+1, ps(p), modulo(ia,Na)/zs, modulo(ig,Nb)/zs)
       do ia = 1, na_u
          ph = -dot_product(kv, tau_uv(:,ia))
          Y(oa1(ia):oa2(ia),p) = Y(oa1(ia):oa2(ia),p)*cmplx(cos(ph), sin(ph), kind=4)
       end do
    end do
    ! reciprocal lattice vectors that can bring q0 into the window
    rmax = norm2(Qwin) + norm2(b1p) + norm2(b2p)
    nmax = ceiling(rmax/min(norm2(b1p), norm2(b2p))) + 2
    allocate(Gl(2,(2*nmax+1)**2), eG(na_u,(2*nmax+1)**2)); ng = 0
    do n1 = -nmax, nmax
       do n2 = -nmax, nmax
          q = n1*b1p + n2*b2p
          if (norm2(q) > rmax) cycle
          ng = ng + 1; Gl(:,ng) = q
          do ia = 1, na_u
             ph = -dot_product(q, tau_uv(:,ia)); eG(ia,ng) = cmplx(cos(ph), sin(ph), kind=4)
          end do
       end do
    end do
    dx = 2*Qwin(1)/(nx-1); dy = 2*Qwin(2)/(ny-1)
    qm = 0.0_dp
    !$omp parallel do schedule(dynamic) private(j,ia,QA,q0,ig,q,Mq,mm,w,fx,fy,ix,iy) reduction(+:qm)
    do i2 = 1, npt
       do j = 1, npt
          if (ps(j) /= ps(i2)) cycle
          do ia = 1, na_u
             QA(ia) = dot_product(Y(oa1(ia):oa2(ia),i2), Y(oa1(ia):oa2(ia),j))
          end do
          q0 = px(:,j) - px(:,i2)
          w = pw(i2)*pw(j)
          do ig = 1, ng
             q = q0 + Gl(:,ig)
             if (abs(q(1)) >= Qwin(1) .or. abs(q(2)) >= Qwin(2)) cycle
             Mq = sum(eG(:,ig)*QA)
             mm = real(Mq)**2 + aimag(Mq)**2
             fx = (q(1) + Qwin(1))/dx; fy = (q(2) + Qwin(2))/dy
             ix = int(fx); iy = int(fy); fx = fx - ix; fy = fy - iy
             ix = ix + 1; iy = iy + 1
             if (ix >= nx .or. iy >= ny) cycle
             qm(ix, iy)     = qm(ix, iy)     + w*mm*(1-fx)*(1-fy)
             qm(ix+1, iy)   = qm(ix+1, iy)   + w*mm*fx*(1-fy)
             qm(ix, iy+1)   = qm(ix, iy+1)   + w*mm*(1-fx)*fy
             qm(ix+1, iy+1) = qm(ix+1, iy+1) + w*mm*fx*fy
          end do
       end do
    end do
    !$omp end parallel do
    qm = qm/(dx*dy)/(2.0_dp*pi)**4
  end subroutine qpi_weighted

  !==================================================================
  ! QPI: contours at the main energies (written to cec.dat), contour JDOS
  ! (periodic, qpi.dat), Bloch-weighted QPI (extended zone, qpi_bloch.dat),
  ! and an energy sweep of both for E(q) line cuts (binary files).
  subroutine contours_and_qpi()
    integer, parameter :: Nq = 256, Nqs = 128, nxm = 321, nxs = 161
    real(dp), allocatable :: px(:,:), pw(:), segs(:,:), qm(:,:), qmap(:,:,:), bm(:,:), bmap(:,:,:)
    integer, allocatable :: ps(:), pbd(:), pva(:), pvb(:), segsb(:,:)
    real(dp) :: tt, hx, hy, mm, e0, elo, ehi, de, kr(3), et(1), dEm, dkm, tpair = 0.0_dp
    real(dp), allocatable :: ebuf(:)
    integer :: ie, npt, nseg, ucec, a, b, j, js, nchk, us, mpts_s, mpts_me
    logical :: ok, okall(3)
    real(4), allocatable :: buf(:,:)
    tt = wtime()
    ! energies
    if (user_qpi) then
       nqe = nqpi_user; qpiE(1:nqe) = qpi_user(1:nqe) + Ef_eV
    else if (gap > 0.02_dp) then
       nqe = 3; qpiE = [ecbm(s_c) + 0.10_dp, evbm(s_v) - 0.10_dp, ecbm(s_c) + 0.30_dp]
    else if (any(dirac(1:nspin)) .and. .not. any(metal(1:nspin))) then
       nqe = 3; qpiE = [evbm(1) + 0.50_dp, evbm(1) - 0.50_dp, evbm(1) + 1.00_dp]
    else
       nqe = 3; qpiE = [EfOwn, EfOwn + 0.30_dp, EfOwn - 0.30_dp]
    end if
    if (.not. user_qpi) call qpi_energies_in_bands()
    ! extended-zone q window (same as the figures): 2.05 x BZ, aspect within 0.8..1.25
    hx = 2.05_dp*bz_ext(1); hy = 2.05_dp*bz_ext(2); mm = max(hx, hy)
    Qwin = [max(hx, mm/1.25_dp), max(hy, mm/1.25_dp)]      ! near-square q window
    allocate(qmap(0:Nq-1, 0:Nq-1, nqe), qm(0:Nq-1, 0:Nq-1), bmap(nxm, nxm, nqe), bm(nxm, nxm))
    bmap = 0.0_dp; okall = .true.
    open(newunit=ucec, file=trim(outdir)//'/cec.dat', status='replace')
    write(ucec,'(a,3f12.6)') '# energies (eV - E_F): ', qpiE(1:nqe) - Ef_eV
    write(ucec,'(a)') '# ie spin band  u1 v1 u2 v2  weight=1/|gradE| (A/eV)   exact linear-triangle contours'
    dEm = 0.0_dp; dkm = 0.0_dp; nchk = 0
    do ie = 1, nqe
       e0 = qpiE(ie)
       call contour_samples(e0, 14000, npt, px, pw, ps, pbd, pva, pvb, nseg, segs, segsb, ucec, ie)
       call jdos_periodic(npt, px, pw, ps, Nq, qm); qmap(:,:,ie) = qm
       if (soc .and. me_on) call spin_texture(ie, nseg, segs, segsb)
       ! accuracy of the contour: exact eigenvalue at sampled contour points
       do js = 1, min(nseg, 16)
          j = 1 + (js-1)*max(1, nseg/16)
          kr = kc([0.5_dp*(segs(1,j)+segs(3,j)), 0.5_dp*(segs(2,j)+segs(4,j)), 0.0_dp])
          call eigk(kr, segsb(1,j), map_all, nbas, segsb(2,j), segsb(2,j), et)
          dEm = max(dEm, abs(et(1) - e0)); dkm = max(dkm, abs(et(1) - e0)*segs(5,j)); nchk = nchk + 1
       end do
       deallocate(px, pw, ps, pbd, pva, pvb, segs, segsb)
       if (me_on) then
          tpair = wtime()
          call contour_samples(e0, 2200, npt, px, pw, ps, pbd, pva, pvb, nseg, segs, segsb, 0, ie)
          call qpi_weighted(npt, px, pw, ps, pbd, pva, pvb, nxm, nxm, bm, ok)
          bmap(:,:,ie) = bm; okall(ie) = ok
          tpair = (wtime() - tpair)/max(real(npt,dp)**2, 1.0_dp)
          deallocate(px, pw, ps, pbd, pva, pvb, segs, segsb)
       end if
       write(*,'(a,f8.3,a,i0,a)') '  contour at E-E_F =', e0-Ef_eV, ' eV : ', nseg, ' segments'
    end do
    close(ucec)
    cec_dE = dEm; cec_dk = dkm
    open(newunit=u, file=trim(outdir)//'/qpi.dat', status='replace')
    write(u,'(a,i0,a,i0)') '# Nq ', Nq, ' nE ', nqe
    write(u,'(a,4f14.8)') '# B1 B2 in-plane (1/Ang): ', b1p, b2p
    write(u,'(a,3f12.6)') '# energies (eV - E_F): ', qpiE(1:nqe) - Ef_eV
    write(u,'(a)') '# i j  q = i/Nq b1 + j/Nq b2 : JDOS(E1) [JDOS(E2) JDOS(E3)]  (states^2 / eV^2 / Ang^-2 per cell^2)'
    do b = 0, Nq-1
       do a = 0, Nq-1
          write(u,'(2i5,*(es13.5))') a, b, qmap(a,b,:)
       end do
    end do
    close(u)
    if (me_on) then
       open(newunit=u, file=trim(outdir)//'/qpi_bloch.dat', status='replace')
       write(u,'(a,i0,a,i0,a,i0,a,2f12.6,a,i0)') '# nx ', nxm, ' ny ', nxm, ' nE ', nqe, ' Qx Qy ', Qwin, &
            ' complete ', merge(1, 0, all(okall(1:nqe)))
       write(u,'(a,3f12.6)') '# energies (eV - E_F): ', qpiE(1:nqe) - Ef_eV
       write(u,'(a)') '# Bloch-overlap weighted JDOS, extended zone; q_x = -Qx..Qx (fastest), q_y = -Qy..Qy'
       do b = 1, nxm
          do a = 1, nxm
             write(u,'(*(es12.4))') bmap(a,b,1:nqe)
          end do
       end do
       close(u)
    end if
    write(*,'(a,f8.1,a)') '  contours + QPI (3 E)     : ', wtime()-tt, ' s'

    ! ---------------- energy sweep for E(q) line cuts and the HTML slider
    tt = wtime()
    elo = min(Ef_eV - 1.0_dp, minval(qpiE(1:nqe)) - 0.2_dp)
    ehi = max(Ef_eV + 1.0_dp, maxval(qpiE(1:nqe)) + 0.2_dp)
    elo = max(elo, Ef_eV - 2.0_dp); ehi = min(ehi, Ef_eV + 2.0_dp)
    nsw = min(81, nint((ehi - elo)/0.025_dp) + 1)
    de = (ehi - elo)/(nsw - 1)
    allocate(Esw(nsw))
    do ie = 1, nsw
       Esw(ie) = elo + (ie-1)*de
    end do
    ! points per contour from the measured pair cost (target ~ 0.4 x budget)
    mpts_s = 2500; mpts_me = 700
    if (me_on .and. tpair > 0.0_dp) then
       mpts_me = int(sqrt(max(0.15_dp*budget, 20.0_dp)/(nsw*tpair)))
       mpts_me = max(250, min(1200, mpts_me))
    end if
    open(newunit=ucec, file=trim(outdir)//'/cec_sweep.dat', status='replace')
    write(ucec,'(a,i0,a,2f12.6)') '# nE ', nsw, '  E range (eV - E_F): ', elo - Ef_eV, ehi - Ef_eV
    write(ucec,'(a)') '# ie spin band  u1 v1 u2 v2  weight'
    open(newunit=us, file=trim(outdir)//'/qpi_sweep.bin', access='stream', form='unformatted', status='replace')
    if (me_on) open(newunit=uv2, file=trim(outdir)//'/qpi_bloch_sweep.bin', access='stream', &
                    form='unformatted', status='replace')
    allocate(buf(0:Nqs-1, 0:Nqs-1))
    deallocate(qm); allocate(qm(0:Nqs-1, 0:Nqs-1))
    deallocate(bm); allocate(bm(nxs, nxs))
    do ie = 1, nsw
       call contour_samples(Esw(ie), mpts_s, npt, px, pw, ps, pbd, pva, pvb, nseg, segs, segsb, ucec, ie)
       call jdos_periodic(npt, px, pw, ps, Nqs, qm)
       buf = real(qm, 4); write(us) buf
       deallocate(px, pw, ps, pbd, pva, pvb, segs, segsb)
       if (me_on) then
          call contour_samples(Esw(ie), mpts_me, npt, px, pw, ps, pbd, pva, pvb, nseg, segs, segsb, 0, ie)
          call qpi_weighted(npt, px, pw, ps, pbd, pva, pvb, nxs, nxs, bm, ok)
          write(uv2) real(bm, 4)
          deallocate(px, pw, ps, pbd, pva, pvb, segs, segsb)
       end if
       call progress(ie, nsw, 'QPI energy sweep         ')
    end do
    close(ucec); close(us)
    if (me_on) close(uv2)
    open(newunit=us, file=trim(outdir)//'/qpi_sweep.txt', status='replace')
    write(us,'(a,i0)') 'nE ', nsw
    write(us,'(a,i0)') 'Nq_periodic ', Nqs
    write(us,'(a,i0)') 'n_bloch ', merge(nxs, 0, me_on)
    write(us,'(a,2f12.6)') 'Qwin ', Qwin
    write(us,'(a,4f14.8)') 'B1B2 ', b1p, b2p
    write(us,'(a,*(f10.5))') 'energies ', Esw - Ef_eV
    close(us)
    call progress_done('QPI energy sweep         ', tt)
  end subroutine contours_and_qpi

  subroutine grow(a, b)
    real(dp), allocatable, intent(inout) :: a(:,:)
    integer, allocatable, intent(inout) :: b(:,:)
    real(dp), allocatable :: ta(:,:)
    integer, allocatable :: tb(:,:)
    allocate(ta(size(a,1), 2*size(a,2)), tb(size(b,1), 2*size(b,2)))
    ta(:,1:size(a,2)) = a; tb(:,1:size(b,2)) = b
    call move_alloc(ta, a); call move_alloc(tb, b)
  end subroutine grow

  !==================================================================
  subroutine validation_checks()
    real(dp) :: kr(3), e1(nbas), e2(nbas), dmax, hdev, smin, smn, x, dEF
    integer :: t, ss, p, j
    real(dp), parameter :: rnd(3,4) = reshape([0.1234_dp,0.3711_dp,0.2113_dp, 0.4127_dp,0.0573_dp,0.3301_dp, &
                                                0.2718_dp,0.6931_dp,0.1417_dp, 0.8172_dp,0.4142_dp,0.7071_dp], [3,4])
    hdev = 0.0_dp; smn = huge(1.0_dp)
    do t = 1, 2
       kr = rnd(1,t)*Bp(:,1) + rnd(2,t)*Bp(:,2) + rnd(3,t)*Bp(:,3)
       call herm_check(kr, 1, x, smin); hdev = max(hdev, x); smn = min(smn, smin)
    end do
    call vcheck(merge('PASS','FAIL', hdev < 1e-6_dp), 'Hermiticity max|H(k)-H(k)^+| (eV)', hdev, &
                'HSX file complete and consistent')
    if (smn > 1e-4_dp) then
       call vcheck('PASS', 'smallest eigenvalue of S(k) (> 1e-6: positive definite)', smn, 'overlap well conditioned')
    else
       call vcheck(merge('PASS', 'FAIL', smn > 1e-6_dp), 'smallest eigenvalue of S(k) (> 1e-6: positive definite)', smn, &
            'diffuse orbitals, nearly dependent but positive definite: H c = E S c is solved stably')
    end if
    dmax = 0.0_dp; x = 0.0_dp
    do t = 1, 4
       kr = rnd(1,t)*Bp(:,1) + rnd(2,t)*Bp(:,2) + rnd(3,t)*Bp(:,3)
       do ss = 1, nspin
          call eigk(kr, ss, map_all, nbas, 1, nbas, e1)
          call eigk(kr + Bp(:,1) - 2*Bp(:,2) + Bp(:,3), ss, map_all, nbas, 1, nbas, e2)
          dmax = max(dmax, maxval(abs(e1-e2)))
          call eigk(-kr, ss, map_all, nbas, 1, nbas, e2); x = max(x, maxval(abs(e1-e2)))
       end do
    end do
    call vcheck(merge('PASS','FAIL', dmax < 1e-6_dp), 'periodicity max|E(k+G)-E(k)| (eV)', dmax, &
                'Bloch phases and supercell correct')
    if (trs) then
       call vcheck(merge('PASS','FAIL', x < 1e-5_dp), 'time reversal max|E(-k)-E(k)| (eV)', x, &
                   'justifies computing only half of each grid')
    else
       call vcheck('INFO', 'time reversal max|E(-k)-E(k)| (eV)', x, &
                   'broken by spin-orbit + magnetism (physical); full grids were computed')
    end if
    ! electron count and Fermi level
    x = sum(Nel(1:nspin))*occf - qtot
    dEF = EfOwn - Ef_eV
    if (gap > 0.02_dp) then
       call vcheck(merge('PASS','FAIL', abs(x) < 0.02_dp), 'electrons below E_F - qtot (insulator)', x, &
                   'band filling consistent with the SIESTA charge')
    else
       call vcheck('INFO', 'electrons below SIESTA E_F (T=0) - qtot', x, &
                   'metal: SIESTA uses a smeared occupation, so this need not be 0')
       if (abs(dEF) < max(3.0_dp*kT, 0.03_dp)) then
          call vcheck('PASS', 'own T=0 Fermi level - SIESTA E_F (eV)', dEF, 'agrees within the SIESTA smearing kT')
       else if (abs(dEF) < max(8.0_dp*kT, 0.10_dp)) then
          call vcheck('WARN', 'own T=0 Fermi level - SIESTA E_F (eV)', dEF, &
                      'larger than the smearing: use --fine, or the SCF k-grid was coarse')
       else
          call vcheck('FAIL', 'own T=0 Fermi level - SIESTA E_F (eV)', dEF, 'inconsistent Fermi level')
       end if
    end if
    ! band edges cannot be beaten by any path point (path is a subset of the BZ)
    dmax = 0.0_dp
    do ss = 1, nspin
       if (metal(ss)) cycle
       do p = 1, npath_pts
          dmax = max(dmax, Epath(nocc(ss),ss,p) - evbm(ss), ecbm(ss) - Epath(nocc(ss)+1,ss,p))
       end do
    end do
    call vcheck(merge('PASS','WARN', dmax < 2e-3_dp), 'band edges vs path (eV, <=0 is ideal)', dmax, &
                'refined VBM/CBM are global extrema of the whole BZ')
    if (soc .and. mz_soc_n > 0) then
       call vcheck(merge('PASS', 'FAIL', mz_soc_dg <= max(0.01_dp*abs(mz_soc), 0.005_dp)), &
            'spin moment m_z (all occupied states): k-grid n vs n/2 (muB; limit 1 %)', mz_soc_dg, &
            'Mulliken spin of every occupied spinor state, uniform k-grids')
    end if
    if (dos_conv >= 0.0_dp) then
       if (dos_conv < 0.01_dp) then
          call vcheck('PASS', 'DOS error of the k-grid (Richardson N, N/2, N/4; rel. L1, 50 meV; limit 1 %)', dos_conv, &
                      'k-grid converged for the DOS at experimental resolution')
       else
          call vcheck('FAIL', 'DOS error of the k-grid (Richardson N, N/2, N/4; rel. L1, 50 meV; limit 1 %)', dos_conv, &
                      'above 1 % after automatic refinement: rerun with --fine or --budget=1200')
       end if
    end if
    if (dimk >= 2 .and. cec_dk > 0.0_dp) then
       call vcheck(merge('PASS', 'WARN', cec_dk < 0.01_dp), 'contour position error vs exact H(k) (1/Ang)', cec_dk, &
                   'exact eigenvalue at contour points; dk = |E(k)-E0|/|grad E|; WARN: use --fine')
    end if
    if (planar) then
       do ss = 1, nspin
          if (spin_degenerate .and. ss == 2) cycle
          if (diraco(ss) .and. nspin == 1) call vcheck('INFO', 'pi-pi* gap at the Dirac point (eV)', edgo(ss), &
                                      'mirror-odd sector touches: Dirac cone')
          if (diraco(ss) .and. nspin == 2) call vcheck('INFO', 'pi-pi* gap at the Dirac point, spin '// &
                                      trim(merge('up  ', 'down', ss == 1))//' (eV)', edgo(ss), 'mirror-odd sector touches: Dirac cone')
       end do
    end if
  end subroutine validation_checks

  !==================================================================
  subroutine write_summary()
    integer :: ss
    character(len=64) :: stype
    open(newunit=u, file=trim(outdir)//'/summary.txt', status='replace')
    write(u,'(a,a)')     'system_label : ', trim(slabel)
    write(u,'(a)')       'engine_version : 1.0'
    write(u,'(a,a)')     'fdf_dir : ', trim(fdfdir)
    write(u,'(a,a)')     'fdf_file : ', trim(fdfname)
    write(u,'(a,i0)')    'n_atoms : ', na_u
    write(u,'(a,i0)')    'n_orbitals : ', no_u
    write(u,'(a,i0)')    'nspin : ', nspin
    write(u,'(a,i0)')    'spin_degenerate : ', merge(1, 0, spin_degenerate)
    if (nspin == 2) write(u,'(a,es12.4)') 'spin_split_max_eV : ', spin_split_max
    write(u,'(a,i0)')    'soc : ', merge(nsraw, 0, soc)
    write(u,'(a,i0)')    'time_reversal : ', merge(1, 0, trs)
    write(u,'(a,i0)')    'periodic_dims : ', dimk
    write(u,'(a,f14.6)') 'E_F_eV : ', Ef_eV
    write(u,'(a,f14.6)') 'E_F_own_minus_siesta_eV : ', EfOwn - Ef_eV
    write(u,'(a,f14.6)') 'kT_siesta_eV : ', kT
    write(u,'(a,f14.6)') 'n_electrons : ', qtot
    write(u,'(a,4f14.8)') 'B1B2_inplane : ', b1p, b2p
    write(u,'(a,i0)')    'planar_pi : ', merge(1,0,planar)
    write(u,'(a,3i6)')   'grids_plane_3d_1d : ', max(Na, Nb), N3, N1
    write(u,'(a,2i6)')   'grid_plane_b1_b2 : ', Na, Nb
    write(*,'(/a)') '  ── physics ────────────────────────────────────────────────────────────────'
    if (dimk == 3) write(*,'(a)') '  (band edges and DOS from the full 3D Brillouin zone)'
    if (gap <= 0.02_dp) write(*,'(a,f10.4,a)') '  own Fermi level - SIESTA : ', EfOwn - Ef_eV, ' eV'
    if (nspin == 2) then
       write(*,'(a,2f10.4,a)') '  electrons up / down      : ', NelOwn(1:2), '   (at the charge-neutral E_F)'
       write(*,'(a,f10.4,a)')  '  magnetic moment          : ', NelOwn(1)-NelOwn(2), ' mu_B / cell'
       write(u,'(a,2f12.5)')   'electrons_up_down : ', NelOwn(1:2)
       write(u,'(a,f12.5)')    'magnetic_moment_muB : ', NelOwn(1)-NelOwn(2)
       write(u,'(a,f12.5)')    'magnetic_moment_at_siesta_EF_muB : ', Nel(1)-Nel(2)
    end if
    if (soc .and. mz_soc_n > 0) then
       write(u,'(a,f12.5)')    'magnetic_moment_muB : ', mz_soc
       write(u,'(a,f12.6)')    'soc_moment_grid_change_muB : ', mz_soc_dg
       write(u,'(a,i0)')       'soc_moment_grid : ', mz_soc_n
    end if
    if (soc .and. .not. trs) then
       do ss = 1, 2
          if (.not. sc_ok(ss)) cycle
          write(u,'(a,a,a,4f12.6)') 'soc_', trim(merge('up', 'dn', ss == 1)), '_VBM_E_k : ', sce(1,ss)-Ef_eV, sck(:,1,ss)
          write(u,'(a,a,a,4f12.6)') 'soc_', trim(merge('up', 'dn', ss == 1)), '_CBM_E_k : ', sce(2,ss)-Ef_eV, sck(:,2,ss)
          write(u,'(a,a,a,4f12.6)') 'soc_', trim(merge('up', 'dn', ss == 1)), '_direct_gap_k : ', sce(3,ss), sck(:,3,ss)
          write(u,'(a,a,a,f12.6)')  'soc_', trim(merge('up', 'dn', ss == 1)), '_gap_eV : ', sce(2,ss)-sce(1,ss)
          write(u,'(a,a,a,2f9.4)')  'soc_', trim(merge('up', 'dn', ss == 1)), '_sz_vbm_cbm : ', scsz(:,ss)
       end do
    end if
    do ss = 1, nspin
       if (spin_degenerate .and. ss == 2) cycle
       if (metal(ss) .and. dirac(ss)) then; stype = 'metal: doped Dirac cone (E_F away from the Dirac point)'
       else if (metal(ss)) then; stype = 'metal (bands cross E_F)'
       else if (dirac(ss)) then; stype = 'Dirac semimetal (gap < 20 meV)'
       else; stype = 'semiconductor / insulator'; end if
       write(*,'(a,i0,a,a)')          '  spin ', ss, ' : ', trim(stype)
       if (nocc(ss)+1 > nbas) then
          write(*,'(a,i0,a)')         '     single band (index ', nocc(ss), ') crosses E_F (no gap to report)'
       else
          write(*,'(a,i0,a,i0)')      '     VB / CB band index    : ', nocc(ss), ' / ', nocc(ss)+1
       end if
       if (nocc(ss)+1 <= nbas) write(*,'(a,f10.4,a)') '     VBM  E-E_F            : ', evbm(ss)-Ef_eV, trim(katxt(kvbm(:,ss)))
       if (nocc(ss)+1 <= nbas) write(*,'(a,f10.4,a)') '     CBM  E-E_F            : ', ecbm(ss)-Ef_eV, trim(katxt(kcbm(:,ss)))
       if (.not. metal(ss)) write(*,'(a,f10.4,a)') '     band gap              : ', ecbm(ss)-evbm(ss), ' eV'
       if (nocc(ss)+1 <= nbas) write(*,'(a,f10.4,a)') '     smallest direct gap   : ', edg(ss), trim(katxt(kdg(:,ss)))
       if (dirac(ss)) write(*,'(a,es11.4,a)') '     Dirac velocity v_F    : ', vF(ss), ' m/s'
       if (metal(ss)) write(*,'(a,f10.4,a)')  '     DOS at E_F            : ', dosEF(ss), ' states/eV/cell'
       if (.not. metal(ss) .and. .not. dirac(ss)) then
          write(*,'(a,*(f9.4))') '     m* holes (m_e)        : ', mvb(1:dimk,ss)
          write(*,'(a,*(f9.4))') '     m* electrons (m_e)    : ', mcb(1:dimk,ss)
       end if
       if (planar) then
          write(*,'(a,f10.4,a)') '     pi-pi* smallest gap   : ', edgo(ss), trim(katxt(kdgo(:,ss)))
          if (diraco(ss)) write(*,'(a,es11.4,a)') '     pi Dirac velocity     : ', vFo(ss), ' m/s'
          write(u,'(a,i0,a,4f12.6)') 'spin', ss, '_pi_gap_k : ', edgo(ss), kdgo(:,ss)
          write(u,'(a,i0,a,es14.6)') 'spin', ss, '_pi_velocity_ms : ', vFo(ss)
       end if
       write(u,'(a,i0,a,a)')      'spin', ss, '_type : ', trim(stype)
       write(u,'(a,i0,a,i0)')     'spin', ss, '_nocc : ', nocc(ss)
       write(u,'(a,i0,a,4f12.6)') 'spin', ss, '_VBM_E_k : ', evbm(ss)-Ef_eV, kvbm(:,ss)
       write(u,'(a,i0,a,4f12.6)') 'spin', ss, '_CBM_E_k : ', ecbm(ss)-Ef_eV, kcbm(:,ss)
       write(u,'(a,i0,a,4f12.6)') 'spin', ss, '_direct_gap_k : ', edg(ss), kdg(:,ss)
       ! fundamental gap of this spin channel (CBM - VBM, may be indirect); 0 for a metallic channel,
       ! whose VBM / CBM lines above are then only the extremes of bands nocc and nocc+1
       write(u,'(a,i0,a,f12.6)')  'spin', ss, '_gap_eV : ', &
            merge(max(ecbm(ss) - evbm(ss), 0.0_dp), 0.0_dp, index(stype, 'insulator') > 0)
       write(u,'(a,i0,a,es14.6)') 'spin', ss, '_dirac_velocity_ms : ', vF(ss)
       write(u,'(a,i0,a,3f12.5)') 'spin', ss, '_mass_holes : ', mvb(:,ss)
       write(u,'(a,i0,a,3f12.5)') 'spin', ss, '_mass_electrons : ', mcb(:,ss)
       write(u,'(a,i0,a,f12.6)')  'spin', ss, '_DOS_EF : ', dosEF(ss)
    end do
    write(*,'(a,f10.4,a)') '  overall gap (both spins) : ', max(gap,0.0_dp), ' eV'
    write(u,'(a,f12.6)') 'overall_gap_eV : ', max(gap,0.0_dp)
    if (gap > 0.02_dp .and. s_v /= s_c) then
       write(u,'(a)') 'gap_type : between spin channels'
    else if (gap > 0.02_dp) then
       ! direct if VBM and CBM share k, or if the smallest direct gap exceeds the fundamental gap by
       ! less than 10 meV (far below kT(300 K) and DFT accuracy: "quasi-direct", reported as direct)
       write(u,'(a,f12.6)') 'gap_direct_minus_fundamental_eV : ', edg(s_v) - gap
       write(u,'(a,f12.6)') 'gap_vbm_cbm_k_separation : ', norm2(fold3(kvbm(:,s_v)) - fold3(kcbm(:,s_c)))
       if (norm2(fold3(kvbm(:,s_v)) - fold3(kcbm(:,s_c))) < 1e-3_dp .or. edg(s_v) - gap < 0.010_dp) then
          write(u,'(a)') 'gap_type : direct'; write(*,'(a)') '  gap type                 : direct'
          if (edg(s_v) - gap > 1e-4_dp) write(*,'(a,f7.2,a)') '  (quasi-direct: direct gap only ', &
               1000*(edg(s_v) - gap), ' meV above the fundamental gap)'
       else
          write(u,'(a)') 'gap_type : indirect'; write(*,'(a)') '  gap type                 : indirect'
       end if
    end if
    write(u,'(a,3f12.6)') 'zoom_center : ', kzc
    if (dimk >= 2) write(u,'(a,3f12.6)') 'qpi_energies : ', qpiE(1:nqe) - Ef_eV
    if (user_qpi .and. nqpi_user > 0) write(u,'(a,*(f12.6))') 'qpi_user_energies : ', qpi_user(1:nqpi_user)
    if (dimk >= 2) write(u,'(a,i0)') 'qpi_bloch : ', merge(1, 0, me_on)
    if (dimk >= 2) write(u,'(a,2f12.6)') 'qpi_window : ', Qwin
    write(u,'(a,f12.6)') 'dos_convergence : ', dos_conv
    write(u,'(a,f12.6)') 'dos_diff_N_halfN : ', dos_diff
    write(u,'(a,f12.6)') 'dos_error_order : ', dos_order
    write(u,'(a,f12.6)') 'contour_error_dE_eV : ', cec_dE
    close(u)
    write(*,'(/a)')       '  ── output ─────────────────────────────────────────────────────────────────'
    write(*,'(a,a)')      '  output folder            : ', trim(outdir)
    write(*,'(a,f8.1,a)') '  total compute time       : ', wtime()-twall0, ' s'
  end subroutine write_summary

  subroutine bad_arg(msg)
    character(len=*), intent(in) :: msg
    write(*,'(a)') ' BandQPI-3D: '//msg
    write(*,'(a)') ' usage: siesta3d.x system.fdf [-np N] [--fast|--fine|--budget=S] [--qpi=E1,E2,E3] ...'
    call exit(2)
  end subroutine bad_arg

  subroutine exe_path(link, dir)
    ! folder of the running executable, resolved through symbolic links (Linux /proc; empty elsewhere)
    character(len=*), intent(in) :: link
    character(len=*), intent(out) :: dir
    character(len=1024) :: tmp, line
    integer :: uu, ios
    dir = ''
    write(tmp,'(a,i0,a)') '/tmp/.bandqpi3d_exe_', getpid(), '.txt'
    call execute_command_line('readlink -f '//link//' > '//trim(tmp)//' 2>/dev/null', exitstat=ios)
    open(newunit=uu, file=trim(tmp), status='old', action='read', iostat=ios)
    if (ios /= 0) return
    read(uu,'(a)', iostat=ios) line
    close(uu, status='delete')
    if (ios == 0 .and. len_trim(line) > 0) dir = dirname(trim(line))
  end subroutine exe_path

  subroutine run_plotter()
    character(len=1024) :: d
    integer :: st1, st2
    d = exedir
    inquire(file=trim(d)//'/bandqpi3d.py', exist=ex)
    if (.not. ex) then
       call get_environment_variable('SIESTA3D_HOME', d)
       inquire(file=trim(d)//'/bandqpi3d.py', exist=ex)
    end if
    if (.not. ex) then
       write(*,'(a)') '  bandqpi3d.py not found next to siesta3d.x ('//trim(exedir)//') nor in $SIESTA3D_HOME:'
       write(*,'(a)') '  keep both files in one folder, then run  python3 /that/folder/bandqpi3d.py '//trim(outdir)
       call exit(1)
    end if
    write(*,'(/a)') '  ── figures, interactive page and checks ───────────────────────────────────'
    write(cmd,'(a,i0,a,i0,a)') 'OMP_NUM_THREADS=', nthreads, ' OPENBLAS_NUM_THREADS=', nthreads, ' '
    cmd = trim(cmd)//' python3 "'//trim(d)//'/bandqpi3d.py" "'//trim(outdir)//'"'
    if (.not. show) cmd = trim(cmd)//' --noshow'
    cmd = trim(cmd)//' '//trim(plotargs)
    st1 = 0; st2 = 0
    call execute_command_line(trim(cmd), exitstat=st1)
    ! independent audit (re-reads the raw files; does not use the plotting code)
    call execute_command_line('python3 "'//trim(d)//'/bandqpi3d.py" --audit "'//trim(outdir)//'"', exitstat=st2)
    write(*,'(a)') '  ───────────────────────────────────────────────────────────────────────────'
    if (st1 /= 0 .or. st2 /= 0) then
       write(*,'(a)') '  BandQPI-3D v1.0  ·  finished: an error or a FAILED check (see above)'
       call exit(1)
    end if
    write(*,'(a)') '  BandQPI-3D v1.0  ·  Satish Kumar  ·  done'
    write(*,'(a)') ''
  end subroutine run_plotter
end program siesta3d

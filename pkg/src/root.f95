
!*****************************************************************************************
!>
!  Root solver methods for:
!
!  * Bracked interval
!  * Without derivatives
!
!### Author
!  * Jacob Williams
!
!@note The default real kind (`wp`) can be
!      changed using optional preprocessor flags.
!      This library was built with real kind:

! Copyright (c) 2021-2022, Jacob Williams All rights reserved.

! Roots-Fortran: Root solvers for modern Fortran https://github.com/jacobwilliams/roots-fortran
! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
! * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
! * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
! * The names of its contributors may not be used to endorse or promote products derived from this software without specific prior written permission.
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
!FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
!STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


    module root_module

    use iso_fortran_env

    implicit none

    private


    integer,parameter,public :: root_module_rk = real64   !! real kind used by this module [8 bytes]

    integer,parameter :: wp = root_module_rk  !! local copy of `root_module_rk` with a shorter name
    integer,parameter :: name_len = 32  !! max length of the method names
    real(wp),parameter,private :: golden_ratio = (1.0_wp + sqrt(5.0_wp))/2.0_wp !! golden ratio `phi`

    type,public :: root_method
        !! a type to define enums for the different methods
        private
        integer :: id = 0 !! unique ID code for the method
        character(len=name_len) :: name = '' !! name of the method
    end type root_method

    type(root_method),parameter,public :: root_method_null                 = root_method(0,  '')                     !! enum type for an invalid method
    type(root_method),parameter,public :: root_method_brent                = root_method(1,  'brent')                !! enum type for `brent` method
    type(root_method),parameter,public :: root_method_bisection            = root_method(2,  'bisection')            !! enum type for `bisection` method
    type(root_method),parameter,public :: root_method_regula_falsi         = root_method(3,  'regula_falsi')         !! enum type for `regula_falsi` method
    type(root_method),parameter,public :: root_method_illinois             = root_method(4,  'illinois')             !! enum type for `illinois` method
    type(root_method),parameter,public :: root_method_anderson_bjorck      = root_method(5,  'anderson_bjorck')      !! enum type for `anderson_bjorck` method
    type(root_method),parameter,public :: root_method_ridders              = root_method(6,  'ridders')              !! enum type for `ridders` method
    type(root_method),parameter,public :: root_method_pegasus              = root_method(7,  'pegasus')              !! enum type for `pegasus` method
    type(root_method),parameter,public :: root_method_bdqrf                = root_method(8,  'bdqrf')                !! enum type for `bdqrf` method
    type(root_method),parameter,public :: root_method_muller               = root_method(9,  'muller')               !! enum type for `muller` method
    type(root_method),parameter,public :: root_method_brenth               = root_method(10, 'brenth')               !! enum type for `brenth` method
    type(root_method),parameter,public :: root_method_brentq               = root_method(11, 'brentq')               !! enum type for `brentq` method
    type(root_method),parameter,public :: root_method_chandrupatla         = root_method(12, 'chandrupatla')         !! enum type for `chandrupatla` method
    type(root_method),parameter,public :: root_method_toms748              = root_method(13, 'toms748')              !! enum type for `toms748` method
    type(root_method),parameter,public :: root_method_zhang                = root_method(14, 'zhang')                !! enum type for `zhang` method
    type(root_method),parameter,public :: root_method_anderson_bjorck_king = root_method(15, 'anderson_bjorck_king') !! enum type for `anderson_bjorck_king` method
    type(root_method),parameter,public :: root_method_blendtf              = root_method(16, 'blendtf')              !! enum type for `blendtf` method
    type(root_method),parameter,public :: root_method_barycentric          = root_method(17, 'barycentric')          !! enum type for `barycentric` method
    type(root_method),parameter,public :: root_method_itp                  = root_method(18, 'itp')                  !! enum type for `itp` method

    type(root_method),parameter,dimension(*),public :: set_of_root_methods = &
        [ root_method_brent,                &
          root_method_bisection,            &
          root_method_regula_falsi,         &
          root_method_illinois,             &
          root_method_anderson_bjorck,      &
          root_method_ridders,              &
          root_method_pegasus,              &
          root_method_bdqrf,                &
          root_method_muller,               &
          root_method_brenth,               &
          root_method_brentq,               &
          root_method_chandrupatla,         &
          root_method_toms748,              &
          root_method_zhang,                &
          root_method_anderson_bjorck_king, &
          root_method_blendtf,              &
          root_method_barycentric,          &
          root_method_itp                   ]  !! list of the available methods (see [[root_scalar]])

    type,abstract,public :: root_solver
        !! abstract class for the root solver methods
        private
        procedure(func),pointer :: f => null()  !! user function to find the root of
        real(wp) :: ftol = 0.0_wp       !! absolute tolerance for `f=0`
        real(wp) :: rtol = 1.0e-2_wp    !! relative tol for x
        real(wp) :: atol = 1.0e-4_wp   !! absolute tol for x [not used by all methods]
        integer  :: maxiter = 2000      !! maximum number of iterations [not used by all methods]
    contains
        private
        procedure,public :: initialize => initialize_root_solver !! initialize the class [must be called first]
        procedure,public :: solve !! main routine for finding the root
        procedure(root_f),deferred :: find_root !! root solver function. Meant to be
                                                !! called from [[solve]], which handles some common
                                                !! startup tasks. All these routines assume that
                                                !! \( f(a_x) \) and \( f(b_x) \) have opposite signs.
        procedure :: get_fa_fb
        procedure :: converged
        procedure :: solution !! to check `f` value against `ftol`
    end type root_solver

    type,extends(root_solver),public :: brent_solver
    !! Classic brent (zeroin) root solver
    private
    contains
    private
    procedure,public :: find_root => brent
    end type brent_solver

    type,extends(root_solver),public :: bisection_solver
    !! Classic bisection root solver
    private
    contains
    private
    procedure,public :: find_root => bisection
    end type bisection_solver

    type,extends(root_solver),public :: regula_falsi_solver
    !! Classic Regula Falsi root solver
    private
    contains
    private
    procedure,public :: find_root => regula_falsi
    end type regula_falsi_solver

    type,extends(root_solver),public :: illinois_solver
    !! Illinois (modified Regula Falsi) root solver
    private
    contains
    private
    procedure,public :: find_root => illinois
    end type illinois_solver

    type,extends(root_solver),public :: anderson_bjorck_solver
    !! anderson bjorck root solver
    private
    contains
    private
    procedure,public :: find_root => anderson_bjorck
    end type anderson_bjorck_solver

    type,extends(root_solver),public :: ridders_solver
    !! ridders root solver
    private
    contains
    private
    procedure,public :: find_root => ridders
    end type ridders_solver

    type,extends(root_solver),public :: pegasus_solver
    !! pegasus root solver
    private
    contains
    private
    procedure,public :: find_root => pegasus
    end type pegasus_solver

    type,extends(root_solver),public :: bdqrf_solver
    !! bdqrf root solver
    private
    contains
    private
    procedure,public :: find_root => bdqrf
    end type bdqrf_solver

    type,extends(root_solver),public :: muller_solver
    !! muller root solver
    private
    contains
    private
    procedure,public :: find_root => muller
    end type muller_solver

    type,extends(root_solver),public :: brenth_solver
    !! brenth root solver
    private
    contains
    private
    procedure,public :: find_root => brenth
    end type brenth_solver

    type,extends(root_solver),public :: brentq_solver
    !! brentq root solver
    private
    contains
    private
    procedure,public :: find_root => brentq
    end type brentq_solver

    type,extends(root_solver),public :: chandrupatla_solver
    !! chandrupatla root solver
    private
    contains
    private
    procedure,public :: find_root => chandrupatla
    end type chandrupatla_solver

    type,extends(root_solver),public :: toms748_solver
    !! TOMS748 root solver
    private
    contains
    private
    procedure,public :: find_root => toms748
    end type toms748_solver

    type,extends(root_solver),public :: zhang_solver
    !! zhang root solver
    private
    contains
    private
    procedure,public :: find_root => zhang
    end type zhang_solver

    type,extends(root_solver),public :: anderson_bjorck_king_solver
    !! anderson-bjorck-king root solver
    private
    contains
    private
    procedure,public :: find_root => anderson_bjorck_king
    end type anderson_bjorck_king_solver

    type,extends(root_solver),public :: blendtf_solver
    !! blendtf root solver
    private
    contains
    private
    procedure,public :: find_root => blendtf
    end type blendtf_solver

    type,extends(root_solver),public :: barycentric_solver
    !! barycentric root solver
    private
    contains
    private
    procedure,public :: find_root => barycentric
    end type barycentric_solver

    type,extends(root_solver),public :: itp_solver
    !! ITP root solver
    private

    ! tuning parameters for ITP:
    real(wp) :: k1 = 0.1_wp  !! from (0, inf)
    real(wp) :: k2 = 0.98_wp * (1.0_wp + golden_ratio)  !! from [1, 1+phi]
    integer  :: n0 = 1

    contains
    private
    procedure,public :: find_root => itp
    procedure,public :: set_optional_inputs => itp_optional_inputs
    end type itp_solver

    abstract interface
        function func(me,x) result(f)
            !! Interface to the function to be minimized
            !! (Object-oriented version).
            !! It should evaluate f(x) for any x in the interval (ax,bx)
            import :: root_solver, wp
            implicit none
            class(root_solver),intent(inout) :: me
            real(wp),intent(in) :: x !! independant variable
            real(wp) :: f !! f(x)
        end function func
        function func2(x) result(f)
            !! Interface to the function to be minimized
            !! (Functional version).
            !! It should evaluate f(x) for any x in the interval (ax,bx)
            import :: wp
            implicit none
            real(wp),intent(in) :: x !! independant variable
            real(wp) :: f !! f(x)
        end function func2
        subroutine root_f(me,ax,bx,fax,fbx,xzero,fzero,iflag)
            !! Root solver function interface
            import :: root_solver, wp
            implicit none
            class(root_solver),intent(inout) :: me
            real(wp),intent(in) :: ax !! left endpoint of initial interval
            real(wp),intent(in) :: bx !! right endpoint of initial interval
            real(wp),intent(in) :: fax !! `f(ax)`
            real(wp),intent(in) :: fbx !! `f(bx)`
            real(wp),intent(out) :: xzero !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
            real(wp),intent(out) :: fzero !! value of `f` at the root (`f(xzero)`)
            integer,intent(out) :: iflag !! status flag
        end subroutine root_f
    end interface

    interface root_scalar
        module procedure :: root_scalar_by_name, &
                            root_scalar_by_type
    end interface root_scalar
    public :: root_scalar

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Initialize the [[root_solver]] class.
!
!  Note that all optional inputs are not used for all methods.

    subroutine initialize_root_solver(me,f,ftol,rtol,atol,maxiter)

    implicit none

    class(root_solver),intent(out) :: me
    procedure(func)               :: f        !! user function `f(x)` to find the root of
    real(wp),intent(in),optional  :: ftol     !! absolute tolerance for `f=0`
    real(wp),intent(in),optional  :: rtol     !! relative tol for x
    real(wp),intent(in),optional  :: atol     !! absolute tol for x
    integer,intent(in),optional   :: maxiter  !! maximum number of iterations

    me%f => f
    if (present(ftol))    me%ftol    = abs(ftol)
    if (present(rtol))    me%rtol    = abs(rtol)
    if (present(atol))    me%atol    = abs(atol)
    if (present(maxiter)) me%maxiter = abs(maxiter)

    end subroutine initialize_root_solver
!*****************************************************************************************

!*****************************************************************************************
!>
!  Convert a root method name to the corresponding [[root_method]] enum type.

    pure function root_name_to_method(name) result(r)

    implicit none

    character(len=*),intent(in) :: name !! method name
    type(root_method) :: r !! the [[root_method]] enum type for that method

    integer :: i !! counter
    character(len=len(name)) :: name_lowercase !! lowercase version of `name`

    ! convert to lowercase:
    name_lowercase = lowercase(name)

    ! find the name in the list:
    do i = 1, size(set_of_root_methods)
        if (name_lowercase == set_of_root_methods(i)%name) then
            r = set_of_root_methods(i)
            return
        end if
    end do

    ! if the name was not found
    r = root_method_null

    end function root_name_to_method
!*****************************************************************************************

!*****************************************************************************************
!>
!  Non-object-oriented wrapper.

    subroutine root_scalar_by_name(method,fun,ax,bx,xzero,fzero,iflag,&
                                   ftol,rtol,atol,maxiter,fax,fbx,&
                                   bisect_on_failure)

    implicit none

    character(len=*),intent(in)   :: method   !! the method to use
    procedure(func2)              :: fun      !! user function to find the root of
    real(wp),intent(in)           :: ax       !! left endpoint of initial interval
    real(wp),intent(in)           :: bx       !! right endpoint of initial interval
    real(wp),intent(out)          :: xzero    !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out)          :: fzero    !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)           :: iflag    !! status flag (`-1`=error, `0`=root found, `-999`=invalid method)
    real(wp),intent(in),optional  :: ftol     !! absolute tolerance for `f=0`
    real(wp),intent(in),optional  :: rtol     !! relative tol for x
    real(wp),intent(in),optional  :: atol     !! absolute tol for x
    integer,intent(in),optional   :: maxiter  !! maximum number of iterations
    real(wp),intent(in),optional  :: fax      !! if `f(ax)` is already known, it can be input here
    real(wp),intent(in),optional  :: fbx      !! if `f(bx)` is already known, it can be input here
    logical,intent(in),optional   :: bisect_on_failure  !! if true, then if the specified method fails,
                                                        !! it will be retried using the bisection method.
                                                        !! (default is False). Note that this can use up
                                                        !! to `maxiter` additional function evaluations.

    type(root_method) :: r

    r = root_name_to_method(method)

    if (r%id /= 0) then
        call root_scalar(r,fun,ax,bx,xzero,fzero,iflag,&
                         ftol,rtol,atol,maxiter,fax,fbx,&
                         bisect_on_failure)
    else
        iflag = -999    ! invalid method
        return
    end if

    end subroutine root_scalar_by_name
!*****************************************************************************************

!*****************************************************************************************
!>
!  Non-object-oriented wrapper.

    subroutine root_scalar_by_type(method,fun,ax,bx,xzero,fzero,iflag,&
                                   ftol,rtol,atol,maxiter,fax,fbx,&
                                   bisect_on_failure)

    implicit none

    type(root_method),intent(in)  :: method   !! the method to use
    procedure(func2)              :: fun      !! user function to find the root of
    real(wp),intent(in)           :: ax       !! left endpoint of initial interval
    real(wp),intent(in)           :: bx       !! right endpoint of initial interval
    real(wp),intent(out)          :: xzero    !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out)          :: fzero    !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)           :: iflag    !! status flag (`-1`=error, `0`=root found, `-999`=invalid method)
    real(wp),intent(in),optional  :: ftol     !! absolute tolerance for `f=0`
    real(wp),intent(in),optional  :: rtol     !! relative tol for x
    real(wp),intent(in),optional  :: atol     !! absolute tol for x
    integer,intent(in),optional   :: maxiter  !! maximum number of iterations
    real(wp),intent(in),optional  :: fax      !! if `f(ax)` is already known, it can be input here
    real(wp),intent(in),optional  :: fbx      !! if `f(bx)` is already known, it can be input here
    logical,intent(in),optional   :: bisect_on_failure  !! if true, then if the specified method fails,
                                                        !! it will be retried using the bisection method.
                                                        !! (default is False). Note that this can use up
                                                        !! to `maxiter` additional function evaluations.

    class(root_solver),allocatable :: s

    select case (method%id)

    case(root_method_brent%id);                allocate(brent_solver                :: s)
    case(root_method_bisection%id);            allocate(bisection_solver            :: s)
    case(root_method_regula_falsi%id);         allocate(regula_falsi_solver         :: s)
    case(root_method_illinois%id);             allocate(illinois_solver             :: s)
    case(root_method_anderson_bjorck%id);      allocate(anderson_bjorck_solver      :: s)
    case(root_method_ridders%id);              allocate(ridders_solver              :: s)
    case(root_method_pegasus%id);              allocate(pegasus_solver              :: s)
    case(root_method_bdqrf%id);                allocate(bdqrf_solver                :: s)
    case(root_method_muller%id);               allocate(muller_solver               :: s)
    case(root_method_brenth%id);               allocate(brenth_solver               :: s)
    case(root_method_brentq%id);               allocate(brentq_solver               :: s)
    case(root_method_chandrupatla%id);         allocate(chandrupatla_solver         :: s)
    case(root_method_toms748%id);              allocate(toms748_solver              :: s)
    case(root_method_zhang%id);                allocate(zhang_solver                :: s)
    case(root_method_anderson_bjorck_king%id); allocate(anderson_bjorck_king_solver :: s)
    case(root_method_blendtf%id);              allocate(blendtf_solver              :: s)
    case(root_method_barycentric%id);          allocate(barycentric_solver          :: s)
    case(root_method_itp%id);                  allocate(itp_solver                  :: s)

    case default
        iflag = -999    ! invalid method
        return
    end select

    call s%initialize(func_wrapper,ftol,rtol,atol,maxiter)
    call s%solve(ax,bx,xzero,fzero,iflag,fax,fbx,bisect_on_failure)

    contains

        function func_wrapper(me,x) result(f)
            implicit none
            class(root_solver),intent(inout) :: me
            real(wp),intent(in) :: x
            real(wp) :: f
            f = fun(x)
        end function func_wrapper

    end subroutine root_scalar_by_type
!*****************************************************************************************

!*****************************************************************************************
!>
!  Main wrapper routine for all the methods.

    subroutine solve(me,ax,bx,xzero,fzero,iflag,fax,fbx,bisect_on_failure)

    implicit none

    class(root_solver),intent(inout) :: me
    real(wp),intent(in)              :: ax      !! left endpoint of initial interval
    real(wp),intent(in)              :: bx      !! right endpoint of initial interval
    real(wp),intent(out)             :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out)             :: fzero   !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)              :: iflag   !! status flag (`-1`=error, `0`=root found, `-4`=ax must be /= bx)
    real(wp),intent(in),optional     :: fax     !! if `f(ax)` is already known, it can be input here
    real(wp),intent(in),optional     :: fbx     !! if `f(bx)` is already known, it can be input here
    logical,intent(in),optional      :: bisect_on_failure !! if true, then if the specified method fails,
                                                          !! it will be retried using the bisection method.
                                                          !! (default is False). Note that this can use up
                                                        !! to `maxiter` additional function evaluations.

    real(wp) :: fa !! `f(ax)` passed to the lower level routine
    real(wp) :: fb !! `f(bx)` passed to the lower level routine

    if (ax==bx) then
        ! ax must be /= bx
        iflag = -4
        xzero = ax  ! just to return something
        fzero = fa  !
    else

        call me%get_fa_fb(ax,bx,fax,fbx,fa,fb)

        ! check trivial cases first:
        if (me%solution(ax,fa,xzero,fzero)) then

            iflag = 0

        else if (me%solution(bx,fb,xzero,fzero)) then

            iflag = 0

        else if (fa*fb>0.0_wp) then

            ! f(ax) and f(bx) do not have different signs
            iflag = -1
            xzero = ax  ! just to return something
            fzero = fa  !

        else

            ! call the root solver.
            ! make sure order is correct.
            if (ax<bx) then
                call me%find_root(ax,bx,fa,fb,xzero,fzero,iflag)
            else
                call me%find_root(bx,ax,fb,fa,xzero,fzero,iflag)
            end if

            ! if it failed, then we have the option to then try bisection
            if (iflag /= 0) then
                if (present(bisect_on_failure)) then
                    if (bisect_on_failure) then
                        ! use the wrapper routine for that with the input class
                        call root_scalar(root_method_bisection,&
                                         func_wrapper,ax,bx,xzero,fzero,iflag,&
                                         me%ftol,me%rtol,me%atol,me%maxiter,fa,fb,&
                                         bisect_on_failure = .false.)
                    end if
                end if
            end if

        end if

    end if

    contains

    function func_wrapper(x) result(f)
        !! wrapper function to use bisection
        implicit none
        real(wp),intent(in) :: x
        real(wp) :: f
        f = me%f(x)
    end function func_wrapper

    end subroutine solve
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns the function values at `ax` and `bx` to start the root finding algorithm.

    subroutine get_fa_fb(me,ax,bx,fax,fbx,fa,fb)

    implicit none

    class(root_solver),intent(inout) :: me
    real(wp),intent(in)              :: ax      !! left endpoint of initial interval
    real(wp),intent(in)              :: bx      !! right endpoint of initial interval
    real(wp),intent(in),optional     :: fax     !! if `f(ax)` is already known, it can be input here
    real(wp),intent(in),optional     :: fbx     !! if `f(bx)` is already known, it can be input here
    real(wp),intent(out)             :: fa      !! `f(ax)` to use
    real(wp),intent(out)             :: fb      !! `f(ax)` to use

    if (present(fax)) then
        fa = fax
    else
        fa = me%f(ax)
    end if

    if (present(fbx)) then
        fb = fbx
    else
        fb = me%f(bx)
    end if

    end subroutine get_fa_fb
!*****************************************************************************************

!*****************************************************************************************
!>
!  Find a zero of the function \( f(x) \) in the given interval
!  \( [a_x,b_x] \) to within a tolerance \( 4 \epsilon |x| + tol \),
!  where \( \epsilon \) is the relative machine precision defined as
!  the smallest representable number such that \( 1.0 + \epsilon > 1.0 \).
!
!### References
!  * R. P. Brent, "[An algorithm with guaranteed convergence for
!    finding a zero of a function](http://maths-people.anu.edu.au/~brent/pd/rpb005.pdf)",
!    The Computer Journal, Vol 14, No. 4., 1971.
!  * R. P. Brent, "[Algorithms for minimization without derivatives](http://maths-people.anu.edu.au/~brent/pub/pub011.html)",
!    Prentice-Hall, Inc., 1973.
!
!### See also
!  * [zeroin.f](http://www.netlib.org/go/zeroin.f) from Netlib

    subroutine brent(me,ax,bx,fax,fbx,xzero,fzero,iflag)

    implicit none

    class(brent_solver),intent(inout) :: me
    real(wp),intent(in)    :: ax    !! left endpoint of initial interval
    real(wp),intent(in)    :: bx    !! right endpoint of initial interval
    real(wp),intent(in)    :: fax   !! `f(ax)`
    real(wp),intent(in)    :: fbx   !! `f(ax)`
    real(wp),intent(out)   :: xzero !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out)   :: fzero !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)    :: iflag !! status flag (`0`=root found)

    real(wp),parameter :: eps = epsilon(1.0_wp)  !! d1mach(4) in original code
    real(wp) :: a,b,c,d,e,fa,fb,fc,tol1,xm,p,q,r,s
    integer :: i !! iteration counter

    ! initialize:
    iflag = 0
    tol1  = eps+1.0_wp
    a     = ax
    b     = bx
    fa    = fax
    fb    = fbx
    c     = a
    fc    = fa
    d     = b-a
    e     = d

    do i=1,me%maxiter

        if (abs(fc)<abs(fb)) then
            a=b
            b=c
            c=a
            fa=fb
            fb=fc
            fc=fa
        end if

        tol1 = 2.0_wp*eps*abs(b)+0.5_wp*me%rtol
        xm = 0.5_wp*(c-b)
        if (abs(xm)<=tol1) exit

        ! see if a bisection is forced
        if ((abs(e)>=tol1) .and. (abs(fa)>abs(fb))) then
            s=fb/fa
            if (a/=c) then
                ! inverse quadratic interpolation
                q=fa/fc
                r=fb/fc
                p=s*(2.0_wp*xm*q*(q-r)-(b-a)*(r-1.0_wp))
                q=(q-1.0_wp)*(r-1.0_wp)*(s-1.0_wp)
            else
                ! linear interpolation
                p=2.0_wp*xm*s
                q=1.0_wp-s
            end if
            if (p<=0.0_wp) then
                p=-p
            else
                q=-q
            end if
            s=e
            e=d
            if (((2.0_wp*p)>=(3.0_wp*xm*q-abs(tol1*q))) .or. (p>=abs(0.5_wp*s*q))) then
                d=xm
                e=d
            else
                d=p/q
            end if
        else
            d=xm
            e=d
        end if

        a=b
        fa=fb
        if (abs(d)<=tol1) then
            if (xm<=0.0_wp) then
                b=b-tol1
            else
                b=b+tol1
            end if
        else
            b=b+d
        end if
        fb=me%f(b)
        if (abs(fb)<=me%ftol) exit  ! absolute convergence in f
        if ((fb*(fc/abs(fc)))>0.0_wp) then
            c=a
            fc=fa
            d=b-a
            e=d
        end if

        if (i==me%maxiter) iflag = -2  ! max iterations reached

    end do

    xzero = b
    fzero = fb

    end subroutine brent
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the zero of the function f(x) in the interval ax,bx using the bisection method.
!
!### See also
!  * G.E. Mullges & F. Uhlig, "Numerical Algorithms with Fortran",
!    Springer, 1996. Section 2.8.1, p 32-34.

    subroutine bisection(me,ax,bx,fax,fbx,xzero,fzero,iflag)

    implicit none

    class(bisection_solver),intent(inout) :: me
    real(wp),intent(in)    :: ax      !! left endpoint of initial interval
    real(wp),intent(in)    :: bx      !! right endpoint of initial interval
    real(wp),intent(in)    :: fax     !! `f(ax)`
    real(wp),intent(in)    :: fbx     !! `f(ax)`
    real(wp),intent(out)   :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out)   :: fzero   !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)    :: iflag   !! status flag (`0`=root found, `-2`=max iterations reached)

    real(wp) :: x1,x2,x3,f1,f2,f3
    integer :: i !! iteration counter
    logical :: root_found !! convergence in x

    ! initialize:
    iflag = 0
    x1    = ax
    x2    = bx
    f1    = fax
    f2    = fbx

    ! main loop
    do i=1,me%maxiter

        ! bisection of the inclusion interval:
        !  x1------x3------x2
        x3 = bisect(x1,x2)

        ! calculate the new function value:
        f3 = me%f(x3)
        if (me%solution(x3,f3,xzero,fzero)) return

        ! determine new inclusion interval:
        if (f2*f3<0.0_wp) then
            ! root lies between x2 and x3
            x1 = x3
            x2 = x2
            f1 = f3
            f2 = f2
        else
            ! root lies between x1 and x3
            x2 = x3
            f2 = f3
        end if

        ! check for convergence:
        root_found = me%converged(x1,x2)
        if (root_found .or. i==me%maxiter) then
            call choose_best(x1,x2,f1,f2,xzero,fzero)
            if (.not. root_found) iflag = -2  ! max iterations reached
            exit
        end if

    end do

    end subroutine bisection
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the zero of the function f(x) in the interval ax,bx using the regula falsi method.

    subroutine regula_falsi(me,ax,bx,fax,fbx,xzero,fzero,iflag)

    implicit none

    class(regula_falsi_solver),intent(inout) :: me
    real(wp),intent(in)    :: ax      !! left endpoint of initial interval
    real(wp),intent(in)    :: bx      !! right endpoint of initial interval
    real(wp),intent(in)    :: fax     !! `f(ax)`
    real(wp),intent(in)    :: fbx     !! `f(ax)`
    real(wp),intent(out)   :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out)   :: fzero   !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)    :: iflag   !! status flag (`0`=root found, `-2`=max iterations reached)

    real(wp) :: x1,x2,x3,f1,f2,f3
    integer :: i !! iteration counter
    logical :: root_found !! convergence in x

    ! initialize:
    iflag = 0
    x1    = ax
    x2    = bx
    f1    = fax
    f2    = fbx

    ! main loop
    do i=1,me%maxiter

        x3 = regula_falsi_step(x1,x2,f1,f2,ax,bx)

        ! calculate the new function value:
        f3 = me%f(x3)
        if (me%solution(x3,f3,xzero,fzero)) return

        ! determine new inclusion interval:
        if (f2*f3<0.0_wp) then
            ! root lies between x2 and x3
            x1 = x3
            x2 = x2
            f1 = f3
            f2 = f2
        else
            ! root lies between x1 and x3
            x2 = x3
            f2 = f3
        end if

        ! check for convergence:
        root_found = me%converged(x1,x2)
        if (root_found .or. i==me%maxiter) then
            call choose_best(x1,x2,f1,f2,xzero,fzero)
            if (.not. root_found) iflag = -2  ! max iterations reached
            exit
        end if

    end do

    end subroutine regula_falsi
!*****************************************************************************************

!*****************************************************************************************
!>
!  Illinois method.
!
!### Reference
!  * M. Dowell, P. Jarratt, "A modified regula falsi method for computing the root
!    of an equation', BIT 11 (1971), 168-174.

    subroutine illinois(me,ax,bx,fax,fbx,xzero,fzero,iflag)

    implicit none

    class(illinois_solver),intent(inout) :: me
    real(wp),intent(in)    :: ax      !! left endpoint of initial interval
    real(wp),intent(in)    :: bx      !! right endpoint of initial interval
    real(wp),intent(in)    :: fax     !! `f(ax)`
    real(wp),intent(in)    :: fbx     !! `f(ax)`
    real(wp),intent(out)   :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out)   :: fzero   !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)    :: iflag   !! status flag (`0`=root found, `-2`=max iterations reached)

    real(wp) :: x1,x2,x3,f1,f2,f3,delta,f1tmp
    integer :: i !! iteration counter
    logical :: root_found !! convergence in x

    ! initialize:
    iflag = 0
    x1    = ax
    x2    = bx
    f1    = fax
    f2    = fbx

    ! main loop
    do i=1,me%maxiter

        x3 = regula_falsi_step(x1,x2,f1,f2,ax,bx)

        ! calculate the new function value:
        f3 = me%f(x3)
        if (me%solution(x3,f3,xzero,fzero)) return

        ! determine new inclusion interval:
        if (f2*f3<0.0_wp) then
            ! root lies between x2 and x3
            x1 = x2
            x2 = x3
            f1 = f2
            f1tmp = f1
            f2 = f3
        else
            ! root lies between x1 and x3
            x2 = x3
            f2 = f3
            f1tmp = f1 ! actual function eval
            f1 = 0.5_wp * f1
        end if

        ! check for convergence:
        root_found = me%converged(x1,x2)
        if (root_found .or. i==me%maxiter) then
            call choose_best(x1,x2,f1tmp,f2,xzero,fzero)
            if (.not. root_found) iflag = -2  ! max iterations reached
            exit
        end if

    end do

    end subroutine illinois
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the zero of the function f(x) in the interval ax,bx using the Anderson-Bjorck method.
!
!### See also
!  * G.E. Mullges & F. Uhlig, "Numerical Algorithms with Fortran",
!    Springer, 1996. Section 2.8.2, p 36.

    subroutine anderson_bjorck(me,ax,bx,fax,fbx,xzero,fzero,iflag)

    implicit none

    class(anderson_bjorck_solver),intent(inout) :: me
    real(wp),intent(in)    :: ax      !! left endpoint of initial interval
    real(wp),intent(in)    :: bx      !! right endpoint of initial interval
    real(wp),intent(in)    :: fax     !! `f(ax)`
    real(wp),intent(in)    :: fbx     !! `f(ax)`
    real(wp),intent(out)   :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out)   :: fzero   !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)    :: iflag   !! status flag (`0`=root found, `-2`=max iterations reached)

    integer :: i !! counter
    logical :: root_found !! convergence in x
    real(wp) :: x1,x2,x3,f1,f2,f3,g,f1tmp

    ! initialize:
    iflag = 0
    x1    = ax
    x2    = bx
    f1    = fax
    f2    = fbx

    ! main loop:
    do i = 1,me%maxiter

        x3 = secant(x1,x2,f1,f2,ax,bx)
        f3 = me%f(x3)
        if (me%solution(x3,f3,xzero,fzero)) return

        ! determine a new inclusion interval:
        if (f2*f3<0.0_wp) then
            ! zero lies between x2 and x3
            x1 = x2
            x2 = x3
            f1 = f2
            f2 = f3
            f1tmp = f1
        else
            ! zero lies between x1 and x3
            g = 1.0_wp - f3/f2
            if (g<=0.0_wp) g = 0.5_wp
            x2 = x3
            f1tmp = f1
            f1 = g*f1
            f2 = f3
        end if

        ! check for convergence:
        root_found = me%converged(x1,x2)
        if (root_found .or. i == me%maxiter) then
            call choose_best(x1,x2,f1tmp,f2,xzero,fzero)
            if (.not. root_found) iflag = -2  ! max iterations reached
            exit
        end if

    end do

    end subroutine anderson_bjorck
!*****************************************************************************************

!*****************************************************************************************
!>
!  Ridders method to find a root of f(x).
!
!### See also
!  * Ridders, C., "A new algorithm for computing a single root of a real continuous function",
!    IEEE Trans. on Circuits and Systems, Vol 26, Issue 11, Nov 1979.

    subroutine ridders(me,ax,bx,fax,fbx,xzero,fzero,iflag)

    implicit none

    class(ridders_solver),intent(inout) :: me
    real(wp),intent(in)  :: ax      !! left endpoint of initial interval
    real(wp),intent(in)  :: bx      !! right endpoint of initial interval
    real(wp),intent(in)  :: fax     !! `f(ax)`
    real(wp),intent(in)  :: fbx     !! `f(ax)`
    real(wp),intent(out) :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out) :: fzero   !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)  :: iflag   !! status flag (`0`=root found, `-2`=max iterations reached, `-3`=singularity in the algorithm)

    integer  :: i !! counter
    real(wp) :: fh,fl,fm,fnew,denom,xh,xl,xm,xnew

    ! initialize:
    iflag = 0
    fl    = fax
    fh    = fbx
    xl    = ax
    xh    = bx
    xzero = huge(1.0_wp)

    do i = 1, me%maxiter

        xm = bisect(xl,xh)
        fm = me%f(xm)
        if (me%solution(xm,fm,xzero,fzero)) return

        denom = sqrt(fm**2-fl*fh)
        if (denom == 0.0_wp) then
            xzero = xm
            fzero = fm
            iflag = -3        ! can't proceed: denominator is zero [TODO: add a bisection if this happens]
            exit
        end if

        xnew = xm + (xm-xl)*(sign(1.0_wp,fl-fh)*fm/denom)
        if (me%converged(xzero,xnew)) then  ! relative convergence in x
            ! additional check to prevent false convergence
            if (me%converged(xl,xm) .or. me%converged(xm,xh)) exit
        end if

        xzero = xnew
        fnew  = me%f(xzero)
        fzero = fnew
        if (abs(fnew) <= me%ftol) exit    ! abs convergence in f

        ! to keep the root bracketed:
        if (sign(fm,fnew) /= fm) then
            xl = xm
            fl = fm
            xh = xzero
            fh = fnew
        else if (sign(fl,fnew) /= fl) then
            xh = xzero
            fh = fnew
        else if (sign(fh,fnew) /= fh) then
            xl = xzero
            fl = fnew
        end if

        if (me%converged(xl,xh)) exit    ! relative convergence in x
        if (i == me%maxiter) iflag = -2  ! max iterations exceeded

    end do

    end subroutine ridders
!*****************************************************************************************

!*****************************************************************************************
!>
!  Pegasus method to find a root of f(x).
!
!### See also
!  * G.E. Mullges & F. Uhlig, "Numerical Algorithms with Fortran",
!    Springer, 1996. Section 2.8.2, p 35.

    subroutine pegasus(me,ax,bx,fax,fbx,xzero,fzero,iflag)

    implicit none

    class(pegasus_solver),intent(inout) :: me
    real(wp),intent(in)  :: ax      !! left endpoint of initial interval
    real(wp),intent(in)  :: bx      !! right endpoint of initial interval
    real(wp),intent(in)  :: fax     !! `f(ax)`
    real(wp),intent(in)  :: fbx     !! `f(ax)`
    real(wp),intent(out) :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out) :: fzero   !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)  :: iflag   !! status flag (`0`=root found, `-2`=max iterations reached)

    integer :: i !! counter
    real(wp) :: x1,x2,x3,f1,f2,f3,f1tmp,denom

    ! initialize:
    iflag = 0
    x1    = ax
    x2    = bx
    f1    = fax
    f2    = fbx

    ! main loop:
    do i = 1, me%maxiter

        ! secant step
        x3 = secant(x1,x2,f1,f2,ax,bx)

        f3  = me%f(x3)  ! calculate f3
        if (me%solution(x3,f3,xzero,fzero)) return

        ! determine a new inclusion interval:
        if (f2*f3<0.0_wp) then  ! root on (x2,x3)
            x1 = x2
            f1 = f2
            f1tmp = f1
        else  ! root on (x1,x3)
            f1tmp = f1
            denom = f2 + f3
            if (denom /= 0.0_wp) then
                ! proceed as normal
                f1 = f1 * f2 / denom
            else
                ! can't proceed, keep as is.
                ! [need a find a test case where this happens -TODO]
            end if
        end if

        x2 = x3
        f2 = f3

        call choose_best(x1,x2,f1tmp,f2,xzero,fzero)

        if (me%converged(x1,x2)) exit   ! check for convergence
        if (i == me%maxiter) iflag = -2 ! max iterations exceeded

    end do

    end subroutine pegasus
!*****************************************************************************************

!*****************************************************************************************
!>
!  Bisected Direct Quadratic Regula Falsi (BDQRF) root solver method
!  to find the root of a 1D function.
!
!### See also
!  * R. G. Gottlieb, B. F. Thompson, "Bisected Direct Quadratic Regula Falsi",
!    Applied Mathematical Sciences, Vol. 4, 2010, no. 15, 709-718.

    subroutine bdqrf(me,ax,bx,fax,fbx,xzero,fzero,iflag)

    implicit none

    class(bdqrf_solver),intent(inout) :: me
    real(wp),intent(in)  :: ax      !! left endpoint of initial interval
    real(wp),intent(in)  :: bx      !! right endpoint of initial interval
    real(wp),intent(in)  :: fax     !! `f(ax)`
    real(wp),intent(in)  :: fbx     !! `f(ax)`
    real(wp),intent(out) :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out) :: fzero   !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)  :: iflag   !! status flag (`0`=root found, `-2`=max iterations reached)

    real(wp) :: xdn,ydn,xup,yup,d,xm,ym,a,b,y2
    integer :: i !! counter

    ! initialize:
    iflag = 0
    xzero = ax
    fzero = fax
    y2    = fbx

    if (fzero<0.0_wp) then
        xdn = ax
        ydn = fzero
        xup = bx
        yup = y2
    else
        xup = ax
        yup = fzero
        xdn = bx
        ydn = y2
    end if

    ! main loop:
    do i = 1, me%maxiter

        xm = bisect(xup,xdn)
        ym = me%f(xm)
        if (me%solution(xm,ym,xzero,fzero)) return ! Convergence

        d = (xup - xdn) / 2.0_wp
        a = (yup + ydn - 2.0_wp*ym)/(2.0_wp*d**2)
        b = (yup - ydn)/(2.0_wp*d)

        xzero = xm - 2.0_wp*ym / (b * (1.0_wp + sqrt(1.0_wp - 4.0_wp*a*ym/b**2)))
        fzero = me%f(xzero)
        if (me%solution(xzero,fzero,xzero,fzero)) return ! Convergence

        if (fzero>0.0_wp) then
            yup = fzero
            xup = xzero
            if (ym<0.0_wp) then
                ydn = ym
                xdn = xm
            end if
        else
            ydn = fzero
            xdn = xzero
            if (ym>0.0_wp) then
                yup = ym
                xup = xm
            end if
        end if

        if (me%converged(xdn,xup) .or. i==me%maxiter) then
            call choose_best(xdn,xup,ydn,yup,xzero,fzero)
            if (i==me%maxiter) iflag = -2 ! maximum number of iterations
            exit
        end if

    end do

    end subroutine bdqrf
!*****************************************************************************************

!*****************************************************************************************
!>
!  Improved Muller method (for real roots only).
!  Will fall back to bisection if any step fails.
!
!### Reference
!  * D. E. Muller, "A Method for Solving Algebraic Equations Using an Automatic Computer",
!    Mathematical Tables and Other Aids to Computation, 10 (1956), 208-215.
!    [link](https://www.ams.org/journals/mcom/1956-10-056/S0025-5718-1956-0083822-0/S0025-5718-1956-0083822-0.pdf)
!  * [Roots.jl](https://github.com/JuliaMath/Roots.jl/blob/master/src/simple.jl), 
!    Julia version of standard Muller

    subroutine muller (me,ax,bx,fax,fbx,xzero,fzero,iflag)

    implicit none

    class(muller_solver),intent(inout) :: me
    real(wp),intent(in)  :: ax      !! left endpoint of initial interval
    real(wp),intent(in)  :: bx      !! right endpoint of initial interval
    real(wp),intent(in)  :: fax     !! `f(ax)`
    real(wp),intent(in)  :: fbx     !! `f(ax)`
    real(wp),intent(out) :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out) :: fzero   !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)  :: iflag   !! status flag (`0`=root found, `-2`=max iterations reached)

    real(wp) :: a,b,c,cx,fa,fb,fc,fcx,x,q,q2,q1,aa,bb,cc,delta,dp,dm,denon,bprev,fbprev
    integer  :: i    !! iteration counter
    logical  :: x_ok !! the new estimate `x` is ok to use

    iflag = 0

    ! pick a third point in the middle [this could also be an optional input]
    cx  = bisect(ax,bx)
    fcx = me%f(cx)
    if (me%solution(cx,fcx,xzero,fzero)) return

    ! [a,b,c]
    a = ax; fa = fax
    b = cx; fb = fcx
    c = bx; fc = fbx

    bprev = huge(1.0_wp)
    fbprev = huge(1.0_wp)

    do i = 1, me%maxiter

        ! muller step:
        q     = (c - b)/(b - a)
        q2    = q**2
        q1    = q + 1.0_wp
        aa    = q*fc - q*q1*fb + q2*fa
        bb    = (q1+q)*fc - q1**2*fb + q2*fa
        cc    = q1*fc
        delta = sqrt(max(0.0_wp, bb**2 - 4.0_wp * aa*cc)) ! to avoid complex roots
        dp    = bb + delta
        dm    = bb - delta
        if (abs(dp) > abs(dm)) then
            denon = dp
        else
            denon = dm
        end if
        x_ok = denon /= 0.0_wp
        if (x_ok) x = c - 2.0_wp*(c - b)*cc/denon

        ! make sure that x is ok, in the correct interval, and distinct.
        ! if not, fall back to bisection on that interval
        if (fa*fb < 0.0_wp) then  ! root in (a,b)
            if (.not. x_ok .or. x<=a .or. x>=b) x = bisect(a,b)
            c  = b
            fc = fb
            b  = x
        else  ! root in (b,c)
            if (.not. x_ok .or. x<=b .or. x>=c) x = bisect(b,c)
            a  = b
            fa = fb
            b  = x
        end if
        ! values are now [a,b,c], with b being the new estimate

        ! function evaluation for next estimate:
        fb = me%f(b)
        if (abs(fb)<=me%ftol) exit

        ! stopping criterion
        if (me%converged(a,c) .or. i == me%maxiter) then
            if ( i == me%maxiter ) iflag = -2 ! max iterations exceeded
            exit
        end if

        bprev = b
        fbprev = fb

    end do

    call choose_best(b,bprev,fb,fbprev,xzero,fzero)

    end subroutine muller
!*****************************************************************************************

!*****************************************************************************************
!>
!  Brent's method with hyperbolic extrapolation.
!
!  A variation on the classic Brent routine to find a zero of the function f
!  between the arguments ax and bx that uses hyperbolic extrapolation instead
!  of inverse quadratic extrapolation.
!
!### Reference
!  * SciPy `brenth.c`

    subroutine brenth(me,ax,bx,fax,fbx,xzero,fzero,iflag)

    implicit none

    class(brenth_solver),intent(inout) :: me
    real(wp),intent(in)  :: ax      !! left endpoint of initial interval
    real(wp),intent(in)  :: bx      !! right endpoint of initial interval
    real(wp),intent(in)  :: fax     !! `f(ax)`
    real(wp),intent(in)  :: fbx     !! `f(ax)`
    real(wp),intent(out) :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out) :: fzero   !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)  :: iflag   !! status flag (`0`=root found, `-2`=max iterations reached)

    real(wp) :: xpre,xcur,xblk,fpre,fcur,fblk,spre,&
                scur,sbis,delta,stry,dpre,dblk,xdelta
    integer :: i !! iteration counter

    iflag = 0
    xpre = ax
    xcur = bx
    fpre = fax
    fcur = fbx

    do i = 1, me%maxiter

        if (fpre*fcur < 0.0_wp) then
            xblk = xpre
            fblk = fpre
            scur = xcur - xpre
            spre = scur
        end if
        if (abs(fblk) < abs(fcur)) then
            xpre = xcur
            xcur = xblk
            xblk = xpre
            fpre = fcur
            fcur = fblk
            fblk = fpre
        end if

        delta = (me%atol + me%rtol*abs(xcur))/2.0_wp
        sbis = (xblk - xcur)/2.0_wp
        if (abs(fcur) <= me%ftol .or. abs(sbis) < delta) exit ! converged

        if (abs(spre) > delta .and. abs(fcur) < abs(fpre)) then
            if (xpre == xblk) then
                ! interpolate
                stry = -fcur*(xcur - xpre)/(fcur - fpre)
            else
                ! extrapolate
                dpre = (fpre - fcur)/(xpre - xcur)
                dblk = (fblk - fcur)/(xblk - xcur)
                stry = -fcur*(fblk - fpre)/(fblk*dpre - fpre*dblk)  ! only difference from brentq
            end if

            if (2.0_wp*abs(stry) < min(abs(spre), 3.0_wp*abs(sbis) - delta)) then
                ! accept step
                spre = scur
                scur = stry
            else
                ! bisect
                spre = sbis
                scur = sbis
            end if
        else
            ! bisect
            spre = sbis
            scur = sbis
        end if

        xpre = xcur
        fpre = fcur
        if (abs(scur) > delta) then
            xcur = xcur + scur
        else
            if (sbis > 0.0_wp) then
                xdelta = delta
            else
                xdelta = -delta
            end if
            xcur = xcur + xdelta
        end if

        fcur = me%f(xcur)
        if (abs(fcur) <= me%ftol) exit ! converged
        if (i == me%maxiter) iflag = -2 ! max iterations reached

    end do

    xzero = xcur
    fzero = fcur

    end subroutine brenth
!*****************************************************************************************

!*****************************************************************************************
!>
!  Classic Brent's method to find a zero of the function f on the sign
!  changing interval [ax, bx], but with a different formula for the extrapolation step.
!
!### Reference
!  * SciPy brentq.c

    subroutine brentq(me,ax,bx,fax,fbx,xzero,fzero,iflag)

    implicit none

    class(brentq_solver),intent(inout) :: me
    real(wp),intent(in)  :: ax      !! left endpoint of initial interval
    real(wp),intent(in)  :: bx      !! right endpoint of initial interval
    real(wp),intent(in)  :: fax     !! `f(ax)`
    real(wp),intent(in)  :: fbx     !! `f(ax)`
    real(wp),intent(out) :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out) :: fzero   !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)  :: iflag   !! status flag (`0`=root found, `-2`=max iterations reached)

    real(wp) :: xpre,xcur,xblk,fpre,fcur,fblk,spre,&
                scur,sbis,delta,stry,dpre,dblk,xdelta
    integer :: i !! iteration counter

    iflag = 0
    xpre = ax
    xcur = bx
    fpre = fax
    fcur = fbx

    do i = 1, me%maxiter

        if (fpre*fcur < 0.0_wp) then
            xblk = xpre
            fblk = fpre
            scur = xcur - xpre
            spre = scur
        end if
        if (abs(fblk) < abs(fcur)) then
            xpre = xcur
            xcur = xblk
            xblk = xpre
            fpre = fcur
            fcur = fblk
            fblk = fpre
        end if

        delta = (me%atol + me%rtol*abs(xcur))/2.0_wp
        sbis = (xblk - xcur)/2.0_wp
        if (abs(fcur) <= me%ftol .or. abs(sbis) < delta) exit ! converged

        if (abs(spre) > delta .and. abs(fcur) < abs(fpre)) then
            if (xpre == xblk) then
                ! interpolate
                stry = -fcur*(xcur - xpre)/(fcur - fpre)
            else
                ! extrapolate
                dpre = (fpre - fcur)/(xpre - xcur)
                dblk = (fblk - fcur)/(xblk - xcur)
                stry = -fcur*(fblk*dblk - fpre*dpre)/(dblk*dpre*(fblk - fpre))  ! only difference from brenth
            end if

            if (2.0_wp*abs(stry) < min(abs(spre), 3.0_wp*abs(sbis) - delta)) then
                ! accept step
                spre = scur
                scur = stry
            else
                ! bisect
                spre = sbis
                scur = sbis
            end if
        else
            ! bisect
            spre = sbis
            scur = sbis
        end if

        xpre = xcur
        fpre = fcur
        if (abs(scur) > delta) then
            xcur = xcur + scur
        else
            if (sbis > 0.0_wp) then
                xdelta = delta
            else
                xdelta = -delta
            end if
            xcur = xcur + xdelta
        end if

        fcur = me%f(xcur)
        if (abs(fcur) <= me%ftol) exit ! converged
        if (i == me%maxiter) iflag = -2 ! max iterations reached

    end do

    xzero = xcur
    fzero = fcur

    end subroutine brentq
!*****************************************************************************************

!*****************************************************************************************
!>
!  Chandrupatla's method.
!
!### Reference
!  * T.R. Chandrupatla, "A new hybrid quadratic/bisection algorithm for
!    finding the zero of a nonlinear function without derivatives," Advances in
!    Engineering Software, Vol 28, 1997, pp. 145-149.
!  * P. Scherer, "Computational Physics: Simulation of Classical and Quantum Systems",
!    Section 6.1.7.3. [this routine was coded from that description]
!  * Python version: https://www.embeddedrelated.com/showarticle/855.php

    subroutine chandrupatla(me,ax,bx,fax,fbx,xzero,fzero,iflag)

    implicit none

    class(chandrupatla_solver),intent(inout) :: me
    real(wp),intent(in)  :: ax      !! left endpoint of initial interval
    real(wp),intent(in)  :: bx      !! right endpoint of initial interval
    real(wp),intent(in)  :: fax     !! `f(ax)`
    real(wp),intent(in)  :: fbx     !! `f(ax)`
    real(wp),intent(out) :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out) :: fzero   !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)  :: iflag   !! status flag (`0`=root found, `-2`=max iterations reached)

    real(wp) :: a,b,c,fa,fb,fc,t,xt,ft,tol,tl,xi,phi,xm,fm
    integer :: i !! iteration counter

    ! initialization:
    iflag = 0
    b  = ax
    a  = bx
    c  = bx
    fa = fbx
    fb = fax
    fc = fb
    t  = 0.5_wp

    ! main loop:
    do i = 1, me%maxiter

        xt = a + t*(b-a)
        ft = me%f(xt)
        if (me%solution(xt,ft,xzero,fzero)) return

        if (ft*fa>0.0_wp) then
            c = a
            fc = fa
        else
            c = b
            b = a
            fc = fb
            fb = fa
        end if
        a = xt
        fa = ft

        if (abs(fb) < abs(fa)) then
            xm = b
            fm = fb
        else
            xm = a
            fm = fa
        end if

        if (i == me%maxiter) then
            iflag = -2 ! max iterations reached
            exit
        end if

        tol = 2.0_wp*me%rtol*abs(xm) + me%atol
        tl = tol/abs(b-c)
        if (tl > 0.5_wp) exit
        t = 0.5_wp ! use bisection unless we can use inverse quadratic below

        if (fa/=fb .and. fb/=fc) then
            xi  = (a-b)/(c-b)
            phi = (fa-fb)/(fc-fb)
            if (1.0_wp - sqrt(1.0_wp - xi) < phi .and. phi < sqrt(xi)) then
                ! inverse quadratic interpolation
                t = (fa/(fb-fa)) * (fc/(fb-fc)) + ((c-a)/(b-a)) * (fa/(fc-fa)) * (fb/(fc-fb))
            end if
        end if

        t = min(1.0_wp-tl, max(tl, t))

    end do

    xzero = xm
    fzero = fm

    end subroutine chandrupatla
!*****************************************************************************************

!*****************************************************************************************
!>
!  TOMS748 rootfinding method.
!
!  Finds either an exact solution or an approximate solution of the
!  equation `f(x)=0` in the interval [ax,bx]. At the begining of each
!  iteration, the current enclosing interval is recorded as [a0,b0].
!  The first iteration is simply a secant step. Starting with the
!  second iteration, three steps are taken in each iteration. First
!  two steps are either quadratic interpolation or cubic inverse
!  interpolation. The third step is a double-size secant step. If the
!  diameter of the enclosing interval obtained after those three steps
!  is larger than 0.5*(b0-a0), then an additional bisection step will
!  be taken.
!
!### References
!  * http://www.netlib.org/toms/748
!  * G. E. Alefeld, F. A. Potra and Yixun Shi,
!    "Algorithm 748: Enclosing Zeros of Continuous Functions",
!    ACM Transactions on Mathematical Software,
!    Vol. 21. No. 3. September 1995. Pages 327-344.

    subroutine toms748(me,ax,bx,fax,fbx,xzero,fzero,iflag)

    implicit none

    class(toms748_solver),intent(inout) :: me
    real(wp),intent(in)  :: ax      !! left endpoint of initial interval
    real(wp),intent(in)  :: bx      !! right endpoint of initial interval
    real(wp),intent(in)  :: fax     !! `f(ax)`
    real(wp),intent(in)  :: fbx     !! `f(ax)`
    real(wp),intent(out) :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out) :: fzero   !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)  :: iflag   !! status flag (`0`=root found, `-2`=max iterations reached)

    integer  :: itnum
    real(wp) :: a,b,fa,fb,c,u,fu,a0,b0,tol,d,fd
    real(wp) :: prof,e,fe,tmpc

    a = ax
    b = bx
    fa = fax
    fb = fbx

    ! initialization. set the number of iteration as 0.
    ! set dumb values for the variables "e" and "fe".
    e  = huge(1.0_wp)
    fe = huge(1.0_wp)

    ! iteration starts. the enclosing interval before executing the
    ! iteration is recorded as [a0, b0].
    do itnum = 1, me%maxiter

        a0 = a
        b0 = b

        ! calculates the termination criterion. stops the procedure if the
        ! criterion is satisfied.
        if (abs(fb) <= abs(fa)) then
            tol = get_tolerance(b)
        else
            tol = get_tolerance(a)
        end if
        if ((b-a)<=tol) exit

        ! for the first iteration, secant step is taken.
        if (itnum == 1) then

            c=a-(fa/(fb-fa))*(b-a)

            ! call subroutine "bracket" to get a shrinked enclosing interval as
            ! well as to update the termination criterion. stop the procedure
            ! if the criterion is satisfied or the exact solution is obtained.
            call bracket(a,b,c,fa,fb,tol,d,fd)
            if ((abs(fa)<=me%ftol) .or. ((b-a)<=tol)) exit

            cycle

        end if

        ! starting with the second iteration, in the first two steps, either
        ! quadratic interpolation is used by calling the subroutine "newqua"
        ! or the cubic inverse interpolation is used by calling the subroutine
        ! "pzero". in the following, if "prof" is not equal to 0, then the
        ! four function values "fa", "fb", "fd", and "fe" are distinct, and
        ! hence "pzero" will be called.
        prof=(fa-fb)*(fa-fd)*(fa-fe)*(fb-fd)*(fb-fe)*(fd-fe)
        if ((itnum == 2) .or. (prof == 0.0_wp)) then
            call newqua(a,b,d,fa,fb,fd,c,2)
        else
            c = pzero(a,b,d,e,fa,fb,fd,fe)
            if ((c-a)*(c-b) >= 0.0_wp) then
                call newqua(a,b,d,fa,fb,fd,c,2)
            end if
        end if
        e=d
        fe=fd

        ! call subroutine "bracket" to get a shrinked enclosing interval as
        ! well as to update the termination criterion. stop the procedure
        ! if the criterion is satisfied or the exact solution is obtained.
        call bracket(a,b,c,fa,fb,tol,d,fd)
        if ((abs(fa)<=me%ftol) .or. ((b-a)<=tol)) exit
        prof=(fa-fb)*(fa-fd)*(fa-fe)*(fb-fd)*(fb-fe)*(fd-fe)
        if (prof == 0.0_wp) then
            call newqua(a,b,d,fa,fb,fd,c,3)
        else
            c = pzero(a,b,d,e,fa,fb,fd,fe)
            if ((c-a)*(c-b) >= 0.0_wp) then
                call newqua(a,b,d,fa,fb,fd,c,3)
            end if
        end if

        ! call subroutine "bracket" to get a shrinked enclosing interval as
        ! well as to update the termination criterion. stop the procedure
        ! if the criterion is satisfied or the exact solution is obtained.
        call bracket(a,b,c,fa,fb,tol,d,fd)
        if ((abs(fa)<=me%ftol) .or. ((b-a)<=tol)) exit
        e=d
        fe=fd

        ! takes the double-size secant step.
        if (abs(fa) < abs(fb)) then
            u=a
            fu=fa
        else
            u=b
            fu=fb
        end if
        c=u-2.0_wp*(fu/(fb-fa))*(b-a)
        if (abs(c-u) > (0.5_wp*(b-a))) then
            c=a+0.5_wp*(b-a)
        end if

        ! call subroutine bracket to get a shrinked enclosing interval as
        ! well as to update the termination criterion. stop the procedure
        ! if the criterion is satisfied or the exact solution is obtained.
        call bracket(a,b,c,fa,fb,tol,d,fd)
        if ((abs(fa)<=me%ftol) .or. ((b-a)<=tol)) exit

        ! determines whether an additional bisection step is needed. and takes
        ! it if necessary.
        if ((b-a) < (0.5_wp*(b0-a0))) cycle
        e=d
        fe=fd

        ! call subroutine "bracket" to get a shrinked enclosing interval as
        ! well as to update the termination criterion. stop the procedure
        ! if the criterion is satisfied or the exact solution is obtained.
        tmpc = a+0.5_wp*(b-a)
        call bracket(a,b,tmpc,fa,fb,tol,d,fd)
        if ((abs(fa)<=me%ftol) .or. ((b-a)<=tol)) exit

        if (itnum == me%maxiter) iflag = -2    ! maximum iterations reached

    end do

    !return result:
    xzero = a
    fzero = fa

    contains
!**********************************************************************************

    !************************************************************************
    subroutine bracket(a,b,c,fa,fb,tol,d,fd)

    !!  Given current enclosing interval [a,b] and a number c in (a,b), if
    !!  f(c)=0 then sets the output a=c. Otherwise determines the new
    !!  enclosing interval: [a,b]=[a,c] or [a,b]=[c,b]. Also updates the
    !!  termination criterion corresponding to the new enclosing interval.

    implicit none

    real(wp),intent(inout)  :: a    !! input as the current left point of the
                                    !! enclosing interval and output as the shrinked
                                    !! new enclosing interval
    real(wp),intent(inout)  :: b    !! input as the current right point of the
                                    !! enclosing interval and output as the shrinked
                                    !! new enclosing interval
    real(wp),intent(inout)  :: c    !! used to determine the new enclosing interval
    real(wp),intent(inout)  :: fa   !! f(a)
    real(wp),intent(inout)  :: fb   !! f(b)
    real(wp),intent(inout)  :: tol  !! input as the current termination
                                    !! criterion and output as the updated termination
                                    !! criterion according to the new enclosing interval
    real(wp),intent(out)    :: d    !! if the new enclosing interval
                                    !! is [a,c] then d=b, otherwise d=a;
    real(wp),intent(out)    :: fd   !! f(d)

    real(wp) :: fc

    ! adjust c if (b-a) is very small or if c is very close to a or b.
    tol = 0.7_wp*tol
    if ((b-a) <= 2.0_wp*tol) then
        c = a+0.5_wp*(b-a)
    else if (c <= a+tol) then
        c = a+tol
    else
        if (c >= b-tol) c = b-tol
    end if

    ! call subroutine to obtain f(c)
    fc = me%f(c)

    ! if c is a root, then set a=c and return. this will terminate the
    ! procedure in the calling routine.
    if (abs(fc) <= me%ftol) then

        a   = c
        fa  = fc
        d   = 0.0_wp
        fd  = 0.0_wp

    else

        ! if c is not a root, then determine the new enclosing interval.
        if ((isign(fa)*isign(fc)) < 0) then
            d   = b
            fd  = fb
            b   = c
            fb  = fc
        else
            d   = a
            fd  = fa
            a   = c
            fa  = fc
        end if

        ! update the termination criterion according to the new enclosing interval.
        if (abs(fb) <= abs(fa)) then
            tol = get_tolerance(b)
        else
            tol = get_tolerance(a)
        end if

    end if

    end subroutine bracket
    !************************************************************************

    !************************************************************************
    pure function isign(x) result(i)

    !! sign of the variable `x` (note: return `0` if `x=0`)

    implicit none

    integer :: i
    real(wp),intent(in) :: x

    if (x > 0.0_wp) then
        i = 1
    else if (x == 0.0_wp) then
        i = 0
    else
        i = -1
    end if

    end function isign
    !************************************************************************

    !************************************************************************
    pure function get_tolerance(b) result(tol)

    !! determines the termination criterion.

    implicit none

    real(wp),intent(in) :: b
    real(wp) :: tol  !! termination criterion: 2*(2*rtol*|b| + atol)

    tol = 2.0_wp * (me%atol + 2.0_wp*abs(b)*me%rtol)

    end function get_tolerance
    !************************************************************************

    !************************************************************************
    pure subroutine newqua(a,b,d,fa,fb,fd,c,k)

    !! uses k newton steps to approximate the zero in (a,b) of the
    !! quadratic polynomial interpolating f(x) at a, b, and d.
    !! safeguard is used to avoid overflow.

    implicit none

    real(wp),intent(in)  :: a
    real(wp),intent(in)  :: b
    real(wp),intent(in)  :: d  !! d lies outside the interval [a,b]
    real(wp),intent(in)  :: fa !! f(a), f(a)f(b)<0
    real(wp),intent(in)  :: fb !! f(b), f(a)f(b)<0
    real(wp),intent(in)  :: fd !! f(d)
    real(wp),intent(out) :: c  !! the approximate zero
                               !! in (a,b) of the quadratic polynomial.
    integer,intent(in)   :: k  !! number of newton steps to take.

    integer  :: ierror,i
    real(wp) :: a0,a1,a2,pc,pdc

    ! initialization
    ! find the coefficients of the quadratic polynomial
    ierror = 0
    a0 = fa
    a1 = (fb-fa)/(b-a)
    a2 = ((fd-fb)/(d-b)-a1)/(d-a)

    do    ! main loop

        ! safeguard to avoid overflow
        if ((a2 == 0.0_wp) .or. (ierror == 1)) then
            c=a-a0/a1
            return
        end if

        ! determine the starting point of newton steps
        if (isign(a2)*isign(fa) > 0) then
            c=a
        else
            c=b
        end if

        ! start the safeguarded newton steps
        do i=1,k
            if (ierror == 0) then
                pc=a0+(a1+a2*(c-b))*(c-a)
                pdc=a1+a2*((2.0_wp*c)-(a+b))
                if (pdc == 0.0_wp) then
                    ierror=1
                else
                    c=c-pc/pdc
                end if
            end if
        end do
        if (ierror/=1) exit

    end do

    end subroutine newqua
    !************************************************************************

    !************************************************************************
    pure function pzero(a,b,d,e,fa,fb,fd,fe) result(c)

    !! uses cubic inverse interpolation of f(x) at a, b, d, and e to
    !! get an approximate root of f(x). this procedure is a slight
    !! modification of aitken-neville algorithm for interpolation
    !! described by stoer and bulirsch in "Intro. to numerical analysis"
    !! springer-verlag. new york (1980).

    implicit none

    real(wp),intent(in) :: a,b,d,e,fa,fb,fd,fe
    real(wp) :: c

    real(wp) :: q11,q21,q31,d21,d31,q22,q32,d32,q33

    q11 = (d-e)*fd/(fe-fd)
    q21 = (b-d)*fb/(fd-fb)
    q31 = (a-b)*fa/(fb-fa)
    d21 = (b-d)*fd/(fd-fb)
    d31 = (a-b)*fb/(fb-fa)

    q22 = (d21-q11)*fb/(fe-fb)
    q32 = (d31-q21)*fa/(fd-fa)
    d32 = (d31-q21)*fd/(fd-fa)
    q33 = (d32-q22)*fa/(fe-fa)

    c = a + q31+q32+q33

    end function pzero
    !************************************************************************

    end subroutine toms748
!*****************************************************************************************

!*****************************************************************************************
!>
!  Zhang's method (with corrections from Stage).
!
!### Reference
!  * A. Zhang, "An Improvement to the Brent's Method",
!    International Journal of Experimental Algorithms (IJEA), Volume (2) : Issue (1) : 2011.
!    https://www.cscjournals.org/download/issuearchive/IJEA/Volume2/IJEA_V2_I1.pdf
!  * S. A. Stage, "Comments on An Improvement to the Brent's Method",
!    International Journal of Experimental Algorithms (IJEA), Volume (4) : Issue (1) : 2013.
!    https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.740.923&rep=rep1&type=pdf

    subroutine zhang(me,ax,bx,fax,fbx,xzero,fzero,iflag)

    implicit none

    class(zhang_solver),intent(inout) :: me
    real(wp),intent(in)  :: ax      !! left endpoint of initial interval
    real(wp),intent(in)  :: bx      !! right endpoint of initial interval
    real(wp),intent(in)  :: fax     !! `f(ax)`
    real(wp),intent(in)  :: fbx     !! `f(ax)`
    real(wp),intent(out) :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out) :: fzero   !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)  :: iflag   !! status flag (`0`=root found, `-2`=max iterations reached)

    real(wp) :: a,b,c,fa,fb,fc,s,fs
    integer :: i !! iteration counter

    iflag = 0
    a  = ax
    b  = bx
    fa = fax
    fb = fbx

    do i = 1, me%maxiter

        c = bisect(a,b)
        fc = me%f(c)
        if (me%solution(c,fc,xzero,fzero)) return

        if (fa/=fc .and. fb/=fc) then
            ! inverse quadratic interpolation
            s = a*fb*fc/((fa-fb)*(fa-fc)) + &
                b*fa*fc/((fb-fa)*(fb-fc)) + &
                c*fa*fb/((fc-fa)*(fc-fb))
            if (a<s .and. s<b) then
                fs = me%f(s)
                if (abs(fs)<=me%ftol) then
                    xzero = s
                    fzero = fs
                    return
                end if
            else
                ! s is not in (a,b)
                s = c ! just use this (there are 3 options in the reference)
                fs = fc
            end if
        else
            ! secant
            if (fa*fc<0.0_wp) then ! root in [a,c]
                s = secant(a,c,fa,fc,ax,bx)
            else ! root in [c,b]
                s = secant(c,b,fc,fb,ax,bx)
            end if
            fs = me%f(s)
            if (me%solution(s,fs,xzero,fzero)) return
        end if

        if (c>s) then
            ! ensures a <= c <= s <= b
            call swap(s,c)
            call swap(fs,fc)
        end if

        if (fc*fs<0.0_wp) then       ! root on [c,s]
            a = c
            b = s
            fa = fc
            fb = fs
        else if (fa*fc<0.0_wp) then  ! root on [a,c]
            b = c
            fb = fc
        else                         ! root on [s,b]
            a = s
            fa = fs
        end if

        if (me%converged(a,b)) exit
        if (i == me%maxiter) iflag = -2 ! max iterations reached

    end do

    ! pick the one closest to the root:
    call choose_best(a,b,fa,fb,xzero,fzero)

    end subroutine zhang
!*****************************************************************************************

!*****************************************************************************************
!>
!  Modified anderson-bjorck-king method. Same as [[anderson_bjorck]], but with
!  an extra initial bisection step.
!
!### See also
!  * Kroger & Torsten, "On-Line Trajectory Generation in Robotic Systems", 2010.
!    https://link.springer.com/content/pdf/bbm%3A978-3-642-05175-3%2F1.pdf

    subroutine anderson_bjorck_king(me,ax,bx,fax,fbx,xzero,fzero,iflag)

    implicit none

    class(anderson_bjorck_king_solver),intent(inout) :: me
    real(wp),intent(in)    :: ax      !! left endpoint of initial interval
    real(wp),intent(in)    :: bx      !! right endpoint of initial interval
    real(wp),intent(in)    :: fax     !! `f(ax)`
    real(wp),intent(in)    :: fbx     !! `f(ax)`
    real(wp),intent(out)   :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out)   :: fzero   !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)    :: iflag   !! status flag (`0`=root found, `-2`=max iterations reached)

    integer :: i !! counter
    logical :: root_found !! convergence in x
    real(wp) :: x1,x2,x3,f1,f2,f3,g,f1tmp

    ! initialize:
    iflag = 0
    x1    = ax
    x2    = bx
    f1    = fax
    f2    = fbx

    ! main loop:
    do i = 1,me%maxiter

        ! bisection step:
        x3 = bisect(x1,x2)

        ! calculate f3:
        f3 = me%f(x3)
        if (me%solution(x3,f3,xzero,fzero)) return

        ! determine a new inclusion interval:
        if (f2*f3<0.0_wp) then
            ! zero lies between x2 and x3
            x1 = x2
            x2 = x3
            f1 = f2
            f2 = f3
        else
            ! zero lies between x1 and x3
            x2 = x3
            f2 = f3
        end if

        ! secant step:
        x3 = secant(x1,x2,f1,f2,ax,bx)

        ! calculate f3:
        f3 = me%f(x3)
        if (me%solution(x3,f3,xzero,fzero)) return

        ! determine a new inclusion interval:
        if (f2*f3<0.0_wp) then
            ! zero lies between x2 and x3
            x1 = x2
            x2 = x3
            f1 = f2
            f2 = f3
            f1tmp = f1
        else
            ! zero lies between x1 and x3
            g = 1.0_wp-f3/f2
            if (g<=0.0_wp) g = 0.5_wp
            x2 = x3
            f1tmp = f1
            f1 = g*f1
            f2 = f3
        end if

        ! check for convergence:
        root_found = me%converged(x1,x2)
        if (root_found .or. i == me%maxiter) then
            call choose_best(x1,x2,f1tmp,f2,xzero,fzero)
            if (.not. root_found) iflag = -2  ! max iterations reached
            exit
        end if

    end do

    end subroutine anderson_bjorck_king
!*****************************************************************************************

!*****************************************************************************************
!>
!  BlendTF blended method of trisection and false position methods.
!
!### Reference
!  * E Badr, S Almotairi, A El Ghamry,
!    "A Comparative Study among New Hybrid Root Finding
!    Algorithms and Traditional Methods", Mathematics 2021, 9, 1306.

    subroutine blendtf(me,ax,bx,fax,fbx,xzero,fzero,iflag)

    implicit none

    class(blendtf_solver),intent(inout) :: me
    real(wp),intent(in)  :: ax      !! left endpoint of initial interval
    real(wp),intent(in)  :: bx      !! right endpoint of initial interval
    real(wp),intent(in)  :: fax     !! `f(ax)`
    real(wp),intent(in)  :: fbx     !! `f(ax)`
    real(wp),intent(out) :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out) :: fzero   !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)  :: iflag   !! status flag (`0`=root found, `-2`=max iterations reached)

    real(wp) :: a1,a2,b1,b2,fa,fb,a,b,xt1,xt2,xf,x,fx,fxt2,&
                fxf,xprev,fxt1,fxprev,fa1,fa2,fb1,fb2
    integer :: i !! iteration counter

    iflag = 0
    a = ax; b = bx
    fa = fax; fb = fbx
    a1 = a; a2 = a
    b1 = b; b2 = b
    fa1 = fa; fa2 = fa
    fb1 = fb; fb2 = fb
    xprev = huge(1.0_wp)
    fxprev = huge(1.0_wp)

    do i = 1, me%maxiter

        if (fa == fb) then
            ! should fallback to bisection if this happens -TODO
            iflag = -3
            return
        end if

        xt1  = (b + 2.0_wp * a) / 3.0_wp
        xt2  = (2.0_wp * b + a) / 3.0_wp
        xf   = a - (fa*(b-a))/(fb-fa)
        x    = xt1
        fxt1 = me%f(xt1); if (me%solution(xt1,fxt1,xzero,fzero)) return
        fxt2 = me%f(xt2); if (me%solution(xt2,fxt2,xzero,fzero)) return
        fxf  = me%f(xf);  if (me%solution(xf,fxf,xzero,fzero))   return
        fx   = fxt1

        if (abs(fxf) < abs(fxt1)) then
            x = xf
            fx = fxf
        elseif (abs(fxt2) < abs(fxt1)) then
            x = xt2
            fx = fxt2
        end if

        ! apply the convergence tols to [a,b]
        if (me%converged(a, b) .or. i == me%maxiter) then
            call choose_best(x,xprev,fx,fxprev,xzero,fzero)
            if (i == me%maxiter) iflag = -2 ! max iterations reached
            exit
        end if
        xprev = x
        fxprev = fx

        if ((fa * fxt1) < 0.0_wp) then
            b1 = xt1; fb1 = fxt1
        else if ((fxt1 * fxt2) < 0.0_wp) then
            a1 = xt1; fa1 = fxt1
            b1 = xt2; fb1 = fxt2
        else
            a1 = xt2; fa1 = fxt2
        end if

        if (fa*fxf < 0.0_wp) then
            b2 = xf; fb2 = fxf
        else
            a2 = xf; fa2 = fxf
        end if

        if (a1>a2) then
            a = a1
            fa = fa1
        else
            a = a2
            fa = fa2
        end if
        if (b1<b2) then
            b = b1
            fb = fb1
        else
            b = b2
            fb = fb2
        end if

    end do

    end subroutine blendtf
!*****************************************************************************************

!*****************************************************************************************
!>
!  Barycentric interpolation method.
!
!### Reference
!  * Enrique Aguilar-Mendez, Roxana Mitzaye Del Castillo,
!    "A highly efficient numerical method to solve non-linear functions using barycentric interpolation"
!    January 2021, Applied Mathematical Sciences 15(7):321-336
!  * Python implemention here: https://github.com/EnriqueAguilarM/

    subroutine barycentric(me,ax,bx,fax,fbx,xzero,fzero,iflag)

        implicit none

        class(barycentric_solver),intent(inout) :: me
        real(wp),intent(in)  :: ax      !! left endpoint of initial interval
        real(wp),intent(in)  :: bx      !! right endpoint of initial interval
        real(wp),intent(in)  :: fax     !! `f(ax)`
        real(wp),intent(in)  :: fbx     !! `f(ax)`
        real(wp),intent(out) :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
        real(wp),intent(out) :: fzero   !! value of `f` at the root (`f(xzero)`)
        integer,intent(out)  :: iflag   !! status flag (`0`=root found, `-2`=max iterations reached)

        integer :: i !! iteration counter
        real(wp) :: k !! a real number in the interval (0, 0.5)
        real(wp) :: x0,x1,x2,x3,x4,x5,x6,&
                    f0,f1,f2,f3,f4,f5,f6,&
                    w0,w1,w2,w3,w22,w32

        real(wp),parameter :: bisection_tol = 1.0e-17_wp  !! if `f0` or `f1` are below this value,
                                                          !! the method switches to bisection.
                                                          !! (should be an input or computed from `epsilon`?)
        real(wp),parameter :: k0 =  1.0_wp / 6.0_wp !! initial value of `k` (should be an input?)

        iflag = 0
        x0 = ax; f0 = fax
        x1 = bx; f1 = fbx
        k = k0

        do i = 1, me%maxiter

            if (i>1) then
                f0 = me%f(x0)
                if (me%solution(x0,f0,xzero,fzero)) return
                f1 = me%f(x1)
                if (me%solution(x1,f1,xzero,fzero)) return
            end if

            if (me%converged(x0, x1) .or. i == me%maxiter) then
                if (i == me%maxiter) iflag = -2 ! max iterations reached
                exit
            end if

            if (abs(f0)<bisection_tol .or. abs(f1)<bisection_tol) then
                ! if f is too small, switch to bisection
                k = 0.5_wp
            end if

            x2 = x0 + k*(x1-x0)   ! first point for interpolation
            f2 = me%f(x2)
            if (me%solution(x2,f2,xzero,fzero)) return

            if (f0*f2 < 0.0_wp) then
                x1 = x2
            else if (f0 == f2) then
                x0 = x2
            else
                x3 = x1 + k*(x0-x1) ! second point for interpolation
                f3 = me%f(x3)
                if (me%solution(x3,f3,xzero,fzero)) return

                if (f3*f1 < 0.0_wp) then
                    x0 = x3
                else if (f3 == f1) then
                    x0 = x2        ! note: typo in python version here?
                    x1 = x3
                else

                    w0 = 1.0_wp/(f0*(f0-f2)*(f0-f3))   ! first quadratic interpolation
                    w2 = 1.0_wp/(f2*(f2-f0)*(f2-f3))
                    w3 = 1.0_wp/(f3*(f3-f0)*(f3-f2))
                    x4 = (w0*x0+w2*x2+w3*x3)/(w0+w3+w2)
                    f4 = me%f(x4)
                    if (me%solution(x4,f4,xzero,fzero)) return

                    w1  = 1.0_wp/(f1*(f1-f2)*(f1-f3))  ! second quadratic interpolation
                    w22 = 1.0_wp/(f2*(f2-f1)*(f2-f3))
                    w32 = 1.0_wp/(f3*(f3-f1)*(f3-f2))
                    x5  = (w1*x1+w22*x2+w32*x3)/(w1+w32+w22)
                    f5  = me%f(x5)
                    if (me%solution(x5,f5,xzero,fzero)) return

                    ! determine if x3,x4 are in the interval [x2,x3]
                    if (x4<x2 .or. x4>x3) then
                        if (x5<x2 .or. x5>x3) then
                            x0 = x2
                            x1 = x3
                        else
                            if (f5*f2<0.0_wp) then
                                x0 = x2
                                x1 = x5
                            else
                                x0 = x5
                                x1 = x3
                            end if
                        end if
                    else
                        if (x5<x2 .or. x5>x3) then
                            if (f4*f2<0.0_wp) then
                                x0 = x2
                                x1 = x4
                            else
                                x0 = x4
                                x1 = x3
                            end if
                        else

                            if (f5*f4>0.0_wp) then
                                if (f4*f2<0.0_wp) then
                                    x0 = x2
                                    x1 = min(x4,x5)
                                else
                                    x0 = max(x4,x5)
                                    x1 = x3
                                end if
                            else
                                ! cubic interpolation if f4 and f5 have opposite signs
                                w0 = w0/(f0-f1)
                                w1 = w1/(f1-f0)
                                w2 = w2/(f2-f1)
                                w3 = w3/(f3-f1)
                                x6 = (w0*x0+w1*x1+w2*x2+w3*x3)/(w0+w1+w3+w2)
                                f6 = me%f(x6)
                                if (me%solution(x6,f6,xzero,fzero)) return

                                ! verify that x6 is in the interval [x2,x3]
                                if (x6<=x2 .or. x6>=x3) then
                                    x0 = min(x4,x5)
                                    x1 = max(x4,x5)
                                else
                                    if (f6*f4<0.0_wp) then
                                        x0 = min(x4,x6)
                                        x1 = max(x4,x6)
                                    else
                                        x0 = min(x5,x6)
                                        x1 = max(x5,x6)
                                    end if
                                end if
                            end if
                        end if
                    end if
                end if
            end if

        end do

        call choose_best(x0,x1,f0,f1,xzero,fzero)

    end subroutine barycentric
!*****************************************************************************************

!*****************************************************************************************
!>
!  For the [[itp]] method, set the optional inputs.

    subroutine itp_optional_inputs(me,k1,k2,n0)

    implicit none

    class(itp_solver),intent(inout) :: me

    real(wp),intent(in),optional :: k1 !! from (0, inf) [Default is 0.1]
    real(wp),intent(in),optional :: k2 !! from [1, 1+phi] [Default is 0.98*(1+phi)]
    integer,intent(in),optional  :: n0 !! [Default is 1

    if (present(k1)) me%k1 = k1
    if (present(k2)) me%k2 = k2
    if (present(n0)) me%n0 = n0

    end subroutine itp_optional_inputs
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the zero of the function f(x) in the interval ax,bx using the
!  Interpolate Truncate and Project (ITP) method.
!
!### See also
!  * Oliveira, I. F. D., Takahashi, R. H. C.,
!    "An Enhancement of the Bisection Method Average Performance Preserving Minmax Optimality",
!    ACM Transactions on Mathematical Software. 47 (1): 5:1-5:24. (2020-12-06)
!
!### Notes
!  * This implementation differs from the reference, in that it has an additional
!    check to make sure the `a,b` values have changed at each iteration. If they
!    haven't, it will perform a basic bisection.

    subroutine itp(me,ax,bx,fax,fbx,xzero,fzero,iflag)

    implicit none

    class(itp_solver),intent(inout) :: me
    real(wp),intent(in)    :: ax      !! left endpoint of initial interval
    real(wp),intent(in)    :: bx      !! right endpoint of initial interval
    real(wp),intent(in)    :: fax     !! `f(ax)`
    real(wp),intent(in)    :: fbx     !! `f(ax)`
    real(wp),intent(out)   :: xzero   !! abscissa approximating a zero of `f` in the interval `ax`,`bx`
    real(wp),intent(out)   :: fzero   !! value of `f` at the root (`f(xzero)`)
    integer,intent(out)    :: iflag   !! status flag (`0`=root found, `-2`=max iterations reached)

    real(wp) :: a,b,ya,yb,x12,r,d,xf,bma,sigma,xt,xitp,yitp,term,aprev,bprev,denom
    integer :: n12, nmax
    integer :: j !! iteration counter
    logical :: root_found !! convergence in x
    logical :: fail !! if we can't interpolate

    real(wp),parameter :: log2 = log(2.0_wp)

    ! initialize:
    if (fax < fbx) then
        a  = ax
        b  = bx
        ya = fax
        yb = fbx
    else
        a  = bx
        b  = ax
        ya = fbx
        yb = fax
    end if
    iflag = 0
    term = (b-a)/(2.0_wp*me%rtol)
    n12 = ceiling ( log(term) / log2 ) ! ceiling(log2(term))
    nmax = n12 + me%n0
    aprev = huge(1.0_wp) ! initialize to unusual values
    bprev = huge(1.0_wp) !

    ! main loop
    do j = 0, min(nmax,me%maxiter)

        x12 = bisect(a,b)
        bma = b - a
        ! note: protect for r<0 as mentioned in paper
        r = max(me%rtol * 2.0_wp ** (nmax-j) - bma/2.0_wp, 0.0_wp)
        d = me%k1 * bma**me%k2

        ! interpolation:
        denom = yb - ya
        fail = abs(denom) <= tiny(1.0_wp)  ! check for divide by zero

        if (.not. fail) then
            xf = (yb*a - ya*b) / denom

            ! truncation:
            sigma = sign(1.0_wp, x12 - xf)
            if (d <= abs(x12-xf)) then
                xt = xf + sigma * d
            else
                xt = x12
            end if

            ! projection:
            if (abs(xt - x12) <= r) then
                xitp = xt
            else
                xitp = x12 - sigma * r
            end if

            ! updating interval:
            yitp = me%f(xitp)
            if (me%solution(xitp,yitp,xzero,fzero)) return
            if (yitp > 0.0_wp) then
                b = xitp
                yb = yitp
            elseif (yitp < 0.0_wp) then
                a = xitp
                ya = yitp
            else
                a = xitp
                b = xitp
            end if

        end if

        if (fail .or. (a==aprev .and. b==bprev)) then
            ! if the interval hasn't changed, then it is stuck.
            ! [this can happen in the test cases].
            ! So just do a bisection
            xitp = bisect(a,b)
            yitp = me%f(xitp)
            if (me%solution(xitp,yitp,xzero,fzero)) return
            if (ya*yitp<0.0_wp) then
                ! root lies between a and xitp
                b = xitp
                yb = yitp
            else
                ! root lies between xitp and b
                a = xitp
                ya = yitp
            end if
        end if
        aprev = a
        bprev = b

        ! check for convergence:
        root_found = me%converged(a,b)
        if (root_found .or. j==me%maxiter) then
            call choose_best(a,b,ya,yb,xzero,fzero)
            if (.not. root_found) iflag = -2  ! max iterations reached
            exit
        end if

    end do

    end subroutine itp
!*****************************************************************************************

!*****************************************************************************************
!>
!  Determines convergence in x based on if the reltol or abstol is satisfied.

    function converged(me,a,b)

    implicit none

    class(root_solver),intent(inout) :: me
    real(wp),intent(in) :: a !! old value
    real(wp),intent(in) :: b !! new value
    logical :: converged

    real(wp) :: d

    ! original way:
    ! converged = (abs(b-a) <= abs(b)*me%rtol + me%atol) exit

    d = abs(b-a)

    if (d <= me%atol) then
        ! absolute
        converged = .true.
    else
        ! relative
        if (a /= 0.0_wp) then
            converged = d / abs(a) <= me%rtol
        else
            converged = .false.
        end if
    end if

    end function converged
!*****************************************************************************************

!*****************************************************************************************
!>
!  Given two points with two function evaluations, choose the best one
!  (the one closest to the root).

    pure subroutine choose_best(x1,x2,f1,f2,xbest,fbest)

    implicit none

    real(wp),intent(in) :: x1,x2
    real(wp),intent(in) :: f1,f2
    real(wp),intent(out) :: xbest
    real(wp),intent(out) :: fbest

    if (abs(f1)<abs(f2)) then
        xbest = x1
        fbest = f1
    else
        xbest = x2
        fbest = f2
    end if

    end subroutine choose_best
!*****************************************************************************************

!*****************************************************************************************
!>
!  Bisection step.

    pure function bisect(x1,x2) result(x3)

    implicit none

    real(wp),intent(in) :: x1,x2
    real(wp) :: x3 !! point half way between x1 and x2

    x3 = (x1 + x2) / 2.0_wp

    end function bisect
!*****************************************************************************************

!*****************************************************************************************
!>
!  Regula Falsi step.
!  With a protection to fall back to bisection if:
!
!   * the computed point is outside the original interval ([ax,bx]).
!   * f2 == f1

    function regula_falsi_step(x1,x2,f1,f2,ax,bx) result(x3)

    implicit none

    real(wp),intent(in) :: x1,x2,f1,f2
    real(wp),intent(in) :: ax !! original interval lower bound
    real(wp),intent(in) :: bx !! original interval upper bound
    real(wp) :: x3 !! intersection of line connecting x1,x2 with x-axis

    real(wp) :: delta

    delta = f2-f1

    if (delta /= 0.0_wp) then
        ! intersection with x-axis of line connecting the two points:
        x3 = x1 - (f1/delta) * (x2-x1)
        if (x3>ax .and. x3<bx) return ! must be a new point in the range
    end if

    ! fall back to bisection for any problem
    x3 = bisect(x1,x2)

    end function regula_falsi_step
!*****************************************************************************************

!*****************************************************************************************
!>
!  Secent step.
!  With a protection to fall back to bisection if:
!
!   * the computed point is outside the original interval ([ax,bx]).
!   * f2 == f1

    pure function secant(x1,x2,f1,f2,ax,bx) result(x3)

    implicit none

    real(wp),intent(in) :: x1,x2,f1,f2
    real(wp),intent(in) :: ax !! original interval lower bound
    real(wp),intent(in) :: bx !! original interval upper bound
    real(wp) :: x3 !! intersection of secant step with x-axis

    if (f2==f1) then
        x3 = bisect(x1,x2)
    else
        ! secant step:
        x3 = x2 - f2 / ( (f2 - f1) / (x2 - x1) )
        if (x3<ax .or. x3>bx) x3 = bisect(x1,x2)
    end if

    end function secant
!*****************************************************************************************

!*****************************************************************************************
!>
!  Swap two real(wp) values.

    pure elemental subroutine swap(a,b)

    implicit none

    real(wp),intent(inout) :: a
    real(wp),intent(inout) :: b

    real(wp) :: tmp

    tmp = a
    a   = b
    b   = tmp

    end subroutine swap
!*****************************************************************************************

!*****************************************************************************************
!>
!  lowercase a string.

    pure function lowercase(str) result(s_lower)

    implicit none

    character(len=*),intent(in) :: str      !! input string
    character(len=(len(str)))   :: s_lower  !! lowercase version of the string

    integer :: i  !! counter
    integer :: j  !! index of uppercase character

    character(len=*),parameter :: upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' !! uppercase characters
    character(len=*),parameter :: lower = 'abcdefghijklmnopqrstuvwxyz' !! lowercase characters

    s_lower = str

    do i = 1, len_trim(str)
        j = index(upper,s_lower(i:i))
        if (j>0) s_lower(i:i) = lower(j:j)
    end do

    end function lowercase
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns true if this is a solution and sets `xzero` and `fzero`.

    logical function solution(me,x,f,xzero,fzero)

    implicit none

    class(root_solver),intent(inout) :: me
    real(wp),intent(in) :: x
    real(wp),intent(in) :: f
    real(wp),intent(inout) :: xzero
    real(wp),intent(inout) :: fzero

    if (abs(f) <= me%ftol) then
        xzero = x
        fzero = f
        solution = .true.
    else
        solution = .false.
    end if

    end function solution
!*****************************************************************************************

!*****************************************************************************************
    end module root_module
!*****************************************************************************************

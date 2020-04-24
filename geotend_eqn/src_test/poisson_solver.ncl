; Poisson equation 3d solver /*{{{*/
undef("poisson_solver")
function poisson_solver(F[*][*][*]:numeric, opt:logical)
;
;  - Input argument
;    F : Forcing [lev] x [lat] x [lon]
;    opt : Option (See below)
;
;  - Description : Compute poisson's equation 
;    Z : [ ( del2 + (d/dz)(f0^2/N^2)*(d/dz) )S = F ]
;    P : [ ( del2 + (d/dp)(f0^2/sig)*(d/dp) )S = F ]
;    del2 : laplacian 2d,   f0 : Corioli at 45N
;    S    : Solution,       N  : Bouyancy frequency
;    F    : Forcing,        sig: static stability in p-coordinate
;                               = -(RT0/p)*(dlnTHE0/dp)
;
;  - Option
;    opt@T0    = T(nprs) or T(nlev)
;       @sig   = sig <- constant value
;       @lev   = vertical coordinate       (default: level from F )
;       @lev_units = "hPa", "Pa", "m"
;       @dt    = time [s] (emprical value) (default: 1.e9  [s]  )
;       @thres = threshold for iteration   (default: 1.e-4      )
;       @f0    = Coriolis parameter        (default: 1.e-4 [s-1])
;       @bound_type  = "N"(neumann),       (default: N)
;                      "D"(dirichlet)
;       @upper_bound = Upper boundary      (default: 0)
;       @lower_bound = Lower boundary      (default: 0)
;       @radius_earth = Earth radius       (default: 6371220 [m])
;
;  - Note! BC unit must be divieded [Pa] when neumann boundary type used.
;  - Laplacian part solve using spherical harmonics.
;   (See google drive 
;    /project_doc/TODO_Changma/etc/Spheical Harmonics.pdf)
;  - Vertical diffusion part is solved using Relaxaion method.
;
;  - Update
;  - 2020.01.29 : ERROR
;                If call `get_sphere_radius()` in a function, 
;                it make segmentation fault.
;                e.g) R = tofloat(get_sphere_radius())
;               
begin
    ;-- Set variables
    ; Define dt instead of relaxation coefficient for convenience
    dt    = 1.e9  ; delta time                 [s]
    thres = 1.e-4 ; threshold for iteration
    f0    = 1.e-4 ; Coriolis parameter         [/s]
    g     = 9.8   ; gravitational acceleration [m/s2]
    R     = 287.  ; Gas constant               [J/kg/K]
    Cp    = 1005. ; Specific heat capacity     [J/kg/K]
    ; a     = get_sphere_radius()    ; <- update (2020.01.29)
    a     = 6371220.   ; Earth radius [m]
    OMG   = 7.2921e-5
    pi    = 4*atan(1.)
    d2r   = pi/180.

    ; TODO: Add if loop for no dimension info
    ;-- Get coordinate
    name_lev = F!0
    name_lat = F!1
    name_lon = F!2
    lev = tofloat(F&$name_lev$) ; aviod error Assignment type mismatch
    lat = F&$name_lat$
    lon = F&$name_lon$
    nlev = dimsizes(lev)
    nlat = dimsizes(lat)
    nlon = dimsizes(lon)
    f    = 2*OMG*sin(lat*d2r) ; Coriolis parameter         [/s]
    f2   = conform_dims((/nlat, nlon/), f*f, 0)

    top = new((/nlat, nlon/), typeof(F))
    top = 0.
    bot = top
    bc_type  = "N" ; boundary condition type

    if (isatt(opt,    "dt")) then dt    = opt@dt    end if    
    if (isatt(opt, "thres")) then thres = opt@thres end if    
    if (isatt(opt,    "f0")) then f0    = opt@f0    end if    
    if (isatt(opt,   "lev")) then lev   = opt@lev   end if    
    if (isatt(opt, "upper_bound")) then top = opt@upper_bound end if
    if (isatt(opt, "lower_bound")) then bot = opt@lower_bound end if
    if (isatt(opt, "bound_type")) then bc_type = opt@bound_type end if 

    ; TODO : how to resolve corioli parameter?
    ;   (f_00*Y_00 + f_20*Y_20) * a_lm*Y_lm = c_lm*Y_lm
    ;   Y_20 * Y_lm = Y_(l+?)m
    ; method 1 : Use shseC in while loop
    ; method 2 : Found Y_(l+?)m mathmatically
    ; For the case of `f = 2*OMG*sin(lat)`
    ; f2_lm = shaeC(f2)
    ; f2_00 = f2_lm(0,0,0)
    ; f2_20 = f2_lm(0,2,0)

    ; Calculate N2(height) or Sigma(pressure)
    T2B = False     ; True: Top-to-Bottom
    if (isatt(opt, "T0")) then
        T0 = opt@T0
        if (opt@lev_units .eq. "m") then ; <= height coordinate
            if (lev(0) .gt. lev(1)) then T2B = True end if    
            dTdz = center_finite_diff(T0, lev, False, 0)
            sig = (g/T0)*( dTdz + (g/Cp) )  ; = N2 (ref. MAD Andrew 122p)
        else ; <= pressure coordinate
            if (lev(1) .gt. lev(0)) then T2B = True end if    
            if (opt@lev_units .eq. "hPa") then lev := lev * 100. end if 
            sig = static_stability(lev, T0, 0, 0)   ; = Sp
            sig := (R/lev)*sig  ; σ ≡ −(RT/p)(dlnθ/dp) (ref. Holton 4th 150p)
        end if  
        print("===== Level  |  Temeprature  |  Sigma =====")
        print(sprintf("%7.0f",lev)+" [Pa] |   "+ \
              sprintf("%6.2f", T0)+" [K]  | "  + sig +" [m2 Pa-2 s-2]")
        print("===========================================")
    else    ; TODO: classify "m" and "Pa"
        sig = conform(lev, opt@sig, 0)
        print("===== Level  |  Sigma =====")
        print(sprintf("%7.0f",lev)+" [Pa] | "+ sig +" [m2 Pa-2 s-2]")
        print("===========================================")
    end if    
    

    ;-- Solver Part
    dlev =   lev(1:nlev-1) - lev(0:nlev-2)         ; =dp or dz
    msig = ( sig(1:nlev-1) + sig(0:nlev-2) ) / 2.  ; index = k + 0.5

    b_lm = shaeC(F)            ; [2] x [lev] x [lat] x [lon]
    a_lm = todouble(b_lm*0.)   ; Because R is double

    bc_top = shaeC(top)
    bc_bot = shaeC(bot)
    
    n    = ispan(0,dimsizes(lat)-1,1)
    nn1  = conform(a_lm(:,0,:,:),n*(n+1),1)  ; [2] x [l] x [m]

    err = 1.d30
    begTime = get_cpu_time()
    do loop = 1, 10000
        ; Boundary value
        if (bc_type .eq. "N") then ; neumann condition
            if (T2B) then    ; True: Top-to-Bottom
                a_lm(:,     0,:,:) = a_lm(:,1,:,:)      \
                                   + bc_top*dlev(0)
                a_lm(:,nlev-1,:,:) = a_lm(:,nlev-2,:,:) \
                                   - bc_bot*dlev(nlev-2)
            else
                a_lm(:,nlev-1,:,:) = a_lm(:,nlev-2,:,:) \
                                   + bc_top*dlev(nlev-2)
                a_lm(:,     0,:,:) = a_lm(:,1,:,:)      \
                                   - bc_bot*dlev(0)
            end if    
        else                       ; dirichlet condition
            if (T2B) then    ; True: Top-to-Bottom
                a_lm(:,     0,:,:) = bc_top
                a_lm(:,nlev-1,:,:) = bc_bot
            else
                a_lm(:,nlev-1,:,:) = bc_top
                a_lm(:,     0,:,:) = bc_bot
            end if    
        end if    
    
        dummy = a_lm
        do k = 1, dimsizes(lev)-2   ; omit boundary
            relax1 = dt/(sig(k  )*dlev(k  )*0.5*(lev(k+1)-lev(k-1)))
            relax2 = dt/(sig(k-1)*dlev(k-1)*0.5*(lev(k+1)-lev(k-1)))

            ; laps: Horizontal laplacian
            laps = -(nn1/a^2)*dt*a_lm(:,k,:,:)
            laps(:,0,0) = 0.        ; No mean value for laplacian

            ; c_lm: Vertical laplacian
            c_lm = - (relax1+relax2)*a_lm(:,k,:,:)        \
                   + relax1*a_lm(:,k+1,:,:)               \
                   + relax2*a_lm(:,k-1,:,:)

            if (dimsizes(f0) .eq. 1) then ; f = f0 (constant)
                c_lm = (f0^2)*c_lm
            else                     ; f = f(lat)
            ; TODO : update using convolution sum instead of 
            ;       transfom method for speed.
                cphy = shseC(c_lm, nlon)
                cphy = cphy*f2
                c_lm = shaeC(cphy)
            end if    

            a_lm(:,k,:,:) = a_lm(:,k,:,:)      \
                            + laps             \  ; Horizontal laplacian
                            + c_lm             \  ; Vertical laplacian
                            - dt*b_lm(:,k,:,:)    ; Forcing
        end do    
        err = avg(abs(a_lm - dummy))
        print(err+"")

        ; TODO: Stop criteria
        ; 1. Absolute error stop < err(1) / 1000.
        ; 2. derr 
        rel_err = err / avg(abs(a_lm))
        if (rel_err .lt. thres) then break end if    
        ; if (rel_err .lt. 1.e-4) then break end if    
        if (loop .eq. 10000) then
            print("Note! Number of iteration over maximum")
        end if    
    end do 
   
    solution = shseC(a_lm, nlon)
    copy_VarCoords(F, solution)
    print("Relaxation time: " + (get_cpu_time() - begTime)+" [s]")

    return solution
end ; /*}}}*/


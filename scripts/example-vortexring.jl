#=##############################################################################
# DESCRIPTION
    Vortex ring simulation.

    This module contains multiple functions that have been copied/pasted
    from FLOWVPM by the author. See
    https://github.com/byuflowlab/FLOWVPM.jl/tree/master/examples/vortexrings

# AUTHORSHIP
  * Author    : Eduardo J. Alvarez
  * Email     : Edo.AlvarezR@gmail.com
  * Created   : Nov 9th, 2024
=###############################################################################


import Elliptic
import Roots
import Cubature
import LinearAlgebra: I, norm, cross
import Printf: @printf

modulepath = splitdir(@__FILE__)[1]         # Path to this module

# Load FastMultipole `vortex.jl` module
using StaticArrays
include(joinpath(modulepath, "..", "test", "vortex.jl"))

# Check whether this is called within a unit test or not
try
    istest
catch
    global istest = false
end



# -------------- USEFUL FUNCTIONS ----------------------------------------------
"Number of particles used to discretized a ring"
number_particles(Nphi, nc; extra_nc=0) = Int( Nphi * ( 1 + 8*sum(1:(nc+extra_nc)) ) )

"Intervals of each ring"
function calc_ring_invervals(nrings, Nphis, ncs, extra_ncs)
    intervals = [0]
    for ri in 1:nrings
        Np = number_particles(Nphis[ri], ncs[ri]; extra_nc=extra_ncs[ri])
        push!(intervals, intervals[end] + Np)
    end
    return intervals
end

"Analytic self-induced velocity of an inviscid ring"
Uring(circulation, R, Rcross, beta) = circulation/(4*pi*R) * ( log(8*R/Rcross) - beta )

"""
  `addvortexring(pfield, circulation, R, AR, Rcross, Nphi, nc, sigma;
extra_nc=0, O=zeros(3), Oaxis=eye(3))`

Adds a vortex ring to the particle field `pfield`. The ring is discretized as
described in Winckelmans' 1989 doctoral thesis (Topics in Vortex Methods...),
where the ring is an ellipse of equivalent radius `R=sqrt(a*b)`, aspect ratio
`AR=a/b`, and cross-sectional radius `Rcross`, where `a` and `b` are the
semi-major and semi-minor axes, respectively. Hence, `AR=1` defines a circle of
radius `R`.

The ring is discretized with `Nphi` cross section evenly spaced and the
thickness of the toroid is discretized with `nc` layers, using particles with
smoothing radius `sigma`. Here, `nc=0` means that the ring is represented only
with particles centered along the centerline, and `nc>0` is the number of layers
around the centerline extending out from 0 to `Rcross`.

Additional layers of empty particles (particle with no strength) beyond `Rcross`
can be added with the optional argument `extra_nc`.

The ring is placed in space at the position `O` and orientation `Oaxis`,
where `Oaxis[:, 1]` is the major axis, `Oaxis[:, 2]` is the minor axis, and
`Oaxis[:, 3]` is the line of symmetry.
"""
function vortexring(circulation::Real,
                            R::Real, AR::Real, Rcross::Real,
                            Nphi::Int, nc::Int, sigma::Real; extra_nc::Int=0,
                            O::Vector{<:Real}=zeros(3), Oaxis=I,
                            verbose=true, v_lvl=0
                            )

    maxnp = number_particles(Nphi, nc; extra_nc=extra_nc)
    pfield = zeros(7, maxnp)

    # ERROR CASE
    if AR < 1
        error("Invalid aspect ratio AR < 1 (AR = $(AR))")
    end

    a = R*sqrt(AR)                      # Semi-major axis
    b = R/sqrt(AR)                      # Semi-minor axis

    fun_S(phi, a, b) = a * Elliptic.E(phi, 1-(b/a)^2) # Arc length from 0 to a given angle
    Stot = fun_S(2*pi, a, b)            # Total perimeter length of centerline

                                        # Non-dimensional arc length from 0 to a given value <=1
    fun_s(phi, a, b) = fun_S(phi, a, b)/fun_S(2*pi, a, b)
                                        # Angle associated to a given non-dimensional arc length
    fun_phi(s, a, b) = abs(s) <= eps() ? 0 :
                     abs(s-1) <= eps() ? 2*pi :
                     Roots.fzero( phi -> fun_s(phi, a, b) - s, (0, 2*pi-eps()), atol=1e-16, rtol=1e-16)

                                        # Length of a given filament in a
                                        # cross section cell
    function fun_length(r, tht, a, b, phi1, phi2)
        S1 = fun_S(phi1, a + r*cos(tht), b + r*cos(tht))
        S2 = fun_S(phi2, a + r*cos(tht), b + r*cos(tht))

        return S2-S1
    end
                                        # Function to be integrated to calculate
                                        # each cell's volume
    function fun_dvol(r, args...)
        return r * fun_length(r, args...)
    end
                                        # Integrate cell volume
    function fun_vol(dvol_wrap, r1, tht1, r2, tht2)
        (val, err) = Cubature.hcubature(dvol_wrap,  [r1, tht1], [r2, tht2];
                                           reltol=1e-8, abstol=0, maxevals=1000)
        return val
    end

    invOaxis = inv(Oaxis)               # Add particles in the global coordinate system
    np = 0
    function addparticle(pfield, X, Gamma, sigma, vol, circulation)
        X_global = Oaxis*X + O
        Gamma_global = Oaxis*Gamma

        np += 1
        pfield[1:3, np] .= X_global
        pfield[4, np] = sigma
        pfield[5:7, np] .= Gamma_global

    end

    rl = Rcross/(2*nc + 1)              # Radial spacing between cross-section layers
    dS = Stot/Nphi                      # Perimeter spacing between cross sections
    ds = dS/Stot                        # Non-dimensional perimeter spacing

    omega = circulation / (pi*Rcross^2) # Average vorticity

    org_np = 0

    # Discretization of torus into cross sections
    for N in 0:Nphi-1

        # Non-dimensional arc-length position of cross section along centerline
        sc1 = ds*N                      # Lower bound
        sc2 = ds*(N+1)                  # Upper bound
        sc = (sc1 + sc2)/2              # Center

        # Angle of cross section along centerline
        phi1 = fun_phi(sc1, a, b)       # Lower bound
        phi2 = fun_phi(sc2, a, b)       # Upper bound
        phic = fun_phi(sc, a, b)        # Center

        Xc = [a*sin(phic), b*cos(phic), 0]  # Center of the cross section
        T = [a*cos(phic), -b*sin(phic), 0]  # Unitary tangent of this cross section
        T ./= norm(T)
        T .*= -1                            # Flip to make +circulation travel +Z
                                        # Local coordinate system of section
        Naxis = hcat(T, cross([0,0,1], T), [0,0,1])

                                        # Volume of each cell in the cross section
        dvol_wrap(X) = fun_dvol(X[1], X[2], a, b, phi1, phi2)


        # Discretization of cross section into layers
        for n in 0:nc+extra_nc

            if n==0           # Particle in the center

                r1, r2 = 0, rl              # Lower and upper radius
                tht1, tht2 = 0, 2*pi        # Left and right angle
                vol = fun_vol(dvol_wrap, r1, tht1, r2, tht2) # Volume
                X = Xc                      # Position
                Gamma = omega*vol*T         # Vortex strength
                                            # Filament length
                length = fun_length(0, 0, a, b, phi1, phi2)
                                            # Circulation
                crcltn = norm(Gamma) / length

                addparticle(pfield, X, Gamma, sigma, vol, crcltn)

            else              # Layers

                rc = (1 + 12*n^2)/(6*n)*rl  # Center radius
                r1 = (2*n-1)*rl             # Lower radius
                r2 = (2*n+1)*rl             # Upper radius
                ncells = 8*n                # Number of cells
                deltatheta = 2*pi/ncells    # Angle of cells

                # Discretize layer into cells around the circumference
                for j in 0:(ncells-1)

                    tht1 = deltatheta*j     # Left angle
                    tht2 = deltatheta*(j+1) # Right angle
                    thtc = (tht1 + tht2)/2  # Center angle
                    vol = fun_vol(dvol_wrap, r1, tht1, r2, tht2) # Volume
                                            # Position
                    X = Xc + Naxis*[0, rc*cos(thtc), rc*sin(thtc)]

                                            # Vortex strength
                    if n<=nc    # Ring particles
                        Gamma = omega*vol*T
                    else        # Particles for viscous diffusion
                        Gamma = eps()*T
                    end
                                            # Filament length
                    length = fun_length(rc, thtc, a, b, phi1, phi2)
                                            # Circulation
                    crcltn = norm(Gamma) / length


                    addparticle(pfield, X, Gamma, sigma, vol, crcltn)
                end

            end

        end

    end

    if verbose
        println("\t"^(v_lvl)*"Number of particles: $(np - org_np)")
    end

    return pfield
end


"""
Calculate centroid, radius, and cross-section radius of all rings from the
position of particles weighted by vortex strength.
"""
function calc_rings_weighted!(outZ, outR, outsgm, vortexparticles, nrings, intervals)

    # Iterate over each ring
    for ri in 1:nrings

        # Calculate centroid
        outZ[ri] .= 0
        magGammatot = 0
        for pi in (intervals[ri]+1):(intervals[ri+1])

            P = vortexparticles.bodies[pi]
            normGamma = norm(P.strength)
            magGammatot += normGamma

            for i in 1:3
                outZ[ri][i] += normGamma*P.position[i]
            end

        end
        outZ[ri] ./= magGammatot

        # Calculate ring radius and cross-section radius
        outR[ri], outsgm[ri] = 0, 0
        for pi in (intervals[ri]+1):(intervals[ri+1])

            P = vortexparticles.bodies[pi]
            normGamma = norm(P.strength)

            outR[ri] += normGamma*sqrt((P.position[1] - outZ[ri][1])^2 + (P.position[2] - outZ[ri][2])^2 + (P.position[3] - outZ[ri][3])^2)
            outsgm[ri] += normGamma*P.sigma

        end
        outR[ri] /= magGammatot
        outsgm[ri] /= magGammatot

    end

    return nothing
end








# -------------- SIMULATION PARAMETERS -----------------------------------------
nsteps    = istest ? 1 : 1000         # Number of time steps
Rtot      = istest ? 0.00001 : 2.0        # (m) run simulation for equivalent
                                        #     time to this many radii

nrings    = 1                           # Number of rings
nc        = 1                           # Number of toroidal layers of discretization
# nc        = 0
dZ        = 0.1                         # (m) spacing between rings
circulations = 1.0*ones(nrings)         # (m^2/s) circulation of each ring
Rs        = 1.0*ones(nrings)            # (m) radius of each ring
ARs       = 1.0*ones(nrings)            # Aspect ratio AR = a/r of each ring
Rcrosss   = 0.15*Rs                     # (m) cross-sectional radii
sigmas    = Rcrosss                     # Particle smoothing of each radius
Nphis     = 100*ones(Int, nrings)       # Number of cross sections per ring
ncs       = nc*ones(Int, nrings)        # Number layers per cross section
extra_ncs = 0*ones(Int, nrings)         # Number of extra layers per cross section
Os        = [[0, 0, dZ*(ri-1)] for ri in 1:nrings]  # Position of each ring
Oaxiss    = [I for ri in 1:nrings]      # Orientation of each ring
nref      = 1                           # Reference ring

beta      = 0.5                         # Parameter for theoretical velocity
faux      = 0.25                        # Shrinks the discretized core by this factor

pfields = []

for ringi in 1:nrings
    pfield = vortexring(circulations[ringi],
                                Rs[ringi], ARs[ringi], Rcrosss[ringi],
                                Nphis[ringi], ncs[ringi], sigmas[ringi]; extra_nc=extra_ncs[ringi],
                                O=Os[ringi], Oaxis=Oaxiss[ringi],
                                verbose=true, v_lvl=0
                                )
    push!(pfields, pfield)
end

pfield = hcat(pfields...)

vortexparticles = VortexParticles(pfield);

# save_vtk("vortexring-000", vortexparticles)



# -------------- RUN SIMULATION ------------------------------------------------
Uref = Uring(circulations[nref], Rs[nref], Rcrosss[nref], beta) # (m/s) reference velocity
dt = (Rtot/Uref) / nsteps         # (s) time step

# Create folder to save simulation
savepath = "example-vortexring"

if !istest
    if isdir(savepath)
        rm(savepath, recursive=true)
    end
    mkpath(savepath)
end

# Run simulation
convect!(vortexparticles, nsteps;
        # integration options
        integrate=Euler(dt),
        # fmm options
        fmm_p=4, fmm_ncrit=50, fmm_multipole_threshold=0.4,
        # direct=false,
        direct=true,
        # save options
        save=!istest, filename=joinpath(savepath, "vortexring"), compress=false,
    )


# --------------- COMPARE TO ANALYTIC SOLUTION ---------------------------------

# Calculate centroid of ring by end of simulation
tend = dt*nsteps                              # (s) simulation end time
Z_vpm = [zeros(3) for ri in 1:nrings]         # Centroid position
R_vpm, sgm_vpm = zeros(nrings), zeros(nrings) # Ring and cross section radii
intervals = calc_ring_invervals(nrings, Nphis, ncs, extra_ncs)

calc_rings_weighted!(Z_vpm, R_vpm, sgm_vpm, vortexparticles, nrings, intervals)

# Calculate resulting ring velocity
U_vpm = norm(Z_vpm[1] - Os[1]) / tend

# Calculate analytic ring velocity
U_ana = Uring(circulations[nref], Rs[nref], Rcrosss[nref], beta)

# Error
err = (U_vpm - U_ana) / U_ana

@printf "%sVortex ring self-induced velocity verification\n"    "\n"*"\t"^1
@printf "%sAnalytical velocity:\t\t%1.3f m/s\n"                 "\t"^2 U_ana
@printf "%sResulting velocity:\t\t%1.3f m/s\n"                  "\t"^2 U_vpm
@printf "%sError:\t\t\t\t%1.8f﹪\n"                              "\t"^2 err*100

# Test result
testresult = abs(err) < 0.01

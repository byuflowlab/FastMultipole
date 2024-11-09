using FLOWMath
using FastMultipole
using BenchmarkTools
using LegendrePolynomials
using PythonPlot

function theta_integral(ρ_max, ΔC, n, ρ̂)
    stuff = (-ρ_max^2 + ΔC^2 + ρ̂^2) / (2 * ΔC * ρ̂)
    θs = range(0, stop=acos(stuff), length=100)
    y2 = zeros(length(θs))
    for i in eachindex(θs)
        y2[i] = abs(Plm(cos(θs[i]),n,0)) * sin(θs[i])
    end
    return trapz(θs, y2)
end

function manual_integral(ρ_max, ΔC, n)
    ρ̂s = range(ΔC-ρ_max, stop=ΔC+ρ_max, length=100)
    y1 = zeros(length(ρ̂s))
    for i in eachindex(y1)
        ρ̂ = ρ̂s[i]
        y1[i] = theta_integral(ρ_max, ΔC, n, ρ̂) / ρ̂^(n-1)
    end
    return trapz(ρ̂s, y1)
end

function manual_sum(r, ρ_max, A, p, ΔC)
    println("Manual Sum: p=$p")
    val = 0.0
    ρ_max_3 = ρ_max^3
    for n in p+1:100
        val += 3*A*r^n/(2*ρ_max_3) * manual_integral(ρ_max, ΔC, n)
    end
    return val
end

function theta_integral2(ρ_max, ΔC, n, ρ̂)
    stuff = (-ρ_max^2 + ΔC^2 + ρ̂^2) / (2 * ΔC * ρ̂)
    θs = range(0, stop=acos(stuff), length=100)
    y2 = zeros(length(θs))
    for i in eachindex(θs)
        y2[i] = sin(θs[i])
    end
    return trapz(θs, y2)
end

function manual_integral2(ρ_max, ΔC, n)
    ρ̂s = range(ΔC-ρ_max, stop=ΔC+ρ_max, length=100)
    y1 = zeros(length(ρ̂s))
    for i in eachindex(y1)
        ρ̂ = ρ̂s[i]
        y1[i] = theta_integral2(ρ_max, ΔC, n, ρ̂) / ρ̂^(n-1)
    end
    return trapz(ρ̂s, y1)
end

function manual_sum2(r, ρ_max, A, p, ΔC)
    println("Manual Sum 2: p=$p")
    val = 0.0
    ρ_max_3 = ρ_max^3
    for n in p+1:100
        val += 3*A*r^n/(2*ρ_max_3) * manual_integral2(ρ_max, ΔC, n)
    end
    return val
end

function manual_integral3(ρ_max, ΔC, n)
    ρ̂s = range(ΔC-ρ_max, stop=ΔC+ρ_max, length=100)
    y1 = zeros(length(ρ̂s))
    for i in eachindex(y1)
        ρ̂ = ρ̂s[i]
        y1[i] = (2*ΔC*ρ̂ + ρ_max^2 - ΔC^2 - ρ̂^2) / (2*ΔC*ρ̂^n)
    end
    return trapz(ρ̂s, y1)
end

function manual_sum3(r, ρ_max, A, p, ΔC)
    println("Manual Sum 3: p=$p")
    val = 0.0
    ρ_max_3 = ρ_max^3
    for n in p+1:100
        val += 3*A*r^n/(2*ρ_max_3) * manual_integral3(ρ_max, ΔC, n)
    end
    return val
end

function manual_integral4(ρ_max, ΔC, n)
    val = 0.0
    x = ΔC + ρ_max
    val += x^(2-n) / (2-n) + (ρ_max^2 - ΔC^2) * x^(1-n) / (1-n) / (2*ΔC) - x^(3-n) / (2*ΔC*(3-n))
    x = ΔC - ρ_max
    val -= x^(2-n) / (2-n) + (ρ_max^2 - ΔC^2) * x^(1-n) / (1-n) / (2*ΔC) - x^(3-n) / (2*ΔC*(3-n))
    return val
end

function manual_sum4(r, ρ_max, A, p, ΔC)
    println("Manual Sum 4: p=$p")
    val = 0.0
    ρ_max_3 = ρ_max^3
    for n in p+1:100
        val += 3*A*r^n/(2*ρ_max_3) * manual_integral4(ρ_max, ΔC, n)
    end
    return val
end

function η(p, ΔC, ρ_max, r)
    val = 0.0
    p < 1 && (val += log((ΔC+ρ_max)/(ΔC-ρ_max)))
    p < 0 && (val += 2*ρ_max/r)
    p < -1 && (val += 2*ΔC*ρ_max/(r*r))
    return val
end

function manual_sum5(r, ρ_max, A, p, ΔC)
    println("Manual Sum 5: p=$p")
    logterm = log((1-r/(ΔC+ρ_max)) / (1-r/(ΔC-ρ_max)))
    t1 = 0.0
    for n in 1:p-2
        t1 += ((r/(ΔC+ρ_max))^n - (r/(ΔC-ρ_max))^n) / n
    end
    t2 = 0.0
    for n in 1:p-1
        t2 += ((r/(ΔC+ρ_max))^n - (r/(ΔC-ρ_max))^n) / n
    end
    t3 = 0.0
    for n in 1:p-3
        t3 += ((r/(ΔC+ρ_max))^n - (r/(ΔC-ρ_max))^n) / n
    end

    val = 0.0
    val += r * (η(p-1, ΔC, ρ_max, r) + logterm + t1)
    val += (ρ_max^2-ΔC^2)/(2*ΔC) * (η(p, ΔC, ρ_max, r) + logterm + t2)
    val += -r^2 / (2*ΔC) * (η(p-2, ΔC, ρ_max, r) + logterm + t3)
    val *= 3*A*r/(2*ρ_max^3)

    return val
end

function manual_sum5_rel(r, ρ_max, A, p, ΔC)
    println("Manual Sum 5 rel: p=$p")
    logterm = log((1-r/(ΔC+ρ_max)) / (1-r/(ΔC-ρ_max)))
    t1 = 0.0
    for n in 1:p-2
        t1 += ((r/(ΔC+ρ_max))^n - (r/(ΔC-ρ_max))^n) / n
    end
    t2 = 0.0
    for n in 1:p-1
        t2 += ((r/(ΔC+ρ_max))^n - (r/(ΔC-ρ_max))^n) / n
    end
    t3 = 0.0
    for n in 1:p-3
        t3 += ((r/(ΔC+ρ_max))^n - (r/(ΔC-ρ_max))^n) / n
    end

    val = 0.0
    val += r * (η(p-1, ΔC, ρ_max, r) + logterm + t1)
    val += (ρ_max^2-ΔC^2)/(2*ΔC) * (η(p, ΔC, ρ_max, r) + logterm + t2)
    val += -r^2 / (2*ΔC) * (η(p-2, ΔC, ρ_max, r) + logterm + t3)
    val *= 3*(ΔC-r)*r/(2*ρ_max^3)

    println("\n--- p=$p ---")
    Γ = 3*r*(ΔC-r)/(2*ρ_max^3)
    r_max = r
    ρ2_ΔC2_2ΔC = (ρ_max^2 - ΔC^2) / (2*ΔC)
    r2_ΔC2 = r^2/(2*ΔC)
    @show Γ, r_max, ρ2_ΔC2_2ΔC, r2_ΔC2
    @show logterm + t1
    @show logterm + t2
    @show logterm + t3

    return val
end

function upper_bound_2(r, ρ_max, A, p, ΔC)
    return A / (ΔC - ρ_max - r) * (r/(ΔC - ρ_max))^(p+1)
end


# manual integral
r, ρ_max, A, ΔC = 1.0, 3.0, 1.0, 5.0
ps = collect(0:20)

# integral of abs(Plm(cos)...
#e_manual = manual_sum.(Ref(r), Ref(ρ_max), Ref(A), ps, Ref(ΔC))

# let abs(Plm(cos(theta))) <= 1; numerical integral
e2_manual = manual_sum2.(Ref(r), Ref(ρ_max), Ref(A), ps, Ref(ΔC))

# analytic integral
e3_manual = manual_sum3.(Ref(r), Ref(ρ_max), Ref(A), ps, Ref(ΔC))

# simplification
e4_manual = manual_sum4.(Ref(r), Ref(ρ_max), Ref(A), ps, Ref(ΔC))

# simplification
e5_manual = manual_sum5.(Ref(r), Ref(ρ_max), Ref(A), ps, Ref(ΔC))
e5_manual_rel = manual_sum5_rel.(Ref(r), Ref(ρ_max), Ref(A), ps, Ref(ΔC))

ub2 = upper_bound_2.(Ref(r), Ref(ρ_max), Ref(A), ps, Ref(ΔC))

p = ps[end]

#@btime upper_bound_1($r,$ρ_max,$A,$p,$ΔC)
#@btime upper_bound_2($r,$ρ_max,$A,$p,$ΔC)

#--- check FastMultipole function ---#

r_min = ΔC - r
r_max = r
ρ_min = ΔC - ρ_max
Pmax = 20
ε_rel = 0.0
this_p = FastMultipole.get_P(r_min, r_max, ρ_min, ρ_max, ΔC^2, Pmax, ε_rel, FastMultipole.UniformUnequalSpheres())

using LegendrePolynomials
using FLOWMath

function theta_integral(n)
    θs = range(0, stop=pi, length=100)
    y2 = zeros(length(θs))
    for i in eachindex(θs)
        y2[i] = Plm(cos(θs[i]),n,0) * sin(θs[i])
    end
    return trapz(θs, y2)
end

function manual_integral(ρ_max, n)
    ρ̂s = range(0, stop=ρ_max, length=100)
    y1 = zeros(length(ρ̂s))
    for i in eachindex(y1)
        ρ̂ = ρ̂s[i]
        y1[i] = ρ̂^(n+2)
    end
    return trapz(ρ̂s, y1) * theta_integral(n)
end

function manual_sum(r, ρ_max, A, p)
    println("Manual Sum: p=$p")
    val = 0.0
    ρ_max_3 = ρ_max^3
    rinv = 1 / r
    for n in p+1:100
        val += manual_integral(ρ_max, n) * rinv^n
    end
    return 3*A/(2*ρ_max_3*r) * val
end

function theta_integral2(n)
    return sqrt(2/(2*n+1))
end

function manual_integral2(ρ_max, n)
    ρ̂s = range(0, stop=ρ_max, length=100)
    y1 = zeros(length(ρ̂s))
    for i in eachindex(y1)
        ρ̂ = ρ̂s[i]
        y1[i] = ρ̂^(n+2)
    end
    return trapz(ρ̂s, y1) * theta_integral2(n)
end

function manual_sum2(r, ρ_max, A, p)
    println("Manual Sum: p=$p")
    val = 0.0
    ρ_max_3 = ρ_max^3
    rinv = 1 / r
    for n in p+1:100
        val += manual_integral2(ρ_max, n) * rinv^n
    end
    return 3*A/(2*ρ_max_3*r) * val
end

function theta_integral3(n)
    return 2.0
end

function manual_integral3(ρ_max, n)
    ρ̂s = range(0, stop=ρ_max, length=100)
    y1 = zeros(length(ρ̂s))
    for i in eachindex(y1)
        ρ̂ = ρ̂s[i]
        y1[i] = ρ̂^(n+2)
    end
    return trapz(ρ̂s, y1) * theta_integral3(n)
end

function manual_sum3(r, ρ_max, A, p)
    println("Manual Sum: p=$p")
    val = 0.0
    ρ_max_3 = ρ_max^3
    rinv = 1 / r
    for n in p+1:100
        val += manual_integral3(ρ_max, n) * rinv^n
    end
    return 3*A/(2*ρ_max_3*r) * val
end

function upper_bound(r, ρ_max, A, p)
    t1 = 3*A / (p+4)
    t2 = sqrt(1/(r*p+6))
    t3 = 1 / (r - ρ_max)
    t4 = (ρ_max / r)^(p+1)
    return t1*t2*t3*t4
end

r, ρ_max, A = 4.0, 3.0, 1.0
ps = collect(0:10)

e_manual = manual_sum.(Ref(r), Ref(ρ_max), Ref(A), ps)
e2_manual = manual_sum2.(Ref(r), Ref(ρ_max), Ref(A), ps)
e3_manual = manual_sum3.(Ref(r), Ref(ρ_max), Ref(A), ps)
e_ub = upper_bound.(Ref(r), Ref(ρ_max), Ref(A), ps)


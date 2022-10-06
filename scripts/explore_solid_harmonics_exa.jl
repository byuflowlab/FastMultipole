using LegendrePolynomials
import FLOWFMM as fmm

function evalMultipole(rho, alpha, beta, P)
    Ylm = zeros(Complex{Float64},(P+1)^2)
    YlmTheta = zeros(Complex{Float64},(P+1)^2)
    y,x = sincos(alpha)
    invY = y == 0 ? 0 : 1 / y
    fact = 1
    pl = 1
    rhom = 1
    ei = exp(im * beta)
    eim = 1.0
    for m=0:P
        p = pl
        lpl = m * m + 2 * m + 1
        lml = m * m + 1
        Ylm[lpl] = rhom * p * eim
        Ylm[lml] = conj(Ylm[lpl])
        p1 = p
        p = x * (2 * m + 1) * p1
        YlmTheta[lpl] = rhom * (p - (m + 1) * x * p1) * invY * eim
        rhom *= rho
        rhol = rhom
        for l=m+1:P
            lpm = l * l + l + m + 1
            lmm = l * l + l - m + 1
            rhol /= -(l + m)
            Ylm[lpm] = rhol * p * eim
            Ylm[lmm] = conj(Ylm[lpm])
            p2 = p1
            p1 = p
            p = (x * (2 * l + 1) * p1 - (l + m) * p2) / (l - m + 1)
            YlmTheta[lpm] = rhol * ((l - m + 1) * p - (l + 1) * x * p1) * invY * eim
            rhol *= rho
        end
        rhom /= -(2 * m + 2) * (2 * m + 1)
        pl = -pl * fact * y
        fact += 2
        eim *= ei
    end
    return Ylm
end

function evalMultipole_man(rho, alpha, beta, P)
    harmonics = zeros(Complex{Float64}, (P+1)^2)
    for l in 0:P
        for m in 0:l
            i_p = l^2 + l + m + 1
            i_m = l^2 + l - m + 1
            this_harmonic = (-1)^(l) * rho^l / factorial(l + m) * Plm(cos(alpha), l, m) * exp(im * m * beta)
            harmonics[i_p] = this_harmonic
            harmonics[i_m] = conj(this_harmonic)
        end
    end
    return harmonics
end

dx = [0.8, pi/6, pi/4]
P = 3
regs_exa = evalMultipole(dx..., P)
regs_man = evalMultipole_man(dx..., P)

@show sum(abs.(regs_exa - regs_man))

function evalLocal(rho, theta, phi, P)
    Ylm = zeros(Complex{Float64},(P+1)^2)
    y, x = sincos(theta)
    fact = 1
    pl = 1
    invR = -1.0 / rho
    rhom = -invR
    ei = exp(im * phi)
    eim = 1.0
    for m=0:P
        p = pl
        npl = m * m + 2 * m + 1
        nml = m * m + 1
        Ylm[npl] = rhom * p * eim
        Ylm[nml] = conj(Ylm[npl])
        p1 = p
        p = x * (2 * m + 1) * p1
        rhom *= invR
        rhon = rhom
        for l=m+1:P
            npm = l * l + l + m + 1
            nmm = l * l + l - m + 1
            Ylm[npm] = rhon * p * eim
            Ylm[nmm] = conj(Ylm[npm])
            p2 = p1
            p1 = p
            p = (x * (2 * l + 1) * p1 - (l + m) * p2) / (l - m + 1)
            rhon *= invR * (l - m + 1)
        end
        pl = -pl * fact * y
        fact += 2
        eim *= ei
    end
    return Ylm
end

function evalLocal_man(rho,theta,phi,P)
    harmonics = zeros(Complex{Float64},(P+1)^2)
    for l=0:P
        for m=0:l
            i_p = l^2 + l + m + 1
            i_m = l^2 + l - m + 1
            harmonics[i_p] = (-1)^l * factorial(l-m) / rho^(l+1) * Plm(cos(theta),l,m) * exp(im * m * phi)
            harmonics[i_m] = conj(harmonics[i_p])
        end
    end
    return harmonics
end

irregs_exa = evalLocal(dx..., P)
irregs_man = evalLocal_man(dx..., P)
irregs_fmm = zeros(Complex{Float64},(P+1)^2)
fmm.irregular_harmonic!(irregs_fmm, dx..., P)

@show sum(abs.(irregs_exa - irregs_man))
@show sum(abs.(irregs_man - irregs_fmm))

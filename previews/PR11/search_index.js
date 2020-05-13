var documenterSearchIndex = {"docs":
[{"location":"generated/potential/#Potential-1","page":"4-Potential","title":"4-Potential","text":"","category":"section"},{"location":"generated/potential/#","page":"4-Potential","title":"4-Potential","text":"function EB(r, t, par)\n    x, y, z = r[1], r[2], r[3]\n    ElectricField = E(x, y, z, par)\n    Ex = ElectricField[1]\n    Ey = ElectricField[2]\n\n    B_x = Bx(Ey, par)\n    B_y = By(Ex, par)\n    B_z = Bz(Ex, Ey, x, y, z, par)\n\n    MagneticField = real(Vec3(B_x, B_y, B_z))\n\n    temporal_part = real(g(z, t, par))\n    ElectricField *= temporal_part\n    MagneticField *= temporal_part\n\n    return ElectricField, MagneticField\nend","category":"page"},{"location":"generated/potential/#","page":"4-Potential","title":"4-Potential","text":"","category":"page"},{"location":"generated/potential/#","page":"4-Potential","title":"4-Potential","text":"This page was generated using Literate.jl.","category":"page"},{"location":"generated/electricfield/#Electric-Field-1","page":"Electric Field","title":"Electric Field","text":"","category":"section"},{"location":"generated/electricfield/#","page":"Electric Field","title":"Electric Field","text":"function Ex(laser, x, y, z, r) end\n\nEy(laser, coords...) = laser.ξy / laser.ξx * Ex(laser, coords...)\nEy(laser, Ex) = laser.ξy / laser.ξx * Ex\n\nfunction Ez(laser, x, y, z, r)\n    E_x = Ex(laser, x, y, z, r)\n    Ez(laser, E_x, Ey(laser, E_x), x, y, z, r)\nend\n\nfunction E(x, y, z, laser)\n    r = hypot(x, y)\n\n    E_x = Ex(laser, x, y, z, r)\n    E_y = Ey(laser, E_x)\n    E_z = Ez(laser, E_x, E_y, x, y, z, r)\n\n    real(Vec3(E_x, E_y, E_z))\nend\n\nE(x, y, z, t, laser) = E(x, y, z, laser) * g(z, t, laser)\n\nE(r, t, laser) = E(r[1], r[2], r[3], t, laser)","category":"page"},{"location":"generated/electricfield/#","page":"Electric Field","title":"Electric Field","text":"","category":"page"},{"location":"generated/electricfield/#","page":"Electric Field","title":"Electric Field","text":"This page was generated using Literate.jl.","category":"page"},{"location":"temporal-profiles/#Temporal-Profiles-1","page":"Temporal Profiles","title":"Temporal Profiles","text":"","category":"section"},{"location":"temporal-profiles/#","page":"Temporal Profiles","title":"Temporal Profiles","text":"The complete expression of the electromaginetic fields describing the laser pulse is given by the combination of the spatial part together with the time dependence.","category":"page"},{"location":"temporal-profiles/#","page":"Temporal Profiles","title":"Temporal Profiles","text":"Modules = [LaserTypes]\nPages = [\"src/envelopes.jl\"]","category":"page"},{"location":"temporal-profiles/#LaserTypes.ConstantProfile","page":"Temporal Profiles","title":"LaserTypes.ConstantProfile","text":"struct ConstantProfile\n\nThis is the trivial profile\n\nenvelope(z t) = 1\n\nwhich gives an infinite duration for the electromagnetic field configuration (we cannot call this a laser pusle since a pulse implicitely has a finite duration).\n\n\n\n\n\n","category":"type"},{"location":"temporal-profiles/#LaserTypes.GaussianProfile","page":"Temporal Profiles","title":"LaserTypes.GaussianProfile","text":"GaussianProfile{V,T,L}\n\nThis envelope provides a finite duration for the laser pulse and thus can provide a more realistic description of an actual laser pulse.\n\nenvelope(z t) = expleft-left(fract-fracz-z_Fctau_0right)^2right\n\nwhere\n\nc is the speed of light\nτ is the duration of the pulse (FWHM) and has the default value 18.02fs\nt₀ is the origin of the time axis and it is 0 by default\nz₀ is the initial position of the intensity peak and has the default value -4*τ*c\n\n\n\n\n\n","category":"type"},{"location":"temporal-profiles/#LaserTypes.QuasiRectangularProfile","page":"Temporal Profiles","title":"LaserTypes.QuasiRectangularProfile","text":"QuasiRectangularProfile{V,T,L}\n\nThis envelope provides a finite duration for the laser pulse and thus can provide a more realistic description of an actual laser pulse.\n\nenvelope(z t) =\n\nwhere\n\nc is the speed of light\nτ is the duration of the pulse (FWHM) and has the default value 18.02fs\nt₀ is the origin of the time axis and it is 0 by default\nz₀ is the initial position of the intensity peak and has the default value -4*τ*c\nΔz is the width of the flat part of the profile and the default value 10*τ*c\n\n\n\n\n\n","category":"type"},{"location":"temporal-profiles/#LaserTypes.g-Tuple{Any,Any,Any}","page":"Temporal Profiles","title":"LaserTypes.g","text":"g(z, t, par)\n\nThe time dependence of the fields is given by\n\ng(z t) = cos(omega t) envelope(z t)\n\nwhere\n\nomega is the angular frequency of the laser pulse\nenvelope(z t) is a function that can be used to control the duration of the pulse\n\n\n\n\n\n","category":"method"},{"location":"generated/envelopes/#Envelopes-1","page":"Envelopes","title":"Envelopes","text":"","category":"section"},{"location":"generated/envelopes/#","page":"Envelopes","title":"Envelopes","text":"\"\"\"\n    struct ConstantProfile\n\nThis is the trivial profile","category":"page"},{"location":"generated/envelopes/#","page":"Envelopes","title":"Envelopes","text":"math envelope(z, t) = 1","category":"page"},{"location":"generated/envelopes/#","page":"Envelopes","title":"Envelopes","text":"which gives an infinite duration for the electromagnetic field configuration (we cannot\ncall this a laser pusle since a pulse implicitely has a finite duration).\n\"\"\"\nstruct ConstantProfile end\n\n@doc \"\"\"\n    GaussianProfile{V,T,L}\n\nThis envelope provides a finite duration for the laser pulse and thus can provide a more\nrealistic description of an actual laser pulse.","category":"page"},{"location":"generated/envelopes/#","page":"Envelopes","title":"Envelopes","text":"math envelope(z, t) = \\exp\\left[-\\left(\\frac{t-\\frac{z-zF}{c}}{\\tau0}\\right)^2\\right],","category":"page"},{"location":"generated/envelopes/#","page":"Envelopes","title":"Envelopes","text":"where\n- `c` is the speed of light\n- `τ` is the duration of the pulse (FWHM) and has the default value 18.02fs\n- `t₀` is the origin of the time axis and it is 0 by default\n- `z₀` is the initial position of the intensity peak and has the default value `-4*τ*c`\n\"\"\"\nGaussianProfile\n\n@with_kw struct GaussianProfile{V,T,L}\n    c::V = c_0\n    τ::T = 18.02u\"fs\"\n    t₀::T = zero(τ)\n    z₀::L = -4*τ*c\nend\n\n@doc \"\"\"\n    QuasiRectangularProfile{V,T,L}\n\nThis envelope provides a finite duration for the laser pulse and thus can provide a more\nrealistic description of an actual laser pulse.","category":"page"},{"location":"generated/envelopes/#","page":"Envelopes","title":"Envelopes","text":"math envelope(z, t) =","category":"page"},{"location":"generated/envelopes/#","page":"Envelopes","title":"Envelopes","text":"where\n- `c` is the speed of light\n- `τ` is the duration of the pulse (FWHM) and has the default value 18.02fs\n- `t₀` is the origin of the time axis and it is 0 by default\n- `z₀` is the initial position of the intensity peak and has the default value `-4*τ*c`\n- `Δz` is the width of the flat part of the profile and the default value `10*τ*c`\n\"\"\"\nQuasiRectangularProfile\n\n@with_kw struct QuasiRectangularProfile{V,T,L}\n    c::V = c_0\n    τ::T = 18.02u\"fs\"\n    t₀::T = zero(τ)\n    z₀::L = -4*τ*c\n    Δz::L = 10*τ*c\nend\n\n\"\"\"\n    g(z, t, par)\n\nThe time dependence of the fields is given by","category":"page"},{"location":"generated/envelopes/#","page":"Envelopes","title":"Envelopes","text":"math g(z, t) = \\cos(\\omega t) envelope(z, t),","category":"page"},{"location":"generated/envelopes/#","page":"Envelopes","title":"Envelopes","text":"where\n- ``\\\\omega`` is the angular frequency of the laser pulse\n- ``envelope(z, t)`` is a function that can be used to control the duration of the pulse\n\n\"\"\"\nfunction g(z, t, par)\n    @unpack profile, ω = par\n\n    cos(ω*t) * envelope(profile, z, t)\nend\n\nenvelope(::ConstantProfile, z, t) = 1\n\nfunction envelope(profile::GaussianProfile, z, t)\n    @unpack c, τ, t₀, z₀ = profile\n    φ = (t - t₀) - (z - z₀) / c\n\n    cosh(uconvert(NoUnits, φ / τ)^2)\nend\n\nfunction envelope(profile::QuasiRectangularProfile, z, t)\n    @unpack c, τ, t₀, z₀, Δz = profile\n    φ = (t - t₀) - (z - z₀) / c\n\n    if φ < zero(φ)\n        cosh(uconvert(NoUnits, φ / τ)^2)\n    elseif φ < Δz\n        one(φ)\n    else\n        cosh(uconvert(NoUnits, (φ - Δz/c) / τ)^2)\n    end\nend","category":"page"},{"location":"generated/envelopes/#","page":"Envelopes","title":"Envelopes","text":"","category":"page"},{"location":"generated/envelopes/#","page":"Envelopes","title":"Envelopes","text":"This page was generated using Literate.jl.","category":"page"},{"location":"generated/magneticfield/#Magnetic-Field-1","page":"Magnetic Field","title":"Magnetic Field","text":"","category":"section"},{"location":"generated/magneticfield/#","page":"Magnetic Field","title":"Magnetic Field","text":"Bx(laser, c1, c2, rest_of_coords...) = -1/laser.c * Ey(laser, c1, c2, rest_of_coords...)\nBx(laser, Ey) = -1/laser.c * Ey\n\nBy(laser, c1, c2, rest_of_coords...) = 1/laser.c * Ex(laser, c1, c2, rest_of_coords...)\nBy(laser, Ex) = 1/laser.c * Ex\n\nfunction Bz(laser, x, y, z, r)\n    E_x = Ex(laser, x, y, z, r)\n    Bz(laser, E_x, Ey(laser, E_x), x, y, z, r)\nend\n\nfunction B(x, y, z, laser)\n    r = hypot(x, y)\n\n    E_x = Ex(laser, x, y, z, r)\n    E_y = Ey(laser, E_x)\n\n    B_x = Bx(laser, E_y)\n    B_y = By(laser, E_x)\n    B_z = Bz(laser, E_x, E_y, x, y, z, r)\n\n    real(Vec3(B_x, B_y, B_z))\nend\n\nB(x, y, z, t, laser) = B(x, y, z, laser) * g(z, t, laser)\n\nB(r, t, laser) = B(r[1], r[2], r[3], t, laser)","category":"page"},{"location":"generated/magneticfield/#","page":"Magnetic Field","title":"Magnetic Field","text":"","category":"page"},{"location":"generated/magneticfield/#","page":"Magnetic Field","title":"Magnetic Field","text":"This page was generated using Literate.jl.","category":"page"},{"location":"generated/laguerre-gauss/#Laguerre-Gauss-profile-1","page":"Laguerre-Gauss profile","title":"Laguerre-Gauss profile","text":"","category":"section"},{"location":"generated/laguerre-gauss/#","page":"Laguerre-Gauss profile","title":"Laguerre-Gauss profile","text":"@doc \"\"\"\n    struct LaguerreGaussLaser{V,Q,M,L,F,C,T,P,I,W,K,E,R}\n\nThe `LaguerreGaussLaser` is defined by the following independent parameters\n- `c` is the speed of light in vaccum, with the default value being in SI (`c_0` from the CODATA2018 in the [PhysicalConstants](https://github.com/JuliaPhysics/PhysicalConstants.jl) package)\n- `q` is the electric charge, with the default value being the one for the electron in SI (`-e` from the CODATA2018 in the [PhysicalConstants](https://github.com/JuliaPhysics/PhysicalConstants.jl) package)\n- `m_q` is the mass of the charge, with the default value being the one for the electron in SI (`m_e` from the CODATA2018 in the [PhysicalConstants](https://github.com/JuliaPhysics/PhysicalConstants.jl) package)\n- `λ` is the laser wavelangth with the default value 0.8μm\n- `a₀` is the normalized vector potential (defined as ``a_0=\\\\frac{eA}{m_e c^2}``)\n- `φ₀` is the initial phase with the default value 0.0\n- `w₀` is the beam radius at the Rayleigh range or [beam waist](https://en.wikipedia.org/wiki/Gaussian_beam#Beam_waist) with the default value 58.0μm\n- `ξx` and `ξy` give the polarization and have the default value `1.0 + 0im` and `0.0 + 0im`\n- `profile` is the temporal profile of the pulse and the default one is a Gaussian one\n\nDuring the initialization of the `LaguerreGaussLaser` type, some useful derived values are\nalso computed\n- `ω` is the angular frequency and is given by ``2\\\\pi \\\\frac{c}{\\\\lambda}``\n- `k` is the wavenumber and is given by ``\\\\frac{2\\\\pi}{\\\\lambda}``\n- `z_R` is the [Rayleigh range](https://en.wikipedia.org/wiki/Rayleigh_length) and is given by ``\\\\frac{k w_0^2}{2}``\n- `T₀` is the laser period and is given by ``\\\\frac{2\\\\pi}{\\\\omega}``\n- `E₀` is the amplitude of the electric field and is given by ``a_0\\\\frac{m_q c \\\\omega}{q}``\n\"\"\"\nLaguerreGaussLaser\n\n@with_kw struct LaguerreGaussLaser{V,Q,M,L,F,C,T,P,I,W,K,E,R}","category":"page"},{"location":"generated/laguerre-gauss/#","page":"Laguerre-Gauss profile","title":"Laguerre-Gauss profile","text":"independent values","category":"page"},{"location":"generated/laguerre-gauss/#","page":"Laguerre-Gauss profile","title":"Laguerre-Gauss profile","text":"    c::V = c_0\n    q::Q = -e\n    m_q::M = m_e\n    λ::L = 0.8u\"μm\"\n    a₀::F = 1.0\n    φ₀::F = 0.0\n    w₀::L = 58.0u\"μm\"\n    ξx::C = 1.0 + 0im\n    ξy::C = 0.0 + 0im\n    @assert hypot(ξx, ξy) ≈ 1\n    profile::P = GaussianProfile(c=c)\n    p::I = 1\n    m::I = 0","category":"page"},{"location":"generated/laguerre-gauss/#","page":"Laguerre-Gauss profile","title":"Laguerre-Gauss profile","text":"derived values","category":"page"},{"location":"generated/laguerre-gauss/#","page":"Laguerre-Gauss profile","title":"Laguerre-Gauss profile","text":"    ω::W = 2π * c / λ; @assert ω ≈ 2π * c / λ\n    k::K = 2π / λ; @assert k ≈ 2π / λ\n    z_R::L = w₀^2 * k / 2; @assert z_R ≈ w₀^2 * k / 2\n    T₀::T = 2π / ω; @assert T₀ ≈ 2π / ω\n    E₀::E = a₀ * m_q * c * ω / abs(q); @assert E₀ ≈ a₀ * m_q * c * ω / abs(q)\n    Nₚₘ::R = √(pochhammer(p+1, abs(m))); @assert Nₚₘ ≈ √(pochhammer(p+1, abs(m)))\nend\n\n\"\"\"\n    GaussLaser(laser::LaguerreGaussLaser)\n\nConvert a `LaguerreGaussLaser` to a `GaussLaser` with the same parameters.\n\"\"\"\nfunction GaussLaser(laser::LaguerreGaussLaser)\n    @unpack c, q, m_q, λ, a₀, φ₀, w₀, ξx, ξy, profile, ω, k, z_R, T₀, E₀ = laser\n    GaussLaser(c, q, m_q, λ, a₀, φ₀, w₀, ξx, ξy, profile, ω, k, z_R, T₀, E₀)\nend\n\nfunction Ex(laser::LaguerreGaussLaser, x, y, z, r)\n    @unpack Nₚₘ, w₀, φ₀, z_R, ξx, p, m = laser\n    wz = w(z, laser)\n    Rz = R(z, z_R)\n    gauss_laser = GaussLaser(laser)\n    Eg = Ex(gauss_laser, z, r)\n    σ = (r/wz)^2\n    mₐ = abs(m)\n    φ = atan(x, y)\n\n    ξx*Eg*Nₚₘ*(r*√2/wz)^mₐ*_₁F₁(-p, mₐ+1, 2σ)*exp(im*((2p+mₐ)*atan(z, z_R)-m*φ-φ₀))\nend\n\nfunction Ez(laser::LaguerreGaussLaser, Ex, Ey, x, y, z, r)\n    @unpack k, z_R, p, m = laser\n    wz = w(z, laser)\n    mₐ = abs(m)\n    ∓ = m > 0 ? (-) : +\n\n    -im / k * (\n        (-2*(1+im*(z/z_R))/wz^2\n        + 4p/(((mₐ+1)*wz^2) * _₁F₁(-p+1, mₐ+2, 2r^2/wz^2))) * (x*Ex + y*Ey)\n        - (!iszero(m) ? mₐ/(x+im*y) * (Ex ∓ im*Ey) : zero(typeof(Ex))/oneunit(typeof(x)))\n        )\nend\n\nfunction Bz(laser::LaguerreGaussLaser, Ex, Ey, x, y, z, r)\n    @unpack k, z_R, p, m, c = laser\n    wz = w(z, laser)\n    mₐ = abs(m)\n    ∓ = m > 0 ? (-) : +\n\n    -im / (c*k) * (\n        (-2*(1+im*(z/z_R))/wz^2\n        + 4p/(((mₐ+1)*wz^2) * _₁F₁(-p+1, mₐ+2, 2r^2/wz^2))) * (y*Ex + x*Ey)\n        - (!iszero(m) ? mₐ/(x+im*y) * (Ey ∓ im*Ex) : zero(typeof(Ex))/oneunit(typeof(x)))\n        )\nend","category":"page"},{"location":"generated/laguerre-gauss/#","page":"Laguerre-Gauss profile","title":"Laguerre-Gauss profile","text":"","category":"page"},{"location":"generated/laguerre-gauss/#","page":"Laguerre-Gauss profile","title":"Laguerre-Gauss profile","text":"This page was generated using Literate.jl.","category":"page"},{"location":"generated/gauss/#Gaussian-profile-1","page":"Gaussian profile","title":"Gaussian profile","text":"","category":"section"},{"location":"generated/gauss/#","page":"Gaussian profile","title":"Gaussian profile","text":"@doc \"\"\"\n    struct GaussLaser{V,Q,M,L,F,C,T,P,W,K,E}\n\nThe `GaussLaser` is defined by the following independent parameters\n- `c` is the speed of light in vaccum, with the default value being in SI (`c_0` from the CODATA2018 in the [PhysicalConstants](https://github.com/JuliaPhysics/PhysicalConstants.jl) package)\n- `q` is the electric charge, with the default value being the one for the electron in SI (`-e` from the CODATA2018 in the [PhysicalConstants](https://github.com/JuliaPhysics/PhysicalConstants.jl) package)\n- `m_q` is the mass of the charge, with the default value being the one for the electron in SI (`m_e` from the CODATA2018 in the [PhysicalConstants](https://github.com/JuliaPhysics/PhysicalConstants.jl) package)\n- `λ` is the laser wavelangth with the default value 0.8μm\n- `a₀` is the normalized vector potential (defined as ``a_0=\\\\frac{eA}{m_e c^2}``)\n- `φ₀` is the initial phase with the default value 0.0\n- `w₀` is the beam radius at the Rayleigh range or [beam waist](https://en.wikipedia.org/wiki/Gaussian_beam#Beam_waist) with the default value 58.0μm\n- `ξx` and `ξy` give the polarization and have the default value `1.0 + 0im` and `0.0 + 0im`\n- `profile` is the temporal profile of the pulse and the default one is a Gaussian one\n\nDuring the initialization of the `GaussLaser` type, some useful derived values are\nalso computed\n- `ω` is the angular frequency and is given by ``2\\\\pi \\\\frac{c}{\\\\lambda}``\n- `k` is the wavenumber and is given by ``\\\\frac{2\\\\pi}{\\\\lambda}``\n- `z_R` is the [Rayleigh range](https://en.wikipedia.org/wiki/Rayleigh_length) and is given by ``\\\\frac{k w_0^2}{2}``\n- `T₀` is the laser period and is given by ``\\\\frac{2\\\\pi}{\\\\omega}``\n- `E₀` is the amplitude of the electric field and is given by ``a_0\\\\frac{m_q c \\\\omega}{q}``\n\"\"\"\nGaussLaser\n\n@with_kw struct GaussLaser{V,Q,M,L,F,C,T,P,W,K,E}","category":"page"},{"location":"generated/gauss/#","page":"Gaussian profile","title":"Gaussian profile","text":"independent values","category":"page"},{"location":"generated/gauss/#","page":"Gaussian profile","title":"Gaussian profile","text":"    c::V = c_0\n    q::Q = -e\n    m_q::M = m_e\n    λ::L = 0.8u\"μm\"\n    a₀::F = 1.0\n    φ₀::F = 0.0\n    w₀::L = 58.0u\"μm\"\n    ξx::C = 1.0 + 0im\n    ξy::C = 0.0 + 0im\n    @assert hypot(ξx, ξy) ≈ 1\n    profile::P = GaussianProfile(c=c)","category":"page"},{"location":"generated/gauss/#","page":"Gaussian profile","title":"Gaussian profile","text":"derived values","category":"page"},{"location":"generated/gauss/#","page":"Gaussian profile","title":"Gaussian profile","text":"    ω::W = 2π * c / λ; @assert ω ≈ 2π * c / λ\n    k::K = 2π / λ; @assert k ≈ 2π / λ\n    z_R::L = w₀^2 * k / 2; @assert z_R ≈ w₀^2 * k / 2\n    T₀::T = 2π / ω; @assert T₀ ≈ 2π / ω\n    E₀::E = a₀ * m_q * c * ω / abs(q); @assert E₀ ≈ a₀ * m_q * c * ω / abs(q)\nend\n\n@inline Ex(laser::GaussLaser, x, y, z, r) = Ex(laser, z, r)\n\nfunction Ex(laser::GaussLaser, z, r)\n    @unpack E₀, w₀, k, z_R, φ₀, ξx = laser\n    wz = w(z, laser)\n    Rz = R(z, z_R)\n\n    ξx * E₀ * w₀/wz * exp(-im*k*z - (r/wz)^2 - im*((k*r^2)/(2Rz) - atan(z, z_R) - φ₀))\nend\n\n@inline Ez(laser::GaussLaser, Ex, Ey, x, y, z, r) = Ez(laser, Ex, Ey, x, y, z)\n\nfunction Ez(laser::GaussLaser, Ex, Ey, x, y, z)\n    @unpack k, z_R = laser\n    wz = w(z, laser)\n\n    2(im - z/z_R) / (k*wz^2) * (x*Ex + y*Ey)\nend\n\n@inline Bz(laser::GaussLaser, Ex, Ey, x, y, z, r) = Bz(laser, Ex, Ey, x, y, z)\n\nfunction Bz(laser::GaussLaser, Ex, Ey, x, y, z)\n    @unpack k, z_R, c = laser\n    wz = w(z, laser)\n\n    2(im - z/z_R) / (c*k*wz^2) * (y*Ex - x*Ey)\nend","category":"page"},{"location":"generated/gauss/#","page":"Gaussian profile","title":"Gaussian profile","text":"","category":"page"},{"location":"generated/gauss/#","page":"Gaussian profile","title":"Gaussian profile","text":"This page was generated using Literate.jl.","category":"page"},{"location":"laguerre-gauss/#Laguerre-Gauss-profile-1","page":"Laguerre-Gauss profile","title":"Laguerre-Gauss profile","text":"","category":"section"},{"location":"laguerre-gauss/#","page":"Laguerre-Gauss profile","title":"Laguerre-Gauss profile","text":"The LaguerreGaussLaser type contains the parameters required for describing a Laguerre-Gauss laser pulse. The spatial profile of the electric field is given by","category":"page"},{"location":"laguerre-gauss/#","page":"Laguerre-Gauss profile","title":"Laguerre-Gauss profile","text":"newcommandiimathrmi                 imaginary unit\nbeginaligned\n    E_x(rzvarphi) = xi_x E_g N_pm left(fracsqrt2rw(z)right)^m ₁F₁left(-p m+1 2left(fracrw(z)right)^2right) expleftiileft((2p+m)arctgfraczz_R-mvarphi - phi_0right)right \n    E_y(rzvarphi) = fracxi_yxi_x E_x(rzvarphi) \n    E_z(rzvarphi) = -fraciik leftleft-2frac1+iifraczz_0w^2(z) + frac4p(m+1)w(z)^2 ₁F₁^-1left(-p m+1 2left(fracrw(z)right)^2right)right(x E_x + y E_y) - fracmx+ii y (E_x mp ii E_y)right\nendaligned","category":"page"},{"location":"laguerre-gauss/#","page":"Laguerre-Gauss profile","title":"Laguerre-Gauss profile","text":"where","category":"page"},{"location":"laguerre-gauss/#","page":"Laguerre-Gauss profile","title":"Laguerre-Gauss profile","text":"E_g is the amplitude of the electric field form the Gaussian profile, E_x(rz)\n₁F₁(-pm+1z) is the confluent hypergeometric function\nxi_x and xi_y give the polarization (choosing xi_x=1 and xi_y=0 for example, would give a linearly polarized field along the x axis while taking both 1sqrt2 would give circular polarization)\nw = w(z) is the beam radius at z\nw_0 is the beam waist\nk is the wavenumber\nR is the radius of curvature\narctgfraczz_R is the Gouy phase\nphi_0 is the initial phase","category":"page"},{"location":"laguerre-gauss/#","page":"Laguerre-Gauss profile","title":"Laguerre-Gauss profile","text":"The magnetic field field is given by","category":"page"},{"location":"laguerre-gauss/#","page":"Laguerre-Gauss profile","title":"Laguerre-Gauss profile","text":"newcommandiimathrmi                 imaginary unit\nbeginaligned\n    B_x(rz) = -frac1cE_y(rz)\n    B_y(rz) = frac1cE_x(rz)\n    B_z (rz) = -fraciick leftleft-2frac1+iifraczz_0w^2(z) + frac4p(m+1)w(z)^2 ₁F₁^-1left(-p m+1 2left(fracrw(z)right)^2right)right(y E_x + x E_y) - fracmx+ii y (E_y mp ii E_x)right\nendaligned","category":"page"},{"location":"laguerre-gauss/#","page":"Laguerre-Gauss profile","title":"Laguerre-Gauss profile","text":"Modules = [LaserTypes]\nPages = [\"src/laguerre-gauss.jl\"]","category":"page"},{"location":"laguerre-gauss/#LaserTypes.GaussLaser-Tuple{LaguerreGaussLaser}","page":"Laguerre-Gauss profile","title":"LaserTypes.GaussLaser","text":"GaussLaser(laser::LaguerreGaussLaser)\n\nConvert a LaguerreGaussLaser to a GaussLaser with the same parameters.\n\n\n\n\n\n","category":"method"},{"location":"laguerre-gauss/#LaserTypes.LaguerreGaussLaser","page":"Laguerre-Gauss profile","title":"LaserTypes.LaguerreGaussLaser","text":"struct LaguerreGaussLaser{V,Q,M,L,F,C,T,P,I,W,K,E,R}\n\nThe LaguerreGaussLaser is defined by the following independent parameters\n\nc is the speed of light in vaccum, with the default value being in SI (c_0 from the CODATA2018 in the PhysicalConstants package)\nq is the electric charge, with the default value being the one for the electron in SI (-e from the CODATA2018 in the PhysicalConstants package)\nm_q is the mass of the charge, with the default value being the one for the electron in SI (m_e from the CODATA2018 in the PhysicalConstants package)\nλ is the laser wavelangth with the default value 0.8μm\na₀ is the normalized vector potential (defined as a_0=fraceAm_e c^2)\nφ₀ is the initial phase with the default value 0.0\nw₀ is the beam radius at the Rayleigh range or beam waist with the default value 58.0μm\nξx and ξy give the polarization and have the default value 1.0 + 0im and 0.0 + 0im\nprofile is the temporal profile of the pulse and the default one is a Gaussian one\n\nDuring the initialization of the LaguerreGaussLaser type, some useful derived values are also computed\n\nω is the angular frequency and is given by 2pi fracclambda\nk is the wavenumber and is given by frac2pilambda\nz_R is the Rayleigh range and is given by frack w_0^22\nT₀ is the laser period and is given by frac2piomega\nE₀ is the amplitude of the electric field and is given by a_0fracm_q c omegaq\n\n\n\n\n\n","category":"type"},{"location":"generated/LaserTypes/#LaserTypes-1","page":"LaserTypes","title":"LaserTypes","text":"","category":"section"},{"location":"generated/LaserTypes/#","page":"LaserTypes","title":"LaserTypes","text":"module LaserTypes\n\nexport E, B, GaussLaser, LaguerreGaussLaser, ConstantProfile, GaussianProfile, QuasiRectangularProfile\n\nusing Parameters\n\nusing Unitful\nusing Parameters\nusing GeometryTypes: Vec3\nimport PhysicalConstants.CODATA2018: c_0, m_e, e\nusing HypergeometricFunctions\n\nconst _₁F₁ = HypergeometricFunctions.drummond1F1\nconst pochhammer = HypergeometricFunctions.pochhammer\n\nfunction w(z, par)\n    @unpack w₀, z_R = par\n\n    w₀ * √(1 + (z/z_R)^2)\nend\n\nR(z, z_R) = z + z_R^2 / z\n\ninclude(\"envelopes.jl\")\ninclude(\"electricfield.jl\")\ninclude(\"magneticfield.jl\")\ninclude(\"potential.jl\")\ninclude(\"gauss.jl\")\ninclude(\"laguerre-gauss.jl\")\n\nend # module","category":"page"},{"location":"generated/LaserTypes/#","page":"LaserTypes","title":"LaserTypes","text":"","category":"page"},{"location":"generated/LaserTypes/#","page":"LaserTypes","title":"LaserTypes","text":"This page was generated using Literate.jl.","category":"page"},{"location":"gauss/#gaussian-1","page":"Gaussian profile","title":"Gaussian profile","text":"","category":"section"},{"location":"gauss/#","page":"Gaussian profile","title":"Gaussian profile","text":"The GaussLaser type contains the parameters required for describing a Gaussian laser pulse. The spatial profile of the electric field is given by","category":"page"},{"location":"gauss/#","page":"Gaussian profile","title":"Gaussian profile","text":"newcommandiimathrmi                 imaginary unit\nbeginaligned\n    E_x (rz) = xi_x E_0 fracw_0w expleft-ii k z -fracr^2w^2 -ii left(frack r^22 R - arctgfraczz_R - phi_0right)right\n    E_y (rz) = fracxi_yxi_x E_x(rz) \n    E_z (rz) = frac2 left(ii -fraczz_Rright)kw^2 (z) leftxE_x (rz) +yE_y(rz)right\nendaligned","category":"page"},{"location":"gauss/#","page":"Gaussian profile","title":"Gaussian profile","text":"where","category":"page"},{"location":"gauss/#","page":"Gaussian profile","title":"Gaussian profile","text":"E_0 is the amplitude of the electric field\nxi_x and xi_y give the polarization (choosing xi_x=1 and xi_y=0 for example, would give a linearly polarized field along the x axis while taking both 1sqrt2 would give circular polarization)\nw = w(z) is the beam radius at z\nw_0 is the beam waist\nk is the wavenumber\nR is the radius of curvature\narctgfraczz_R is the Gouy phase\nphi_0 is the initial phase","category":"page"},{"location":"gauss/#","page":"Gaussian profile","title":"Gaussian profile","text":"The magnetic field field is given by","category":"page"},{"location":"gauss/#","page":"Gaussian profile","title":"Gaussian profile","text":"newcommandiimathrmi                 imaginary unit\nbeginaligned\n    B_x(rz) = -frac1cE_y(rz)\n    B_y(rz) = frac1cE_x(rz)\n    B_z (rz) = frac2 left(ii -fraczz_Rright)ckw^2 (z) yE_x (rz) - xE_y(rz)\nendaligned","category":"page"},{"location":"gauss/#","page":"Gaussian profile","title":"Gaussian profile","text":"Modules = [LaserTypes]\nPages = [\"src/gauss.jl\"]\nOrder = [:type]","category":"page"},{"location":"gauss/#LaserTypes.GaussLaser","page":"Gaussian profile","title":"LaserTypes.GaussLaser","text":"struct GaussLaser{V,Q,M,L,F,C,T,P,W,K,E}\n\nThe GaussLaser is defined by the following independent parameters\n\nc is the speed of light in vaccum, with the default value being in SI (c_0 from the CODATA2018 in the PhysicalConstants package)\nq is the electric charge, with the default value being the one for the electron in SI (-e from the CODATA2018 in the PhysicalConstants package)\nm_q is the mass of the charge, with the default value being the one for the electron in SI (m_e from the CODATA2018 in the PhysicalConstants package)\nλ is the laser wavelangth with the default value 0.8μm\na₀ is the normalized vector potential (defined as a_0=fraceAm_e c^2)\nφ₀ is the initial phase with the default value 0.0\nw₀ is the beam radius at the Rayleigh range or beam waist with the default value 58.0μm\nξx and ξy give the polarization and have the default value 1.0 + 0im and 0.0 + 0im\nprofile is the temporal profile of the pulse and the default one is a Gaussian one\n\nDuring the initialization of the GaussLaser type, some useful derived values are also computed\n\nω is the angular frequency and is given by 2pi fracclambda\nk is the wavenumber and is given by frac2pilambda\nz_R is the Rayleigh range and is given by frack w_0^22\nT₀ is the laser period and is given by frac2piomega\nE₀ is the amplitude of the electric field and is given by a_0fracm_q c omegaq\n\n\n\n\n\n","category":"type"},{"location":"gauss/#","page":"Gaussian profile","title":"Gaussian profile","text":"The electric field is normalized such that in origin (with the default linear polarization) we obtain","category":"page"},{"location":"gauss/#","page":"Gaussian profile","title":"Gaussian profile","text":"vecE(000) = (E_000)","category":"page"},{"location":"gauss/#","page":"Gaussian profile","title":"Gaussian profile","text":"when the temporal profile is ignored (or constant and equal to unity).","category":"page"},{"location":"#LaserTypes.jl-1","page":"LaserTypes.jl","title":"LaserTypes.jl","text":"","category":"section"},{"location":"#","page":"LaserTypes.jl","title":"LaserTypes.jl","text":"This package aims to provide a common interface for different laser types.","category":"page"}]
}

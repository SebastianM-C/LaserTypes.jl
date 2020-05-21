var documenterSearchIndex = {"docs":
[{"location":"covariant/#Covariant-formulation-1","page":"Covariant formulation","title":"Covariant formulation","text":"","category":"section"},{"location":"covariant/#","page":"Covariant formulation","title":"Covariant formulation","text":"Electrodynamics can also be expressed in a covarinat formulaltion. In this case the objects of interest would be 4-vectors and tensors. This package provides a convenience function for computing the electromagnetic tensor F from the corresponding electric and magnetic fields for a given laser type.","category":"page"},{"location":"covariant/#","page":"Covariant formulation","title":"Covariant formulation","text":"Modules = [LaserTypes]\nPages = [\"src/faraday.jl\"]","category":"page"},{"location":"covariant/#LaserTypes.Fμν-Tuple{Any,Any}","page":"Covariant formulation","title":"LaserTypes.Fμν","text":"Fμν(x, laser)\n\nCompute the electromagnetic tensor F in covarinat matrix form at the spacetime 4-vector specified by x. The x = (x^0  ct x^1 x^2 x^3) convention is used.\n\nF^μν =\nleftbeginarraycccc\n    0     -E_x  c   -E_y  c  -E_z  c \nE_x  c        0       -B_z       B_y   \nE_y  c       B_z        0       -B_x   \nE_z  c      -B_y       B_z        0\nendarrayright\n\n\n\n\n\n","category":"method"},{"location":"laguerre-gauss/#Laguerre-Gauss-profile-1","page":"Laguerre-Gauss profile","title":"Laguerre-Gauss profile","text":"","category":"section"},{"location":"laguerre-gauss/#","page":"Laguerre-Gauss profile","title":"Laguerre-Gauss profile","text":"The LaguerreGaussLaser type contains the parameters required for describing a Laguerre-Gauss laser pulse. The spatial profile of the electric field is given by","category":"page"},{"location":"laguerre-gauss/#","page":"Laguerre-Gauss profile","title":"Laguerre-Gauss profile","text":"newcommandiimathrmi                 imaginary unit\nbeginaligned\n    E_x(rzφ) = ξ_x E_g N_pm left(fracsqrt2rw(z)right)^m ₁F₁left(-p m+1 2left(fracrw(z)right)^2right) expleftiileft((2p+m)arctgfraczz_R-mφ - phi_0right)right \n    E_y(rzφ) = fracξ_yξ_x E_x(rzφ) \n    E_z(rzφ) = -fraciik leftleft-2frac1+iifraczz_0w(z)^2 + frac4p(m+1)w(z)^2 ₁F₁^-1left(-p m+1 2left(fracrw(z)right)^2right)right(x E_x + y E_y) - fracmx+ii y (E_x mp ii E_y)right\nendaligned","category":"page"},{"location":"laguerre-gauss/#","page":"Laguerre-Gauss profile","title":"Laguerre-Gauss profile","text":"where","category":"page"},{"location":"laguerre-gauss/#","page":"Laguerre-Gauss profile","title":"Laguerre-Gauss profile","text":"E_g is the amplitude of the electric field form the Gaussian profile, E_x(rz)\nN_pm is a normalization factor given by sqrt(p+1)_m, with (x)_n the Pochhammer symbol\n₁F₁(αβz) is the confluent hypergeometric function of the first kind\np is the radial index of the mode, p  ℤ p  0\nm is the azimuthal index of the mode, m  ℤ\nξ_x and ξ_y give the polarization (choosing ξ_x=1 and ξ_y=0 for example, would give a linearly polarized field along the x axis while taking them 1sqrt2 and mathrmisqrt2 would give right and left-handed circular polarization)\nw(z) is the beam radius at z\nw_0 is the beam waist\nk is the wavenumber\nR is the radius of curvature\narctgfraczz_R is the Gouy phase\nϕ_0 is the initial phase","category":"page"},{"location":"laguerre-gauss/#","page":"Laguerre-Gauss profile","title":"Laguerre-Gauss profile","text":"The magnetic field field is given by","category":"page"},{"location":"laguerre-gauss/#","page":"Laguerre-Gauss profile","title":"Laguerre-Gauss profile","text":"newcommandiimathrmi                 imaginary unit\nbeginaligned\n    B_x(rz) = -frac1cE_y(rz)\n    B_y(rz) = frac1cE_x(rz)\n    B_z (rz) = -fraciick leftleft-2frac1+iifraczz_0w^2(z) + frac4p(m+1)w(z)^2 ₁F₁^-1left(-p m+1 2left(fracrw(z)right)^2right)right(y E_x + x E_y) - fracmx+ii y (E_y mp ii E_x)right\nendaligned","category":"page"},{"location":"laguerre-gauss/#","page":"Laguerre-Gauss profile","title":"Laguerre-Gauss profile","text":"note: Choice of normalization\nThe normalization was chosen such that the normalization factor N_pm is adimensional and the result is consistent with the choice for the Gaussian mode. This way it is easy to obtain the electric field from the mode, as the dimensions are given through E_0.The Laguerre-Gauss modes are expressed via the Gauss ones in order to reuse code and together with the choice for N_pm ensures consistent normalization preferences.","category":"page"},{"location":"laguerre-gauss/#","page":"Laguerre-Gauss profile","title":"Laguerre-Gauss profile","text":"note: Note\nIn the expression of the Laguerre-Gauss modes above, the confluent hypergeometric function of the first kind was used instead of the standard associated Laguerre polynomials. This helps reduce the number of factorials required for the normalization constant. The relation between the two is given byL_n^α(z) = frac(α+1)_nn ₁F₁(-nα+1z)where (x)_n is the Pochhammer symbol","category":"page"},{"location":"temporal-profiles/#Temporal-Profiles-1","page":"Temporal Profiles","title":"Temporal Profiles","text":"","category":"section"},{"location":"temporal-profiles/#","page":"Temporal Profiles","title":"Temporal Profiles","text":"The complete expression of the electromaginetic fields describing the laser pulse is given by the combination of the spatial part together with the time dependence.","category":"page"},{"location":"temporal-profiles/#","page":"Temporal Profiles","title":"Temporal Profiles","text":"LaserTypes.g(::Any, ::Any, ::Any)","category":"page"},{"location":"temporal-profiles/#LaserTypes.g-Tuple{Any,Any,Any}","page":"Temporal Profiles","title":"LaserTypes.g","text":"g(z, t, par)\n\nThe time dependence of the fields defined by this package is given by\n\ng(z t) = cos(ω t) envelope(z t)\n\nwhere\n\nz and t are the position on the Oz axis and the time\npar are the laser parameters which pass the corresponding profile to the envelope\n\nand\n\nω is the angular frequency of the laser pulse\nenvelope(z t) is a function that can be used to control the duration of the pulse\n\n\n\n\n\n","category":"method"},{"location":"temporal-profiles/#","page":"Temporal Profiles","title":"Temporal Profiles","text":"There are several types that provide profiles to the temporal envelope:","category":"page"},{"location":"temporal-profiles/#","page":"Temporal Profiles","title":"Temporal Profiles","text":"Modules = [LaserTypes]\nPages = [\"src/envelopes.jl\"]\nOrder = [:type]","category":"page"},{"location":"temporal-profiles/#LaserTypes.ConstantProfile","page":"Temporal Profiles","title":"LaserTypes.ConstantProfile","text":"struct ConstantProfile\n\nThis is the trivial profile\n\nenvelope(z t) = 1\n\nwhich gives an infinite duration for the electromagnetic field configuration (we cannot call this a laser pusle since a pulse implicitely has a finite duration).\n\n\n\n\n\n","category":"type"},{"location":"temporal-profiles/#LaserTypes.GaussProfile","page":"Temporal Profiles","title":"LaserTypes.GaussProfile","text":"GaussProfile{V,T,L}\n\nThis envelope provides a finite duration for the laser pulse and thus can provide a more realistic description of an actual laser pulse.\n\nenvelope(z t) = coshleftleft(fracφτright)^2right\n\nwhere\n\nvarphi = (t - t_0) - fracz - z_0c\n\nand\n\nc is the speed of light\nτ is the duration of the pulse (FWHM) and has the default value 18.02fs\nt₀ is the origin of the time axis and it is 0 by default\nz₀ is the initial position of the intensity peak and has the default value -4*τ*c\n\n\n\n\n\n","category":"type"},{"location":"temporal-profiles/#LaserTypes.QuasiRectangularProfile","page":"Temporal Profiles","title":"LaserTypes.QuasiRectangularProfile","text":"QuasiRectangularProfile{V,T,L}\n\nThis envelope produces a pulse with a predominant constant part of width Δz which could offer better results than the Gaussian profile in the paraxial limit (which is considered for the spatial profiles). The shape of the envelope is given by\n\nenvelope(z t) =\n    begincases\n    coshleftleft(fracφτright)^2right  textfor  φ  0\n    1  textfor  φ  Δz\n    coshleftleft(fracφ - Δzcτright)^2right  textotherwise\n    endcases\n\nwhere\n\nvarphi = (t - t_0) - fracz - z_0c\n\nand\n\nc is the speed of light\nτ is the duration of the pulse (FWHM) and has the default value 18.02fs\nt₀ is the origin of the time axis and it is 0 by default\nz₀ is the initial position of the intensity peak and has the default value -4*τ*c\nΔz is the width of the flat part of the profile and the default value 10*τ*c\n\n\n\n\n\n","category":"type"},{"location":"advanced/#Advanced-usage-1","page":"Advanced usage","title":"Advanced usage","text":"","category":"section"},{"location":"advanced/#","page":"Advanced usage","title":"Advanced usage","text":"If you need more control over the parameters, you can use the constructors for the laser types directly.","category":"page"},{"location":"advanced/#","page":"Advanced usage","title":"Advanced usage","text":"Modules = [LaserTypes]\nPages = [\"src/gauss.jl\", \"src/laguerre-gauss.jl\"]\nOrder = [:type]","category":"page"},{"location":"advanced/#LaserTypes.GaussLaser","page":"Advanced usage","title":"LaserTypes.GaussLaser","text":"struct GaussLaser{V,Q,M,L,F,C,T,P,W,K,E}\n\nThe GaussLaser is defined by the following independent parameters\n\nc is the speed of light in vaccum, with the default value being in SI (c_0 from the CODATA2018 in the PhysicalConstants package)\nq is the electric charge, with the default value being the one for the electron in SI (-e from the CODATA2018 in the PhysicalConstants package)\nm_q is the mass of the charge, with the default value being the one for the electron in SI (m_e from the CODATA2018 in the PhysicalConstants package)\nλ is the laser wavelangth with the default value 0.8μm\na₀ is the normalized vector potential (defined as a_0=fraceAm_e c^2)\nϕ₀ is the initial phase with the default value 0.0\nw₀ is the beam radius at the Rayleigh range or beam waist with the default value 58.0μm\nξx and ξy give the polarization and have the default value 1.0 + 0im and 0.0 + 0im\nprofile is the temporal profile of the pulse and the default one is a Gaussian one\n\nDuring the initialization of the GaussLaser type, some useful derived values are also computed\n\nω is the angular frequency and is given by 2π c  λ\nk is the wavenumber and is given by 2π  λ\nz_R is the Rayleigh range and is given by k w_0^2  2\nT₀ is the laser period and is given by 2π  ω\nE₀ is the amplitude of the electric field and is given by a_0fracm_q c omegaq\n\n\n\n\n\n","category":"type"},{"location":"advanced/#LaserTypes.LaguerreGaussLaser","page":"Advanced usage","title":"LaserTypes.LaguerreGaussLaser","text":"struct LaguerreGaussLaser{V,Q,M,L,F,C,T,P,I,W,K,E,R}\n\nThe LaguerreGaussLaser is defined by the following independent parameters\n\nc is the speed of light in vaccum, with the default value being in SI (c_0 from the CODATA2018 in the PhysicalConstants package)\nq is the electric charge, with the default value being the one for the electron in SI (-e from the CODATA2018 in the PhysicalConstants package)\nm_q is the mass of the charge, with the default value being the one for the electron in SI (m_e from the CODATA2018 in the PhysicalConstants package)\nλ is the laser wavelangth with the default value 0.8μm\na₀ is the normalized vector potential (defined as a_0=fraceAm_e c^2)\nϕ₀ is the initial phase with the default value 0.0\nw₀ is the beam radius at the Rayleigh range or beam waist with the default value 58.0μm\nξx and ξy give the polarization and have the default value 1.0 + 0im and 0.0 + 0im\nprofile is the temporal profile of the pulse and the default one is a Gaussian one\np is the radial index of the mode, p  ℤ p  0, with the default value 1\nm is the azimuthal index of the mode, m  ℤ, with the default value 0\n\nDuring the initialization of the LaguerreGaussLaser type, some useful derived values are also computed\n\nω is the angular frequency and is given by 2π c  λ\nk is the wavenumber and is given by 2π  λ\nz_R is the Rayleigh range and is given by k w_0^2  2\nT₀ is the laser period and is given by 2π  ω\nE₀ is the amplitude of the electric field and is given by a_0fracm_q c omegaq\nNₚₘ is a normalization factor given by sqrt(p+1)_m, with (x)_n the Pochhammer symbol\n\n\n\n\n\n","category":"type"},{"location":"advanced/#","page":"Advanced usage","title":"Advanced usage","text":"note: Note\nWhen specifying derived values, the corresponding independent parameters must also be given in order to avoid discrepancies. For example if you give k, you must also give the corresponding λ.","category":"page"},{"location":"gauss/#gaussian-1","page":"Gaussian profile","title":"Gaussian profile","text":"","category":"section"},{"location":"gauss/#","page":"Gaussian profile","title":"Gaussian profile","text":"The GaussLaser type contains the parameters required for describing a Gaussian laser pulse. The spatial profile of the electric field is given by","category":"page"},{"location":"gauss/#","page":"Gaussian profile","title":"Gaussian profile","text":"newcommandiimathrmi                 imaginary unit\nbeginaligned\n    E_x (rz) = xi_x E_0 fracw_0w expleft-ii k z -fracr^2w^2 -ii left(frack r^22 R - arctgfraczz_R - phi_0right)right\n    E_y (rz) = fracxi_yxi_x E_x(rz) \n    E_z (rz) = frac2 left(ii -fraczz_Rright)kw^2 (z) leftxE_x (rz) +yE_y(rz)right\nendaligned","category":"page"},{"location":"gauss/#","page":"Gaussian profile","title":"Gaussian profile","text":"where","category":"page"},{"location":"gauss/#","page":"Gaussian profile","title":"Gaussian profile","text":"E_0 is the amplitude of the electric field\nξ_x and ξ_y give the polarization (choosing ξ_x=1 and ξ_y=0 for example, would give a linearly polarized field along the x axis while taking them 1sqrt2 and mathrmisqrt2would give right and left-handed circular polarization)\nw(z) is the beam radius at z\nw_0 is the beam waist\nk is the wavenumber\nR is the radius of curvature\narctgfraczz_R is the Gouy phase\nϕ_0 is the initial phase","category":"page"},{"location":"gauss/#","page":"Gaussian profile","title":"Gaussian profile","text":"The magnetic field field is given by","category":"page"},{"location":"gauss/#","page":"Gaussian profile","title":"Gaussian profile","text":"newcommandiimathrmi                 imaginary unit\nbeginaligned\n    B_x(rz) = -frac1cE_y(rz)\n    B_y(rz) = frac1cE_x(rz)\n    B_z (rz) = frac2 left(ii -fraczz_Rright)ckw^2 (z) yE_x (rz) - xE_y(rz)\nendaligned","category":"page"},{"location":"gauss/#","page":"Gaussian profile","title":"Gaussian profile","text":"note: Choice of normalization\nThe electric field is normalized such that in origin (with the default linear polarization) we obtainvecE(000) = (E_000)when the temporal profile is ignored (or ConstantProfile) and thus E_0 is indeed the amplitude of the electric field.","category":"page"},{"location":"#LaserTypes.jl-1","page":"LaserTypes.jl","title":"LaserTypes.jl","text":"","category":"section"},{"location":"#","page":"LaserTypes.jl","title":"LaserTypes.jl","text":"This package aims to provide a common interface for different laser types, such as the Gaussian laser or the Laguerre-Gauss one.","category":"page"},{"location":"#Instalation-1","page":"LaserTypes.jl","title":"Instalation","text":"","category":"section"},{"location":"#","page":"LaserTypes.jl","title":"LaserTypes.jl","text":"To install this package use the Pkg mode in the Julia REPL","category":"page"},{"location":"#","page":"LaserTypes.jl","title":"LaserTypes.jl","text":"]add https://github.com/SebastianM-C/LaserTypes.jl","category":"page"},{"location":"setup/#Laser-initialization-1","page":"Laser initialization","title":"Laser initialization","text":"","category":"section"},{"location":"setup/#","page":"Laser initialization","title":"Laser initialization","text":"Each supported laser has its own (Julia) type which contains all the parameters required to decribe the electromagnetic field configuration. The simplest way to initialize a laser is via the setup_laser function.","category":"page"},{"location":"setup/#","page":"Laser initialization","title":"Laser initialization","text":"Modules = [LaserTypes]\nPages = [\"src/setup.jl\"]","category":"page"},{"location":"setup/#LaserTypes.setup_laser-Tuple{Any,Any}","page":"Laser initialization","title":"LaserTypes.setup_laser","text":"setup_laser(laser, units; τ=nothing, kwargs...)\n\nInitialize the specified laser using the default parameteres in the specified units. Supported units are\n\n:SI: Values are in SI, the numeric types are regular numbers (Float64)\n:SI_unitful: Values are in SI, but using Unitful.jl\n:atomic: Values are in atomic units, the numeric types are regular numbers (Float64)\n:atomic_unitful: Values are in atomic units, but using UnitfulAtomic.jl\n\nThe keyword arguments can be used to give specific values to the parameters instead of using the defaults. You can specify parameteres such as the wavelength and beam waist via λ and w₀. The duration of the pulse (assuming Gaussian temporal profile) can be given via τ. Be sure to use the same units as the ones provided via units. You can also specify dimensionless parameteres such as the normalized vector potential (a₀), the initial phase (ϕ₀) and the polarization (ξx and ξy). See the docomuntation for each laser type for more details on the supported arguments.\n\n\n\n\n\n","category":"method"}]
}

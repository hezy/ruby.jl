using CSV, DataFrames, Plots, LsqFit, Random


function Gaussian(x::Vector{Float64},
                  fwhm::Float64
                  )::Vector{Float64}
    σ = fwhm / (2√(2log(2)))
    return @. 1 / √(2π) / σ * exp(-x^2 / 2σ^2)
end


function Lorentzian(x::Vector{Float64},
                    fwhm::Float64
                    )::Vector{Float64}
    γ = fwhm / 2
    return @. (γ / pi) / (x^2 + γ^2)
end


function Pseudo_Voigt(x::Vector{Float64},
                      fwhm::Float64),
                      n:Float64
                      )::Vector{Float64}
    return n * Lorentzian(x, fwhm) + (1 - n) * Gaussian(x, fwhm)
end


function peak(λ::Vector{Float64},
              λ₀::Float64,
              w::Float64,
              n::Float64,
              A::Float64
              )::Vector{Float64}
    return @. A * Pseudo_Voigt(λ - λ₀, w, n)
end


function background(λ::Vector{Float64},
                    b₀::Float64,
                    b₁::Float64,
                    b₂::Float64
                    )::Vector{Float64}
    return @. b₀ + b₁ * λ + b₂ * λ^2
end


function ruby(λ::Vector{Float64},
              p::Vector{Float64}
              )::Vector{Float64}
    return peak(λ, p[1], p[2], p[3], p[4]) +
           peak(λ, p[5], p[6], p[7], p[8]) +
           background(λ, p[9], p[10], p[11])
end


function plot_it(λ::Vector{Float64},
                 y::Vector{Float64},
                 Title::String
                 )
    default(show = true)
    return plot(λ, y, title = Title, xlabel = "λ (nm)", ylabel = "Intensity (arb.)")
end

function pressure(λ::Float64)::Float64
    A = 1904.0
    B = 7.665
    λ₀ = 693.29427
    return A / B * (((λ - λ₀) / λ₀ + 1)^B - 1)
end

df = DataFrame(CSV.File("6.6GPa-b.csv"))
λ = df[!, 1];
y = df[!, 2];

max_intensity = findmax(df[!, 2])
guess_a₁ = max_intensity[1]
guess_a₂ = 1.6 * guess_a₁
max_position = max_intensity[2]
guess_λ₂ = df[max_position[], 1]
guess_λ₁ = guess_λ₂ - 1.58

p0 = [guess_λ₁, 2, 0.5, guess_a₁, guess_λ₂, 2, 0.5, guess_a₂, 0, 0, 0]

fit = curve_fit(ruby, λ, y, p0)
param = fit.param

scatter(λ, y)
plot!(λ, ruby(λ, param))

p = pressure(param[5])

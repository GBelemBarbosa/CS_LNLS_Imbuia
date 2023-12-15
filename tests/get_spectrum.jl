using StatsBase
include("../src/proximal_subgradient_method.jl")
include("../plots/plots_util.jl")
include("data_treatment.jl") # Do this include the other way around, if you think is better

function get_s(methods:: Vector{Symbol}, iDFT:: Array{<:Number}, interferogram:: Vector{<:Number}; λ=1.0)
    # Your method here
    spectrum = 0

    for mtd ∈ methods
        func     = eval(mtd)
        spectrum = func(iDFT, interferogram; λ=λ)
    end

    return spectrum
end

methods=[:proximal_gradient_method]
λ=1

# This is where you would use your method(s) to get the spectrum from m samples
# I suggest one of you to modify get_s the way you want (to run multiple tests or not, how to output the spectrum(s), etc)
#spectrumₘ=get_s(methods, iDFTₘ, interferogramₘ; λ=λ)

# This is how you would compare both spectrums and interferograms, given the spectrumₘ
# p₂, p₃ = plot_all(frequency_range, spectrumₘ, spectrum, sampleₘ, [i for i=1:n], iDFTₘ, real.(interferogram))
# plot(p₂, p₃)
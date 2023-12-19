include("../src/data_utils.jl")
include("../tests/tests_utils.jl")
include("../plots/plots_util.jl")

percent = 0.05

# Example CSV assuming real interferogram
frequency_range, spectrum, interferogram, n = read_CSV("63148-62-9-IR_Im_Silicone oil.csv")

# If instead all the spectrum (both real and imaginary components) are needed (still assuming real interferogram)
# frequency_range, spectrum, interferogram, n = read_CSV("63148-62-9-IR_Re_Silicone oil.csv", "63148-62-9-IR_Im_Silicone oil.csv")

m = round(Int64, n*percent)

sampleₘ, interferogramₘ, iDFTₘ = get_m_samples(interferogram, m)
# Can change sampleₘ ⊂ (0:n-1) for any x-axis for the interferogram (this information can't be recovered with only the spectrum)
# To do so you need an estimate for the interferogram range (tₘᵢₙ and tₘₐₓ) and perform sampleₘ.*(tₘₐₓ-tₘᵢₙ)/(n-1).+tₘᵢₙ

# spectrumₘ = ...
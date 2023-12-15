include("../src/data_utils.jl")
include("../tests/tests_utils.jl")
include("../plots/plots_util.jl")

percent = 0.05

# Example CSV
frequency_range, spectrum, interferogram, n = read_CSV("63148-62-9-IR_Im_Silicone oil.csv")

m = round(Int64, n*percent)

# Can change 1:n for any x-axis for the interferogram (this information can't be recovered with only the spectrum)
sampleₘ, interferogramₘ, iDFTₘ = get_m_samples(interferogram, m)

# spectrumₘ = ...
using Plots
using LaTeXStrings

function reconstruct_interferogram(:: Vector{<:Number}, t_sample:: Vector{<:Number}, spectrum_sample:: Vector{<:Number}, iDFT_sample:: Array{<:Number}, interferogram:: Vector{<:Number})
    plot(t, real.(interferogram_data), lw =2, alpha = .5, label="Interferogram data")
    scatter!(t_sample, real.(iDFT_sample*spectrum_sample), label="Interferogram from samples")
    plot!(title="Interferogram from irfft of data vs. reconstructed from inverse problem spectrum", xlabel="Optical path [cm]", ylabel="Volts [V]")
end

function plot_spectrum(frequency_range:: Vector{<:Number}, spectrum_sample:: Vector{<:Number}, spectrum_data:: Vector{<:Number})
    plot(frequency_range, real.(spectrum_data), lw =2, alpha = .5, label="Spectrum data (real)")
    plot(frequency_range, imag.(spectrum_data), lw =2, alpha = .5, label="Spectrum data (imaginary)")
    plot!(frequency_range, real.(spectrum_sample), label="Spectrum from samples (real)")
    plot!(frequency_range, imag.(spectrum_sample), label="Spectrum from samples (imaginary)")
    plot!(title="Spectrum data vs. inverse problem spectrum", xlabel="Ratio of "*L"\lambda", ylabel="Fourier frequency term")
end

function plot_all(frequency_range:: Vector{<:Number}, spectrum_sample:: Vector{<:Number}, spectrum:: Vector{<:Number}, t_sample:: Vector{<:Number}, t:: Vector{<:Number}, iDFT_sample:: Array{<:Number}, interferogram:: Vector{<:Number}) 
    p₁, p₂=plot_spectrum(frequency_range, spectrum_sample, spectrum), reconstruct_interferogram(t, t_sample, spectrum_sample, iDFT_sample, interferogram)
end
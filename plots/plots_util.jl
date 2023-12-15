using Plots
using LaTeXStrings

function reconstruct_interferogram(:: Vector{<:Number}, t_sample:: Vector{<:Number}, spectrum_sample:: Vector{<:Number}, iDFT_sample:: Array{<:Number}, interferogram:: Vector{<:Number})
    interferogram_sample = iDFT_sample*spectrum_sample

    plot(t, real.(interferogram), lw =2, alpha = .5, label="Interferogram (real)")
    plot(t, imag.(interferogram), label="Interferogram (imaginary)")
    scatter!(t_sample, real.(interferogram_sample), label="Interferogram from samples (real)")
    scatter!(t_sample, imag.(interferogram_sample), label="Interferogram from samples (imaginary)")
    plot!(title="Interferogram from iDFT of data vs. reconstructed interferogram from samples", xlabel="Optical path [cm]", ylabel="Volts [V]")
end

function plot_spectrum(frequency_range:: Vector{<:Number}, spectrum_sample:: Vector{<:Number}, spectrum_data:: Vector{<:Number})
    plot(frequency_range, real.(spectrum_data), lw =2, alpha = .5, label="Data")
    plot!(frequency_range, real.(spectrum_sample), lw =2, alpha = .5, label="Spectrum from samples (real)")
    plot!(frequency_range, imag.(spectrum_sample), label="Spectrum from samples (imaginary)")
    plot!(title="Spectrum data vs. inverse problem", xlabel="Ratio of "*L"\lambda", ylabel="Fourier frequency term")
end

function plot_all(frequency_range:: Vector{<:Number}, spectrum_sample:: Vector{<:Number}, spectrum:: Vector{<:Number}, t_sample:: Vector{<:Number}, t:: Vector{<:Number}, iDFT_sample:: Array{<:Number}, interferogram:: Vector{<:Number}) 
    p₁, p₂=plot_spectrum(frequency_range, spectrum_sample, spectrum), reconstruct_interferogram(t, t_sample, spectrum_sample, iDFT_sample, interferogram)
end
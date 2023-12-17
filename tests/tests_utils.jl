using FFTW
using LinearAlgebra

function iDFTt(t:: Vector{<:Number}, range:: Tuple{<:Number, <:Number}, n:: Int64)
    # t for the construction must ∈ [0, n-1]
    # but even as t ∈ range, t[begin] and t[end] might not be range[1] and range[2]
    t=Float64.(t)
    t.-=range[1]
    t.*=(n-1)/(range[2]-range[1])

    # if n even, create a new space for the 0 frequency
    if iseven(n)
        n+=1
    end
    n₂=round(Int, n/2)+1
    m=length(t)

    IDFT=zeros(Complex, m, n)
    IDFT[:, n₂].=1.0
    IDFT[:, n₂+1]=exp.((2*π*im/n).*t)
    
    for i=n₂+2:n
        IDFT[:, i]=IDFT[:, i-1].*IDFT[:, n₂+1]
    end
    IDFT[:, begin:n₂-1]=conj(reverse(IDFT[:, n₂+1:end], dims=2))

    return IDFT./n
end

function get_normally_dist_data(range:: Tuple{<:Number, <:Number}, interferogram_f:: Function, n:: Int64; iDFTt=iDFTt, pl=true)
    # n random samples randomly distributed
    sample        = rand(range, n)
    sample.sort()
    interferogram = interferogram_f.(sample) 
    iDFT          = iDFTt(sample, range, n)

    if pl
        scatter(sample,  real.(interferogram), lw =2, alpha = .5, label="Real part")
        scatter!(sample,  imag.(interferogram), label="Imaginary part")
        display(plot!(title="Sampled interferogram data", xlabel="Optical path [cm]", ylabel="Volts [V]"))
    end

    return sample, interferogram, iDFT
end

function get_m_samples(interferogram:: Vector{<:Number}, m:: Int64; iDFTt=iDFTt, pl=true)
    # discrete interferogram assumed to be evenly spaced, as the one outputted from read_CSV
    n=length(interferogram)
    # m random samples sampled from t
    sampleₘ        = sample(1:n, m, ordered=true, replace=false)
    interferogramₘ = interferogram[sampleₘ] 
    iDFTₘ          = iDFTt(sampleₘ, (1, n), n)

    if pl
        scatter(sampleₘ,  real.(interferogramₘ), lw =2, alpha = .5, label="Real part")
        scatter!(sampleₘ,  imag.(interferogramₘ), label="Imaginary part")
        display(plot!(title="m interferogram samples", xlabel="Optical path [cm]", ylabel="Volts [V]"))
    end

    return sampleₘ, interferogramₘ, iDFTₘ
end

function read_CSV(nist_csv_file:: String; pl=true)
    df_nist_data    = CSV.read("data/nist_data/" * nist_csv_file, DataFrame)
    frequency_range = df_nist_data[:, :x]
    spectrum        = df_nist_data[:, :y]

    interferogram = ifft(spectrum)

    if pl
        p₁=plot(frequency_range, spectrum, title="Data spectrum", xlabel="Frequency", ylabel="Fourier component", label="")
        p₂=plot(real.(interferogram), lw =2, alpha = .5, label="Real part")
        plot!(imag.(interferogram), label="Imaginary part")
        plot!(title="Interferogram data", xlabel="Optical path [cm]", ylabel="Volts [V]")

        display(plot(p₁, p₂))
    end

    return frequency_range, spectrum, interferogram, length(interferogram)
end
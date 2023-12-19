using FFTW
using LinearAlgebra

function iDFTt(t:: Vector{Float64}, range:: Tuple{<:Number, <:Number}, n:: Int64)
    # t for the construction must ⊂ [0, n-1]
    # but even as t ⊂ range, t[begin] and t[end] might not be range[1] and range[2]
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

function iDFTt(t:: Vector{Int64}, range:: Tuple{<:Number, <:Number}, n:: Int64)
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

function iDFTt(t:: Vector{<:Number}, n:: Int64) # range assumed to be (0, n-1)
    # t for the construction must ∈ [0, n-1]
    # but even as t ∈ range, t[begin] and t[end] might not be range[1] and range[2]

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

function get_normally_dist_data(range:: Tuple{<:Number, <:Number}, interferogram_f:: Function, m:: Int64, n:: Int64; iDFTt=iDFTt, pl=true)
    # n random samples randomly distributed
    sampleₘ        = rand(range, m)
    sampleₘ.sort()
    interferogramₘ = interferogram_f.(sampleₘ) 
    iDFTₘ          = iDFTt(sampleₘ, range, n)

    if pl
        display(scatter(sampleₘ,  interferogramₘ, title="m interferogram samples", xlabel="Optical path [cm]", ylabel="Volts [V]", label=""))
    end

    return sampleₘ, interferogramₘ, iDFTₘ
end

function get_m_samples(interferogram:: Vector{<:Number}, m:: Int64; iDFTt=iDFTt, pl=true)
    # discrete interferogram assumed to be evenly spaced, as the one outputted from read_CSV
    n              = length(interferogram)
    # m random samples sampled from t
    sampleₘ        = sample(0:n-1, m, ordered=true, replace=false)
    interferogramₘ = interferogram[sampleₘ.+1] 
    iDFTₘ          = iDFTt(sampleₘ, n)

    if pl
        display(scatter(sampleₘ,  interferogramₘ, title="m interferogram samples", xlabel="Optical path [cm]", ylabel="Volts [V]", label=""))
    end

    return sampleₘ, interferogramₘ, iDFTₘ
end

function read_CSV(nist_csv_file:: String; pl=true)
    df_nist_data    = CSV.read("data/nist_data/" * nist_csv_file, DataFrame)
    frequency_range = df_nist_data[:, :x]
    spectrum        = df_nist_data[:, :y]

    n             = 2*length(frequency_range)-1
    interferogram = irfft(spectrum, n)

    if pl
        p₁=plot(frequency_range, spectrum, title="Data spectrum", xlabel="Frequency", ylabel="Fourier component", label="")
        p₂=plot(interferogram, title="Interferogram data", xlabel="Optical path [cm]", ylabel="Volts [V]", label="")

        display(plot(p₁, p₂))
    end

    return frequency_range, spectrum, interferogram, n
end

function read_CSV(real_nist_csv_file:: String, imag_nist_csv_file:: String; pl=true)
    df_nist_data_re = CSV.read("data/nist_data/" * real_nist_csv_file, DataFrame)
    df_nist_data_im = CSV.read("data/nist_data/" * imag_nist_csv_file, DataFrame)
    frequency_range = df_nist_data_re[:, :x]
    spectrum_re     = df_nist_data_re[:, :y]
    spectrum_im     = df_nist_data_im[:, :y]
    spectrum        = spectrum_re+im.*spectrum_im

    n             = 2*length(frequency_range)-1
    interferogram = irfft(spectrum, n)

    if pl
        p₁=plot(frequency_range, spectrum_re, label="Real part")
        plot!(frequency_range, spectrum_im, label="Imaginary part")
        plot!(title="Data spectrum", xlabel="Frequency", ylabel="Fourier component")
        p₂=plot(interferogram, title="Interferogram data", xlabel="Optical path [cm]", ylabel="Volts [V]", label="")

        display(plot(p₁, p₂))
    end

    return frequency_range, spectrum, interferogram, n
end
PDⱼ(xGⱼ:: Number, j:: Int64)=xGⱼ

dDⱼ(xGⱼ:: Number, j:: Int64)=0

function UATPBTAT(x:: Vector{<:Number}, j:: Int64, n:: Int64)
    y=zeros(ComplexF64, n)
    y[j]=x[j]

    return y
end

function UATPBTAT(x:: Vector{<:Number}, T:: Vector{Int64}, n:: Int64)
    y=zeros(ComplexF64, n)
    for j=T
        y[j]=x[j]
    end

    return y
end

ω(x:: Vector{<:Number})=[norm(x[j], 2)^2 for j=1:length(x)]

ωₛ(x:: Vector{<:Number}, s:: Int64; ω=ω)=partialsort(ω(x), s, rev=true)

function Sₛ(x:: Vector{<:Number}, s:: Int64; ω=ω)
    ωx=ω(x)
    uωx=unique(ωx)
    uωxₛ=partialsort(uωx, min(s, length(uωx)), rev=true)

    return findall(x -> x>=uωxₛ, ωx)
end

T(ωx:: Vector{<:Number}, s:: Int64)=partialsortperm(ωx, 1:s, rev=true)

function proxhL(L:: Number, x:: Vector{<:Number}, m:: Int64, s:: Int64, λ:: Number; T=T, ω=ω)
    ωx=ω(x)
    println(ωx[1:3])
    Tωx=T(ωx, s)

    println(searchsortedlast(ωx[Tωx], 2*λ/L, rev=true, lt=<=))
    return UATPBTAT(x, Tωx[1:searchsortedlast(ωx[Tωx], 2*λ/L, rev=true, lt=<=)], m)
end

TL(L:: Number, x:: Vector{<:Number}; ∇f=∇f)=x.-∇f(x)./L
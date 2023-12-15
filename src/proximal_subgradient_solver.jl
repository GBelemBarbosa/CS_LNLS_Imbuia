using Plots
using LaTeXStrings

function proximal_subgradient(f:: Function, ∂f:: Function, Lₖ:: Function, x₀:: Array{<:Number}, L₀:: Number, s:: Number, k_max:: Int64; ϵ=eps(), p=Inf)
    x=x₀
    L=L₀
    k=0
    solved=false
    
    while !(solved || k>=k_max)
        ∂fx=∂f(x)
        println(∂fx[1:3])
        fx=f(x)
        L, x=Lₖ(L, k, x, fx, ∂fx) #Backtracking mais atualização
        
        solved=(norm(∂fx, p)<ϵ)
        k+=1
        println(k)
    end 

    return x
end
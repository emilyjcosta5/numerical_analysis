# Author: Emily Costa

using Calculus

function newtons_method(f::Function, p_0::Number,
                   E::AbstractFloat=1e-5, N_0::Integer=1000)
    p = 0
    for n in N_0
        p = p_0 - (f(p_0) / derivative(f, p_0))
        if (abs(p - p_0) < E)
            break
        end
        p_0 = p
    end
    return p
end

println(newtons_method(x -> x^3 - 2*x^2 - 5, 1))
  
# 2.4 1.a.
println(newtons_method(x -> (x^2) - 2*x*2.7183^(-x) + 2.7183^(-2x), 0.5))

# 2.6 1.a.
println(newtons_method(x -> (x^3) - 2*x^2 - 5, 20))

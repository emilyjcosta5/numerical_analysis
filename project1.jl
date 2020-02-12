using Calculus

function bisection_method(f::Function, a::Number, b::Number,
                   E::AbstractFloat=1e-5, N_0::Integer=1e6)
    fa = f(a)
    p = 0
    for n in N_0
        p = a + (b-a)/2
        fp = f(p)
        if (fp == 0 || fp < E)
            break
        end
        if (fa*fp > 0)
            a = p
            fa = fp
        else b = p
        end
    end
    return p 
end

function newtons_method(f::Function, fp::Function, p_0::Number,
                   E::AbstractFloat=1e-5, N_0::Integer=10)
    p = p_0
    n = 1
    while(n<N_0)&&(abs(p-p_0)>E)
        fpn = f(p_0)
        dfpn = fp(p_0)
        p = p_0 - (fpn / dfpn) 
        println(p)
        n = n + 1
    end
    return p
end

# 1. a.
#println(newtons_method(x -> x^3+x^2+2x, x -> 3x^2+2x+2, 1))
#println(newtons_method(x -> x^3+x^2+2x, x -> 3x^2+2x+2, 10))
println(newtons_method(x -> x^3+x^2+2x, x -> 3x^2+2x+2, 100))
# b. 
#println(newtons_method(x -> exp(x)-x-1, x -> exp(x)-1, 1))

# 2. a. 
# println(bisection_method(x -> x^3+x^2+2x, 0, 5))
# b.
# println(bisection_method(x -> exp(x)-x-1, 0, 5))
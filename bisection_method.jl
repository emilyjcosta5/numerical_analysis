
function bisection_method(f::Function, a::Number, b::Number,
                   E::AbstractFloat=1e-5, N_0::Integer=1000)
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

bisection_method(p -> p^2 - 4,1,2)



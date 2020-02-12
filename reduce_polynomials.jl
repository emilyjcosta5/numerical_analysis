# First use newtons_method then reuce to polynomials of 
# lower degree to determine any complex zeros.

function horners_method(n::Number, a::Array, x_0::Number)
    y = a[0]
    z = a[0]
    i = 1
    for i = 1:length(a)
        y = x_0*y + a[i]
        z = x_0*z + y
    end
    y = x_*y + a[-1]
    return y,z 
end


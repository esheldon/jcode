module testpmap

function myfn(x)
    val=1.0
    for i=1:100000
        val *= 2*x
    end
    val
end

function test(n=100)

    xvals = rand(n)
    res = pmap(myfn, xvals)

    null
end


end # module

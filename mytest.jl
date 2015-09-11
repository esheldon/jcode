module mytest

type Test
    x
    y

    Test() = new(0.0, 0.0)

    function Test(x=0.0, y=0.0)
        println("setting x,y = $(x), $(y)")
        new(x,y)
    end
end

function test()

    println("trying t=Test()")
    t=Test()

    println("trying t2=Test(x=3,y=4)")
    t2=Test(x=3,y=4)

end

end

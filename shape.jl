module shape

export Shape, shape_set_g!

typealias ShapeFloat Float64

type Shape
    g1::ShapeFloat
    g2::ShapeFloat

    e1::ShapeFloat
    e2::ShapeFloat

    eta1::ShapeFloat
    eta2::ShapeFloat
end

# all zeros constructor
Shape() = Shape(0.0,0.0,0.0,0.0,0.0,0.0)

function shape_set_g!(self::Shape, g1::ShapeFloat, g2::ShapeFloat)

    self.g1=g1
    self.g2=g2

    g=sqrt(g1*g1 + g2*g2)

    if g >= 1
        throw(DomainError())
    elseif g==0
        self.e1=0
        self.e2=0
        self.eta1=0
        self.eta2=0
    else
        eta = 2*atanh(g)
        e = tanh(eta)

        if e >= 1
            throw(DomainError())
        end

        cos2theta = g1/g
        sin2theta = g2/g

        self.e1=e*cos2theta
        self.e2=e*sin2theta
        self.eta1=eta*cos2theta
        self.eta2=eta*sin2theta

    end

end

end

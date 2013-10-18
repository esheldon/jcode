module shape

export Shape, ShapeG, shape_set_g!, shape_set_e!, shape_set_eta!, shear

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

function Shape(g1::ShapeFloat, g2::ShapeFloat)
    return ShapeG(g1, g2)
end

function ShapeG(g1::ShapeFloat, g2::ShapeFloat)
    s=Shape()
    shape_set_g!(s, g1, g2)
    return s
end
function ShapeE(e1::ShapeFloat, e2::ShapeFloat)
    s=Shape()
    shape_set_e!(s, e1, e2)
    return s
end
function ShapeEta(eta1::ShapeFloat, eta2::ShapeFloat)
    s=Shape()
    shape_set_eta!(s, eta1, eta2)
    return s
end


(-)(self::Shape) = ShapeG(-self.g1, -self.g2)

# shear comes second
function (+)(self::Shape, sh::Shape)

    sout = Shape()

    oneplusedot = 1.0 + self.e1*sh.e1 + self.e2*sh.e2

    if sh.e1 == 0 && sh.e2 == 0
        # adding nothing
        sout = deepcopy(self)
    elseif oneplusedot != 0
        # result is not round
        esq = sh.e1^2 + sh.e2^2

        fac = (1.0 - sqrt(1.0-esq))/esq

        e1 = self.e1 + sh.e1 + sh.e2*fac*(self.e2*sh.e1 - self.e1*sh.e2)
        e2 = self.e2 + sh.e2 + sh.e1*fac*(self.e1*sh.e2 - self.e2*sh.e1)

        e1 /= oneplusedot
        e2 /= oneplusedot

        shape_set_e!(sout, e1, e2)
    end

    return sout

end


# set the shape given g1,g2 keeping e and eta consistent
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

# set the shape given e1,e2 keeping g and eta consistent
function shape_set_e!(self::Shape, e1::ShapeFloat, e2::ShapeFloat)

    self.e1=e1
    self.e2=e2

    e=sqrt(e1*e1 + e2*e2)

    if e >= 1
        throw(DomainError())
    elseif e==0
        self.g1=0
        self.g2=0
        self.eta1=0
        self.eta2=0
    else

        eta = atanh(e)
        g = tanh(eta/2)

        if g >= 1
            throw(DomainError())
        end

        cos2theta = e1/e
        sin2theta = e2/e

        self.g1=g*cos2theta
        self.g2=g*sin2theta
        self.eta1=eta*cos2theta
        self.eta2=eta*sin2theta

    end

end

# set the shape given eta1,eta2 keeping g and e consistent
function shape_set_eta!(self::Shape, eta1::ShapeFloat, eta2::ShapeFloat)

    self.eta1=eta1
    self.eta2=eta2

    eta=sqrt(eta1*eta1 + eta2*eta2)

    if eta >= 1
        throw(DomainError())
    elseif eta==0
        self.g1=0
        self.g2=0
        self.e1=0
        self.e2=0
    else

        e=tanh(eta)
        g=tanh(eta/2.)

        if g >= 1
            throw(DomainError())
        end
        if e >= 1
            throw(DomainError())
        end

        cos2theta = eta1/eta
        sin2theta = eta2/eta

        self.g1=g*cos2theta
        self.g2=g*sin2theta
        self.e1=e*cos2theta
        self.e2=e*sin2theta

    end

end




end

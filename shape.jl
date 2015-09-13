# shapes
#
# currently the g1,g2  e1,e2 and eta1,eta2 styles are all carried around and
# updated properly when you use the set_g! set_e! or set_eta! functions
#

module shape

import Base.show, Base.IO, Base.string
using mfloat

export Shape

immutable Shape
    g1::MFloat
    g2::MFloat

    e1::MFloat
    e2::MFloat

    eta1::MFloat
    eta2::MFloat

    # all zeros constructor
    Shape() = new(0.0,0.0,0.0,0.0,0.0,0.0)

    # default constructor is for g1,g2
    function Shape(g1::MFloat, g2::MFloat) 

        e1::MFloat = 0.0
        e2::MFloat = 0.0
        eta1::MFloat = 0.0
        eta2::MFloat = 0.0

        g=sqrt(g1*g1 + g2*g2)

        if g >= 1
            throw(DomainError("g >= 1: $(g)"))
        end

        if g != 0
            eta = 2*atanh(g)
            e = tanh(eta)

            if e >= 1
                throw(DomainError("e >= 1: $(e)"))
            end

            cos2theta = g1/g
            sin2theta = g2/g

            e1=e*cos2theta
            e2=e*sin2theta
            eta1=eta*cos2theta
            eta2=eta*sin2theta

        end
        new(g1,g2,e1,e2,eta1,eta2)
    end

end

string(self::Shape) = "g: ($(self.g1), $(self.g2)) e: ($(self.e1), $(self.e2) eta: ($(self.eta1), $(self.eta2))"

show(io::Base.IO, self::Shape) = print(io,string(self))
show(self::Shape) = show(STDOUT, self)

# a new shape explicitly noting the use of g1,g2
ShapeG(g1::MFloat, g2::MFloat) = Shape(g1,g2)

# a new shape explicitly noting the use of e1,e2
function ShapeE(e1::MFloat, e2::MFloat)
    g1::MFloat=0.0
    g2::MFloat=0.0
    eta1::MFloat=0.0
    eta2::MFloat=0.0

    e=sqrt(e1*e1 + e2*e2)

    if e >= 1
        throw(DomainError("e >= 1: $(e)"))
    else

        eta = atanh(e)
        g = tanh(eta/2)

        if g >= 1
            throw(DomainError("g >= 1: $(g)"))
        end

        cos2theta = e1/e
        sin2theta = e2/e

        g1=g*cos2theta
        g2=g*sin2theta
        eta1=eta*cos2theta
        eta2=eta*sin2theta

    end

    Shape(g1,g2,e1,e2,eta1,eta2)
end

# a new shape explicitly noting the use of eta1,eta2
function ShapeEta(eta1::MFloat, eta2::MFloat)
    g1::MFloat=0.0
    g2::MFloat=0.0
    e1::MFloat=0.0
    e2::MFloat=0.0

    eta=sqrt(eta1*eta1 + eta2*eta2)

    if eta >= 1
        throw(DomainError())
    else

        e=tanh(eta)
        g=tanh(eta/2.)

        if g >= 1
            throw(DomainError("g >= 1: $(g)"))
        end
        if e >= 1
            throw(DomainError("e >= 1: $(e)"))
        end

        cos2theta = eta1/eta
        sin2theta = eta2/eta

        g1=g*cos2theta
        g2=g*sin2theta
        e1=e*cos2theta
        e2=e*sin2theta

    end

    Shape(g1,g2,e1,e2,eta1,eta2)
end


# get a new negated version of the shape
Base.(:-)(self::Shape) = ShapeG(-self.g1, -self.g2)

# get a new shape, sheared by the specified amount
function Base.(:+)(self::Shape, shear::Shape)

    g1,g2 = 0.0,0.0

    if shear.g1 != 0 || shear.g2 != 0
        A = 1 + self.g1*shear.g1 + self.g2*shear.g2
        B = self.g2*shear.g1 - self.g1*shear.g2
        denom_inv = 1./(A*A + B*B)

        g1 = A*(self.g1 + shear.g1) + B*(self.g2 + shear.g2)
        g2 = A*(self.g2 + shear.g2) - B*(self.g1 + shear.g1)

        g1 *= denom_inv
        g2 *= denom_inv
    end

    println("g1,g2: ($(g1),$(g2))")
    Shape(g1, g2)
end

# this version uses the distortion formula
function shear_dist(self::Shape, sh::Shape)

    e1, e2 = 0.0, 0.0

    oneplusedot = 1.0 + self.e1*sh.e1 + self.e2*sh.e2

    if (sh.e1 != 0 || sh.e2 != 0) && (oneplusedot != 0)
        esq = sh.e1^2 + sh.e2^2

        fac = (1.0 - sqrt(1.0-esq))/esq

        e1 = self.e1 + sh.e1 + sh.e2*fac*(self.e2*sh.e1 - self.e1*sh.e2)
        e2 = self.e2 + sh.e2 + sh.e1*fac*(self.e1*sh.e2 - self.e2*sh.e1)

        e1 /= oneplusedot
        e2 /= oneplusedot
    end

    ShapeE(e1,e2)
end

function test(;s1=-0.2, s2=-0.1)

    s=Shape(0.2, 0.1)
    println("shape: ",s)

    shear = Shape(s1, s2)
    println("shear: ",shear)

    newshape = s + shear
    println("sheared: ", newshape)


end



end

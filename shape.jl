# shape
# Defines the Shape type and operations
#
#


module shape

import Base.show, Base.IO, Base.string, Base.isequal
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
    end

    if e != 0
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

    if eta != 0
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

#
# representation
#

string(self::Shape) = "g: ($(self.g1), $(self.g2)) e: ($(self.e1), $(self.e2) eta: ($(self.eta1), $(self.eta2))"

show(io::Base.IO, self::Shape) = print(io,string(self))
show(self::Shape) = show(STDOUT, self)


#
# negation and identity
#

Base.(:-)(self::Shape) = ShapeG(-self.g1, -self.g2)
Base.(:+)(self::Shape) = self

#
# shearing
#

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

    Shape(g1, g2)
end

#
# test for equality
#

isequal(s1::Shape, s2::Shape) = (s1.g1==s2.g1) && (s1.g2==s2.g2)

#
# simple test
#

function test(;shear1=-0.2, shear2=-0.1)

    s=Shape(0.2, 0.1)
    println("shape: ",s)

    shear = Shape(shear1, shear2)
    println("shear: ",shear)

    newshape = s + shear
    println("sheared: ", newshape)

    println()
    shear2 = Shape(0.05, 0.05)
    println("shear2: ",shear2)

    newshape2 = s + shear2
    println("sheared: ", newshape2)


    ssame = Shape(0.2, 0.1)
    println("s==ssame:   ",s==ssame)
    println("s==sheared: ",s==newshape)


end



end

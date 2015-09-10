# shapes
#
# currently the g1,g2  e1,e2 and eta1,eta2 styles are all carried around and
# updated properly when you use the set_g! set_e! or set_eta! functions
#

module shape

using mfloat

export Shape
export $

type Shape
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
        self=Shape()
        set_g!(self, g1, g2)
        return self
    end

end

# a representation for the shape
#function Base.(:$)(self::Shape)
#    "g1: $(self.g1) g2: $(self.g2)"
#end

# a new shape explicitly noting the use of g1,g2
function ShapeG(g1::MFloat, g2::MFloat)
    return Shape(g1,g2)
end

# a new shape explicitly noting the use of e1,e2
function ShapeE(e1::MFloat, e2::MFloat)
    s=Shape()
    set_e!(s, e1, e2)
    return s
end

# a new shape explicitly noting the use of eta1,eta2
function ShapeEta(eta1::MFloat, eta2::MFloat)
    s=Shape()
    set_eta!(s, eta1, eta2)
    return s
end


# get a new negated version of the shape
Base.(:-)(self::Shape) = ShapeG(-self.g1, -self.g2)

# get a new shape, sheared by the specified amount
function shear(self::Shape, sh::Shape)

    # zeros
    out = Shape()

    A = 1 + self.g1*sh.g1 + self.g2*sh.g2
    B = self.g2*sh.g1 - self.g1*sh.g2
    denom_inv = 1./(A*A + B*B)

    if sh.g1 != 0 || sh.g2 != 0

        g1o = A*(self.g1 + sh.g1) + B*(self.g2 + sh.g2)
        g2o = A*(self.g2 + sh.g2) - B*(self.g1 + sh.g1)

        g1o *= denom_inv
        g2o *= denom_inv

        set_g!( out, g1o, g2o )
    end

    return out

end

# this version uses the distortion formula
function shear_dist(self::Shape, sh::Shape)

    new_shape = Shape()

    oneplusedot = 1.0 + self.e1*sh.e1 + self.e2*sh.e2

    if (sh.e1 != 0 || sh.e2 != 0) && (oneplusedot != 0)
        esq = sh.e1^2 + sh.e2^2

        fac = (1.0 - sqrt(1.0-esq))/esq

        e1 = self.e1 + sh.e1 + sh.e2*fac*(self.e2*sh.e1 - self.e1*sh.e2)
        e2 = self.e2 + sh.e2 + sh.e1*fac*(self.e1*sh.e2 - self.e2*sh.e1)

        e1 /= oneplusedot
        e2 /= oneplusedot

        set_e!(new_shape, e1, e2)
    end

    return new_shape

end


# set the shape given g1,g2 keeping e and eta consistent
function set_g!(self::Shape, g1::MFloat, g2::MFloat)

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

    return self
end

# set the shape given e1,e2 keeping g and eta consistent
function set_e!(self::Shape, e1::MFloat, e2::MFloat)

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

    return self
end

# set the shape given eta1,eta2 keeping g and e consistent
function set_eta!(self::Shape, eta1::MFloat, eta2::MFloat)

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

    return self
end




end

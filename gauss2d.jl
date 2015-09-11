module gauss2d

import Base.string, Base.show, Base.IO
import mfloat.MFloat

# export means it can be imported, which is about extending it
# remember all the symbols are available to those that "use" the module

export Gauss2D

#
# Gauss2D
#

type Gauss2D
    p::MFloat

    x::MFloat
    y::MFloat

    ixx::MFloat
    ixy::MFloat
    iyy::MFloat

    # derived quantities
    det::MFloat

    dxx::MFloat
    dxy::MFloat
    dyy::MFloat # iyy/det

    norm::MFloat # 1/( 2*pi*sqrt(det) )

    pnorm::MFloat # p*norm

    Gauss2D() = new(0.0,
                    0.0,0.0,
                    0.0,0.0,0.0,
                    0.0,
                    0.0,0.0,0.0,
                    0.0,
                    0.0)

    function Gauss2D(p::MFloat,
                     x::MFloat,
                     y::MFloat,
                     ixx::MFloat,
                     ixy::MFloat,
                     iyy::MFloat)
        self=Gauss2D()

        fill!(self,p,y,x,iyy,ixy,ixx)

        return self
    end

end


function show(io::Base.IO, self::Gauss2D) 
    println("p: $(self.p) x: $(self.x) y: $(self.y) ixx: $(self.ixx) ixy: $(self.ixy) iyy: $(self.iyy)")
end
show(self::Gauss2D) = show(STDOUT, self)


function fill!(self::Gauss2D,
               p::MFloat,
               x::MFloat,
               y::MFloat,
               ixx::MFloat,
               ixy::MFloat,
               iyy::MFloat)

    self.det = iyy*ixx - ixy*ixy;

    if self.det <= 0
        throw(DomainError()) 
    end

    self.p   = p
    self.x = x
    self.y = y
    self.ixx = ixx
    self.ixy = ixy
    self.iyy = iyy

    self.dxx = self.ixx/self.det
    self.dxy = self.ixy/self.det
    self.dyy = self.iyy/self.det
    self.norm = 1./(2*pi*sqrt(self.det))

    self.pnorm = p*self.norm

    return self
end

function eval(self::Gauss2D, x::MFloat, y::MFloat)
    v = x-self.x
    u = y-self.y

    chi2 = self.dxx*u*u + self.dyy*v*v - 2.0*self.dxy*u*v

    val=0.0
    if chi2 < MAX_CHI2
        val = self.pnorm*exp( -0.5*chi2 )
    end

    return val
end



end # module

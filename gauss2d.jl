module gauss2d

export Gauss2D, set!, get

import mfloat.MFloat

# this is so we can add a new method
# note we don't have to export if a method on
# these existing functions
import Base.show, Base.string

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

        set!(self,p,y,x,iyy,ixy,ixx)

        return self
    end

end

function string(s::Gauss2D)
    "p: $(s.p) x: $(s.x) y: $(s.y) ixx: $(s.ixx) ixy: $(s.ixy) iyy: $(s.iyy)"
end

# stdout
function show(s::Gauss2D; stream=STDOUT)
    println(stream,string(s))
end
function set!(self::Gauss2D,
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

function get(self::Gauss2D, x::MFloat, y::MFloat; max_chi2::MFloat = 100.0)
    u = y-self.y
    v = x-self.x

    chi2 = self.dxx*u*u + self.dyy*v*v - 2.0*self.dxy*u*v

    val=0.0
    if chi2 < max_chi2
        val = self.pnorm*exp( -0.5*chi2 )
    end

    return val
end




end

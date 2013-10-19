module gauss2d

export Gauss2D, gauss2d_set!

using mfloat

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

        gauss2d_set!(self,p,y,x,iyy,ixy,ixx)

        return self
    end

end

function gauss2d_set!(self::Gauss2D,
                      p::MFloat,
                      y::MFloat,
                      x::MFloat,
                      iyy::MFloat,
                      ixy::MFloat,
                      ixx::MFloat)

    self.det = iyy*ixx - ixy*ixy;

    if self.det <= 0
        throw(DomainError()) 
    end

    self.p   = p
    self.y = y
    self.x = x
    self.iyy = iyy
    self.ixy = ixy
    self.ixx = ixx

    self.dyy = self.iyy/self.det
    self.dxy = self.ixy/self.det
    self.dxx = self.ixx/self.det
    self.norm = 1./(2*pi*sqrt(self.det))

    self.pnorm = p*self.norm

    return self
end


function evaluate(self::Gauss2D, x::MFloat, y::MFloat)
    u = y-self.y
    v = x-self.x

    chi2 = self.dxx*u*u + self.dyy*v*v - 2.0*self.dxy*u*v

    val = self.pnorm*exp( -0.5*chi2 )

    return val
end

GAUSS2D_MAX_CHI2 = 25.0

# only evaluate within a range
function evaluate(self::Gauss2D, x::MFloat, y::MFloat, max_chi2::MFloat)
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

module gauss2d

import Base.string, Base.show, Base.IO
import mfloat.MFloat

using point2d
using cov2d

# export means it can be imported, which is about extending it
# remember all the symbols are available to those that "use" the module

export Gauss2D

#
# Gauss2D
#

type Gauss2D
    p::MFloat

    cen::Point2D
    cov::Cov2D

    norm::MFloat # 1/( 2*pi*sqrt(cov.det) )
    pnorm::MFloat # p*norm

    Gauss2D() = new(1.0, Point2D(), Cov2D(), 1.0, 1.0)

    function Gauss2D(p::MFloat, cen::Point2D, cov::Cov2D)

        norm = 1./(2*pi*sqrt(cov.det))
        pnorm = p*norm

        new(p, cen, cov, norm, pnorm)
    end

end


function string(self::Gauss2D) 
    "p: $(self.p) cen: $(self.cen) cov: $(self.cov)"
end
function show(io::Base.IO, self::Gauss2D) 
    print(io,string(self))
end
show(self::Gauss2D) = show(STDOUT, self)


function fill!(self::Gauss2D;
               p::MFloat=1.0,
               cen::Point2D=Point2D(),
               cov::Cov2D=Cov2D(),
               )

    self.p   = p
    self.cen = cen
    self.cov = cov

    self.norm = 1./(2*pi*sqrt(self.cov.det))
    self.pnorm = p*self.norm

    return self
end


function eval(self::Gauss2D, pt::Point2D)
    eval(self, x=pt.x, y=pt.y)
end


function eval(self::Gauss2D; x::MFloat=0.0, y::MFloat=0.0)

    cen = self.cen
    cov = self.cov

    xdiff = x-cen.x
    ydiff = y-cen.y

    chi2 =       cov.dxx*ydiff*ydiff
           +     cov.dyy*xdiff*xdiff
           - 2.0*cov.dxy*ydiff*xdiff

    val=0.0
    if chi2 < MAX_CHI2
        val = self.pnorm*exp( -0.5*chi2 )
    end

    val
end



end # module

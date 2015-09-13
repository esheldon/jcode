module igauss2d

import Base.string, Base.show, Base.IO, Base.one
import mfloat.MFloat

using point2d
using cov2d

# export means it can be imported, which is about extending it
# remember all the symbols are available to those that "use" the module

export Gauss2D

#
# Gauss2D
#

immutable Gauss2D
    p::MFloat

    cen::Point2D
    cov::Cov2D

    norm::MFloat # 1/( 2*pi*sqrt(cov.det) )
    pnorm::MFloat # p*norm

    Gauss2D() = Gauss2D(1.0, Point2D(), Cov2D())

    function Gauss2D(p::MFloat, cen::Point2D, cov::Cov2D)

        norm = 1./(2*pi*sqrt(cov.det))
        pnorm = p*norm

        new(p, cen, cov, norm, pnorm)
    end

end

one(Gauss2D) = Gauss2D()

string(self::Gauss2D) = "p: $(self.p) cen: $(self.cen) cov: $(self.cov)"

show(io::Base.IO, self::Gauss2D) = print(io,string(self))

show(self::Gauss2D) = show(STDOUT, self)


function eval(self::Gauss2D, pt::Point2D)
    eval(self, x=pt.x, y=pt.y)
end


function eval(self::Gauss2D; x::MFloat=0.0, y::MFloat=0.0)

    cen = self.cen
    cov = self.cov

    xd = x-cen.x
    yd = y-cen.y

    chi2 =(       cov.dxx*yd*yd
            +     cov.dyy*xd*xd
            - 2.0*cov.dxy*yd*xd )

    val=0.0
    if chi2 < MAX_CHI2
        val = self.pnorm*exp( -0.5*chi2 )
    end

    val
end



end # module

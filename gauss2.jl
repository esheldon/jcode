module gauss2

import Base.string, Base.show, Base.IO, Base.one
import mfloat.MFloat

using point2
using cov2

# export means it can be imported, which is about extending it
# remember all the symbols are available to those that "use" the module

export Gauss2

#
# Gauss2
#

"""
Gauss2(p, cen, cov)

A 2-dimensional gaussian

parameters
----------
- p::`Mfloat`    The total flux in the gaussian
- cen::`Point2`  The center as a `Point2` object
- cov::`Mfloat`  The covariance matrix as a `Cov2` object
"""
immutable Gauss2
    p::MFloat

    cen::Point2
    cov::Cov2

    norm::MFloat # 1/( 2*pi*sqrt(cov.det) )
    pnorm::MFloat # p*norm

    Gauss2() = Gauss2(1.0, Point2(), Cov2())

    function Gauss2(p::MFloat, cen::Point2, cov::Cov2)

        norm = 1./(2*pi*sqrt(cov.det))
        pnorm = p*norm

        new(p, cen, cov, norm, pnorm)
    end

end

one(Gauss2) = Gauss2()

string(self::Gauss2) = "p: $(self.p) cen: $(self.cen) cov: $(self.cov)"

show(io::Base.IO, self::Gauss2) = print(io,string(self))

show(self::Gauss2) = show(STDOUT, self)


"""
val=gauss2.eval(g, y=, x=)

Evaluate a 2-d gaussian at the specified location

parameters
----------
- self::`Gauss2` The 2-d gaussian
- y::`MFloat` **keyword** The y location
- x::`MFloat` **keyword** The x location

returns
-------
- val::`Mfloat` The value at the specified location

Alternatively, you can call as val=gauss2.eval(g, pt)

- pt::`Point2` The location to evaluate the gaussian as a `Point2` object
"""
function eval(self::Gauss2; x::MFloat=0.0, y::MFloat=0.0)

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

function eval(self::Gauss2, pt::Point2)
    eval(self, x=pt.x, y=pt.y)
end




end # module

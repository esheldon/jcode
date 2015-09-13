# symmetric 2d covariance matrix
# the upcoming FixedSizeArray should replace this at some level

module cov2

import Base.string, Base.show, Base.IO
import mfloat.MFloat

export Cov2

#
# Cov2
#

immutable Cov2
    ixx::MFloat
    ixy::MFloat
    iyy::MFloat

    det::MFloat

    dxx::MFloat
    dxy::MFloat
    dyy::MFloat

    function Cov2(; ixx::MFloat=1.0,
                     ixy::MFloat=0.0,
                     iyy::MFloat=1.0)

        det = ixx*iyy - ixy^2
        if det <= 0.0
            throw(DomainError("singular cov, det = $(det)"))
        end

        idet = 1.0/det
        dxx=ixx*idet
        dxy=ixy*idet
        dyy=iyy*idet

        new(ixx,ixy,iyy,det,dxx,dxy,dyy)
    end

end

Base.(:+)(cov1::Cov2, cov2::Cov2) = Cov2(ixx=cov1.ixx+cov2.ixx,
                                            ixy=cov1.ixy+cov2.ixy,
                                            iyy=cov1.iyy+cov2.iyy)

Base.(:-)(cov1::Cov2, cov2::Cov2) = Cov2(ixx=cov1.ixx-cov2.ixx,
                                            ixy=cov1.ixy-cov2.ixy,
                                            iyy=cov1.iyy-cov2.iyy)


string(self::Cov2) = "ixx: $(self.ixx) ixy: $(self.ixy) iyy: $(self.iyy)"

show(io::Base.IO, self::Cov2) = print(io, string(self))

show(self::Cov2) = show(STDOUT, self)


function test()
     cov1=Cov2(ixx=1.5, ixy=0.1, iyy=2.5)
     cov2=Cov2(ixx=1.0, ixy=0.1, iyy=1.0)

     println("cov1:",cov1)
     println("cov2:",cov2)

     println("cov1 + cov2:",cov1 + cov2)
     println("cov1 - cov2:",cov1 - cov2)

     show(cov1)

end

end # module

# symmetric 2d covariance matrix
# the upcoming FixedSizeArray should replace this at some level

module cov2d

import Base.string, Base.show, Base.IO
import mfloat.MFloat

export Cov2D

#
# Cov2D
#

type Cov2D
    ixx::MFloat
    ixy::MFloat
    iyy::MFloat

    det::MFloat

    dxx::MFloat
    dxy::MFloat
    dyy::MFloat # iyy/det

    #Cov2D() = new(1.0,0.0,1.0,
    #              1.0,
    #              1.0,0.0,1.0)

    function Cov2D(; ixx::MFloat=1.0,
                     ixy::MFloat=0.0,
                     iyy::MFloat=1.0)

        det = ixx*iyy - ixy*ixy
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

Base.(:+)(cov1::Cov2D, cov2::Cov2D) = Cov2D(ixx=cov1.ixx+cov2.ixx,
                                            ixy=cov1.ixy+cov2.ixy,
                                            iyy=cov1.iyy+cov2.iyy)

Base.(:-)(cov1::Cov2D, cov2::Cov2D) = Cov2D(ixx=cov1.ixx-cov2.ixx,
                                            ixy=cov1.ixy-cov2.ixy,
                                            iyy=cov1.iyy-cov2.iyy)


function string(self::Cov2D)
    "ixx: $(self.ixx) ixy: $(self.ixy) iyy: $(self.iyy)"
end
function show(io::Base.IO, self::Cov2D) 
    print(io, string(self))
end
show(self::Cov2D) = show(STDOUT, self)


function test()
     cov1=Cov2D(ixx=1.5, ixy=0.1, iyy=2.5)
     cov2=Cov2D(ixx=1.0, ixy=0.1, iyy=1.0)

     println("cov1:",cov1)
     println("cov2:",cov2)

     println("cov1 + cov2:",cov1 + cov2)
     println("cov1 - cov2:",cov1 - cov2)

     show(cov1)

end

end # module

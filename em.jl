
module em

import Base.string, Base.show, Base.IO
import mfloat.MFloat

immutable EMSum
    gi::MFLoat

    # scratch on a given pixel
    trowsum::MFloat
    tcolsum::MFloat
    tu2sum::MFloat
    tuvsum::MFloat
    tv2sum::MFloat

    # sums over all pixels
    pnew::MFloat
    rowsum::MFloat
    colsum::MFloat
    u2sum::MFloat
    uvsum::MFloat
    v2sum::MFloat

    EMSums() = new(0.0,
                   0.0,0.0,0.0,0.0,0.0,
                   0.0,0.0,0.0,0.0,0.0,0.0)
end

function clear_sums!(sums::Vector{EMSum})
    n=length(sums)
    for i=1:n:
        ums[i] = EMSum()
    end
end

function emrun(image::Array{MFloat,2},
               sky_guess::MFloat,
               counts_guess::MFLoat,
               gm::GMix,
               tol::MFloat,
               maxiter::Int)

    skysum::MFLoat=0;
    T::MFLoat=0
    T_last::MFloat=0
    igrat::MFloat=0

    y::Int
    x::Int

    ny, nx = size(image)
    n_points = nx*ny

    # will be the following when we have a jacobian

    #area = n_points*jacob.scale^2
    area = round(MFloat, n_points)
    nsky = sky/counts
    psky = sky/(counts/area)

    ngauss=length(gm)
    sums::Vector{EMSum}

    for i=1:maxiter
        skysum=0.0
        clear_sums!(sums)

        for ix=1:nx
            for iy=1:ny
                ux = convert(MFloat, ix)
                uy = convert(MFloat, iy)


                gtot::MFloat=0.0
                imnorm = image[iy,ix]

                imnorm /= counts

                for i=1:ngauss
                    ydiff = uy - gm[i].y
                    xdiff = ux - gm[i].x

                    y2 = ydiff^2
                    x2 = xdiff^2
                    yx = ydiff*xdiff;

                    chi2 = (      gm[i].cov.dxx*y2
                            +     gm[i].cov.dyy*x2
                            - 2.0*gm[i].cov.dxy*yx )
                    
                    if chi2 < MAX_CHI2 && chi2 > 0.0
                        sums[i].gi = gm[i].pnorm*exp( -0.5*chi2 )
                    else:
                        sums[i].gi = 0.0
                    end

                    gtot += sums[i].gi
                    sums[i].trowsum = y*sums[i].gi
                    sums[i].tcolsum = x*sums[i].gi
                    sums[i].tu2sum  = y2*sums[i].gi
                    sums[i].tuvsum  = yx*sums[i].gi
                    sums[i].tv2sum  = x2*sums[i].gi

                end # loop over gaussians

                igrat = imnorm/gtot

                for i=1:ngauss
                    wtau = sums[i].gi*igrat

                    sums[i].pnew += wtau

                    # row*gi/gtot*imnorm
                    sums[i].rowsum += sums[i].trowsum*igrat
                    sums[i].colsum += sums[i].tcolsum*igrat
                    sums[i].u2sum  += sums[i].tu2sum*igrat
                    sums[i].uvsum  += sums[i].tuvsum*igrat
                    sums[i].v2sum  += sums[i].tv2sum*igrat

                end # sum over gaussians

                skysum += nsky*imnorm/gtot
                #u += jacob->dudcol
                #v += jacob->dvdcol

            end # sum over y
        end # sum over x
    end
end

Base.(:-)(self::Cov2) = Cov2(ixx=-self.ixx,
                             ixy=-self.ixy,
                             iyy=-self.iyy)
Base.(:+)(self::Cov2) = self

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

     println("-cov1:",-cov1)
     println("+cov1:",+cov1)

     println("cov1 + cov2:",cov1 + cov2)
     println("cov1 - cov2:",cov1 - cov2)

     show(cov1)

end

end # module

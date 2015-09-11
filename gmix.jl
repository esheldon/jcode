"""
Defines the following types

GMixPars
GMix

"""
module gmix

# will add methods to these
import Base.length,
       Base.fill!,
       Base.getindex,
       Base.show,
       Base.start,
       Base.next,
       Base.done,
       Base.IO

import mfloat.MFloat

using shape
using gauss2d

# export means it can be imported, which is about extending it
# remember all the symbols are available to those that "use" the module
export GMix, GMixPars, GMixModel

# no enums in Julia yet
const GAUSS=1
const FULL=2
const COELLIP=3
const TURB=4
const EXP=5
const DEV=6

const MAX_CHI2 = 25.0



#
# GMixPars represent parameters that can be used to construct
# a GMix object
#

typealias GMixModel Int

# note FULL not supported here
type GMixPars
    model::GMixModel

    pars::Vector{MFloat}

    shape::Shape

    function GMixPars(model::GMixModel, parsin)
        tpars = convert(Vector{MFloat}, parsin)

        npars=length(tpars)
        nexpected = GMIX_SIMPLE_NPARS[model]
        if npars != nexpected
            throw(error("""
                        for model $model expected $nexpected pars but 
                        got $pars.model
                        """))
        end

        # the current convention is to always have g1,g2 in the
        # 3 4 slots
        # note for FULL this would not make sense, as they can
        # have different shapes
        tshape = Shape(tpars[3], tpars[4])

        new(model, tpars, tshape)
    end

end

# getindex can also be used like gmix_pars[ind]
getindex(self::GMixPars, ind::Int) = self.pars[ind]
length(self::GMixPars) = length(self.pars)

function show(io::Base.IO, self::GMixPars)
    mname=gmix_names[self.model]
    println(io,"model: $(mname)")
    println(io,"pars:  $(self.pars)")
    println(io,"shape: $(self.shape.g1) $(self.shape.g2)")
end

show(self::GMixPars) = show(STDOUT, self)




#
# GMix gaussian mixture
#

type GMix
    data::Vector{Gauss2D}

    GMix() = new(Array(Gauss2D,0))

    function GMix(num::Int) 
        tmp = Array(Gauss2D, num)

        if num > 0
            for i=1:num
                tmp[i] = Gauss2D()
            end
        end

        new(tmp)
    end

end

function simple_zeros(model::GMixModel) 
    """
    all zeros
    """
    ngauss=GMIX_SIMPLE_NGAUSS[model]
    return GMix(ngauss)
end
function make_model(pars::GMixPars)
    """
    The GMixPars fully describe the mixture
    """
    self=simple_zeros(pars.model)

    fill!(self, pars)

    return self
end
function make_model(model::GMixModel, pars)
    """
    the model and pars describe a GMixPars and in turn a GMix
    """
    gp = GMixPars(model, pars)
    return make_model(gp)
end




# from a parameters object
# currently only for simple
function fill!(self::GMix, pars::GMixPars)
    if pars.model==EXP
        fill_exp6!(self, pars)
    elseif pars.model==GAUSS
        fill_gauss!(self, pars)
    else
        throw(error("Bad GMixModel $(pars.model)"))
    end
    return self
end

# indexing and iteration
# getindex can also be used like gmix[ind]
getindex(self::GMix, ind::Int) = self.data[ind]
length(self::GMix) = length(self.data)
start(self::GMix) = start(self.data)
next(self::GMix, i) = next(self.data,i)
done(self::GMix, i) = done(self.data,i)

function show(io::Base.IO, self::GMix)
    for g in self.data
        show(io, g)
    end
end
show(self::GMix) = show(STDOUT,self)

function gmix_eval(self::GMix, x::MFloat, y::MFloat)
    """
    fully inlined the Gauss2D evaluation, enourmously faster
    
    The compiler does not wish to inline; in the future there may be a inline
    annotation for functions

    """
    val::MFloat = 0.0

    for g in self

        v = x-g.x
        u = y-g.y

        chi2 = g.dxx*u*u + g.dyy*v*v - 2.0*g.dxy*u*v

        if chi2 < MAX_CHI2
            val += g.pnorm*exp( -0.5*chi2 )
        end

    end

    return val
end


function get_cen(self::GMix)
    x::MFloat = 0
    y::MFloat = 0
    psum::MFloat = 0

    for g in self
        psum += g.p
        x += g.p*g.x
        y += g.p*g.y
    end

    x /= psum
    y /= psum

    return x,y
end

function set_cen!(self::GMix, x::MFloat, y::MFloat)
    x_cur, y_cur = get_cen(self)

    x_shift = x-x_cur
    y_shift = y-y_cur

    for g in self
        g.x += x_shift
        g.y += y_shift
    end

    return x,y
end

# only makes sense for co-centric gaussians; would
# need to account for different centers
function get_T(self::GMix)
    T::MFloat = 0
    psum::MFloat = 0

    for g in self
        psum += g.p
        T += g.p*(g.ixx + g.iyy)
    end

    T /= psum
    return T
end
function get_psum(self::GMix)
    psum::MFloat = 0

    for g in self
        psum += g.p
    end
    return psum
end
function set_psum!(self::GMix, psum::MFloat)
    psum_cur = get_psum(self)
    rat = psum/psum_cur

    for g in self
        g.p *= rat
    end

    return psum
end

function convolve(obj::GMix, psf::GMix)
    """
    Convolve a gaussian mixture with a psf to get a new mixture
    """
    ntot = length(obj) * length(psf)
    self = GMix(ntot)

    convolve_fill!(self, obj, psf)

    return self
end
function convolve_fill!(self::GMix, obj::GMix, psf::GMix)
    ntot = length(obj) * length(psf)

    sz=length(self)
    if ntot != sz
        throw(error("convolve expected $ntot gauss but have $sz"))
    end

    psf_xcen,psf_ycen = get_cen(psf)

    psum::MFloat = 0
    for g in psf
        psum += g.p
    end

    iself=1
    for og in obj
        for pg in psf

            ixx = og.ixx + pg.ixx
            ixy = og.ixy + pg.ixy
            iyy = og.iyy + pg.iyy

            x = og.x + (pg.x-psf_xcen)
            y = og.y + (pg.y-psf_ycen)

            fill!(self[iself],
                  og.p*pg.p/psum,
                  x,y,ixx,ixy,iyy)
            iself += 1
        end
    end
end




const GMIX_SIMPLE_NPARS = Dict(GAUSS=>6,
                               EXP=>6,
                               DEV=>6,
                               TURB=>8)
const GMIX_SIMPLE_NGAUSS = Dict(GAUSS=>1,
                                EXP=>6,
                                DEV=>10,
                                TURB=>3)

const gmix_names = Dict(GAUSS   => "GAUSS",
                        FULL    => "FULL",
                        COELLIP => "COELLIP",
                        EXP     => "EXP",
                        DEV     => "DEV",
                        TURB    => "TURB")


# from Hogg & Lang, normalized
const fvals_exp6 = MFloat[0.002467115141477932, 
                          0.018147435573256168, 
                          0.07944063151366336, 
                          0.27137669897479122, 
                          0.79782256866993773, 
                          2.1623306025075739]
const pvals_exp6 = MFloat[0.00061601229677880041, 
                          0.0079461395724623237, 
                          0.053280454055540001, 
                          0.21797364640726541, 
                          0.45496740582554868, 
                          0.26521634184240478]

function fill_exp6!(self::GMix, pars::GMixPars)
    const ngauss_expected=6

    if pars.model != EXP
        throw(error("expected model==EXP got $pars.model"))
    end
    if length(self) != ngauss_expected
        throw(error("expected ngauss=$ngauss_expected got $(length(self))"))
    end

    fill_simple!(self, pars, fvals_exp6, pvals_exp6)
end

const fvals_gauss = MFloat[1.0]
const pvals_gauss = MFloat[1.0]

function fill_gauss!(self::GMix, pars::GMixPars)
    const ngauss_expected=1

    if pars.model != GAUSS
        throw(error("expected model==GAUSS got $pars.model"))
    end
    if length(self) != ngauss_expected
        throw(error("expected ngauss=$ngauss_expected got $(length(self))"))
    end

    fill_simple!(self, pars, fvals_gauss, pvals_gauss)
end



function fill_simple!(self::GMix,
                      pars::GMixPars,
                      fvals::Vector{MFloat},
                      pvals::Vector{MFloat})

    x      = pars[1]
    y      = pars[2]
    T      = pars[5]
    counts = pars[6]

    for i=1:length(self)
        gauss=self[i]

        T_i = T*fvals[i]
        counts_i=counts*pvals[i]

        gauss2d.fill!(gauss,
                      counts_i,
                      x,
                      y,
                      (T_i/2.)*(1+pars.shape.e1), # ixx
                      (T_i/2.)*pars.shape.e2,     # ixy
                      (T_i/2.)*(1-pars.shape.e1)) # iyy

    end
    return self

end


#
# rendering the mixture into an image
#

function make_image(self::GMix, dims)
    """
    Get a new image with a rendering of the mixture
    """
    im=zeros(MFloat,dims)
    draw_image!(self, im)
    return im
end

function draw_image!(self::GMix, image::Array{MFloat,2})
    """
    The mixture is *added* to the input image
    """

    nx,ny = size(image)
    for ix=1:nx
        x = convert(MFloat, ix)
        for iy=1:ny
            y = convert(MFloat, iy)
            image[ix,iy] += gmix_eval(self, x, y)
        end
    end
end


function get_loglike(self::GMix,
                     image::Array{MFloat,2},
                     ivar::MFloat)

                 
    loglike::MFloat = 0.0
    s2n_numer::MFloat = 0.0
    s2n_denom::MFloat = 0.0

    nx,ny = size(image)

    for ix=1:nx
        x = convert(MFloat, ix)
        for iy=1:ny
            y = convert(MFloat, iy)

            model_val = gmix_eval(self, x, y)
            pixval = image[ix,iy]

            diff = model_val-pixval
            loglike += diff*diff*ivar
            s2n_numer += pixval*model_val*ivar
            s2n_denom += model_val*model_val*ivar
        end
    end

    loglike *= (-0.5)
    return loglike, s2n_numer, s2n_denom
end

end # module

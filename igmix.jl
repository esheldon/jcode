"""
uses immutable Gauss2D, but GMix is still mutable so
    that it's vector of Gauss2D does not have to be
    copied all the time

Defines the following types

SimplePars
GMix

"""
module igmix

# will add methods to these
import Base.length,
       Base.fill!,
       Base.getindex,
       Base.setindex!,
       Base.string,
       Base.show,
       Base.start,
       Base.next,
       Base.done,
       Base.IO

import mfloat.MFloat

using shape
using igauss2d
using point2d
using cov2d

# export means it can be imported, which is about extending it
# remember all the symbols are available to those that "use" the module
export GMix, SimplePars

# models
# no enums in Julia yet
const GAUSS=1
const FULL=2
const COELLIP=3
const TURB=4
const EXP=5
const DEV=6

const MAX_CHI2 = 25.0



#
# SimplePars represent parameters that can be used to construct
# a GMix object
#

typealias GMixModel Int

type SimplePars

    model::GMixModel
    model_name::String
    pars::Vector{MFloat}

    function SimplePars(model, parsin)
        """
        parameters
        ----------
        model can either be the string or number version,
        e.g. "exp" or EXP
        """

        if !(model in keys(SIMPLE_MODEL_NAMES))
            throw(error("bad simple gmix model: $model"))
        end

        modnum = SIMPLE_MODEL_NUMS[model]
        modname = SIMPLE_MODEL_NAMES[model]

        tpars = convert(Vector{MFloat}, parsin)

        npars=length(tpars)
        nexpected = GMIX_SIMPLE_NPARS[model]
        if npars != nexpected
            throw(error("for model $model expected $nexpected pars " *
                        "but got $npars"))
        end

        new(modnum, modname, tpars)
    end

end

# getindex can also be used like gmix_pars[ind]
getindex(self::SimplePars, ind::Int) = self.pars[ind]
setindex!(self::SimplePars, val::MFloat, ind::Int) = self.pars[ind] = val
length(self::SimplePars) = length(self.pars)

function show(io::Base.IO, self::SimplePars)
    print(io,"model: $(self.model_name)")
    print(io,"pars:  $(self.pars)")
end

show(self::SimplePars) = show(STDOUT, self)




#
# GMix gaussian mixture
#

type GMix
    data::Vector{Gauss2D}

    function GMix(num::Int) 
        tmp = ones(Gauss2D, num)
        new(tmp)
    end
end

# indexing and iteration
# getindex can also be used like gmix[ind]
getindex(self::GMix, ind::Int) = self.data[ind]
setindex!(self::GMix, g::Gauss2D, ind::Int) = self.data[ind] = g
length(self::GMix) = length(self.data)
start(self::GMix) = start(self.data)
next(self::GMix, i) = next(self.data,i)
done(self::GMix, i) = done(self.data,i)

# stringify and showing
function string(self::GMix)
    s=[]
    for i=1:size(self)
        push!(s, string(self[i]) )
    end

    join(s, " ")
end
function string(io::Base.IO, self::GMix)
    print(io, string(self))
end
show(self::GMix) = show(STDOUT,self)

# evaluation and properties
function gmix_eval(self::GMix; x::MFloat=0.0, y::MFloat=0.0)
    """
    Evaluate the gaussian mixture at the specified location
    """
    val::MFloat = 0.0

    for i=1:size(self)
        cen = self[i].cen
        cov = self[i].cov

        xd = x-cen.x
        yd = y-cen.y

        chi2 = (      cov.dxx * yd^2
                +     cov.dyy * xd^2
                - 2.0*cov.dxy * xd*yd )

        if chi2 < MAX_CHI2
            val += self[i].pnorm*exp( -0.5*chi2 )
        end

    end

    val
end

function gmix_eval(self::GMix, pt::Point2D)
    """
    Evaluate the gaussian mixture at the specified location
    """

    gmix_eval(self, x=pt.x, y=pt.y)
end



function get_cen(self::GMix)
    x::MFloat = 0
    y::MFloat = 0
    psum::MFloat = 0

    for i=1:size(self)
        psum += self[i].p
        x += self[i].p*self[i].cen.x
        y += self[i].p*self[i].cen.y
    end

    x /= psum
    y /= psum

    Point2D(x=x, y=y)
end

function set_cen!(self::GMix, pt::Point2D)
    set_cen!(self, x=pt.x, y=pt.y)
end


function set_cen!(self::GMix; x::MFloat=0.0, y::MFloat=0.0)

    pt_cur = get_cen(self)

    x_shift = x-pt_cur.x
    y_shift = y-pt_cur.y

    pshift = Point2D(x=x_shift, y=y_shift)

    for i=1:size(self)
        newcen = self[i].cen + pshift

        self[i] = Gauss2D(self[i].p, newcen, self[i].cov)
    end
end

function get_T(self::GMix)
    """
    use the parallel axis theorem
    """
    psum::MFloat = 0

    ixx::MFloat = 0
    iyy::MFloat = 0

    cen = get_cen(self)

    for i=1:size(self)
        psum += self[i].p

        xdiff = self[i].x - cen.x
        ydiff = self[i].y - cen.y

        ixx += self[i].p*(self[i].ixx + xdiff*xdiff)
        iyy += self[i].p*(self[i].iyy + ydiff*ydiff)

    end

    ixx /= psum
    iyy /= psum

    ixx + iyy
end

function get_psum(self::GMix)
    psum::MFloat = 0

    for i=1:size(self)
        psum += self[i].p
    end
    psum
end
function set_psum!(self::GMix, psum::MFloat)
    psum_cur = get_psum(self)
    rat = psum/psum_cur

    for i=1:size(self)
        self[i] = Gauss2D(self[i].p*rat, self[i].cen, self[i].cov)
    end

    psum
end


# convolutions
function convolve(obj::GMix, psf::GMix)
    """
    Convolve a gaussian mixture with a psf to get a new mixture
    """

    ntot = length(obj) * length(psf)

    self = GMix(ntot)

    convolve_fill!(self, obj, psf)

    self
end

function convolve_fill!(self::GMix, obj::GMix, psf::GMix)
    """
    Convolve two gaussian mixtures, storing in self
    """
    ntot = length(obj) * length(psf)

    sz=length(self)
    if ntot != sz
        throw(error("convolve expected $ntot gauss but have $sz"))
    end

    psf_psum = get_psum(psf)
    psf_cen = get_cen(psf)

    iself=1
    for oi=1:size(obj)
        for pi=1:size(psf)

            cen = obj[oi].cen + (psf[pi].cen - psf_cen)
            cov = obj[oi].cov + psf[pi].cov

            p = obj[oi].p*psf[pi].p/psf_psum

            self[iself] = Gauss2D(p, cen, cov)
            iself += 1
        end
    end
end



# creating models
function simple_zeros(model) 
    """
    make a gmix model, all zeros
    """
    ngauss=GMIX_SIMPLE_NGAUSS[model]
    return GMix(ngauss)
end
function make_simple(pars::SimplePars)
    """
    make a gmix model

    parameters
    ----------
    pars: SimplePars
        A full SimplePars instance
    """
    self=simple_zeros(pars.model)

    fill_simple!(self, pars)

    return self
end
function make_simple(model, pars)
    """
    make a gmix model

    parameters
    ----------
    model: model indicator
        e.g. EXP or "exp"
    pars: array
        e.g. [cen1,cen2,g1,g2,T,flux]
    """
    gp = SimplePars(model, pars)
    return make_simple(gp)
end



# filling existing gmix objects from parameters
function fill_simple!(self::GMix, pars::SimplePars)
    if pars.model==EXP
        fill_exp!(self, pars)
    elseif pars.model==GAUSS
        fill_gauss!(self, pars)
    else
        throw(error("fill_simple does not yet support $(pars.model)"))
    end
end

function fill_exp!(self::GMix, pars::SimplePars)
    const ngauss_expected=6

    if pars.model != EXP
        throw(error("expected model==EXP got $pars.model"))
    end
    if length(self) != ngauss_expected
        throw(error("expected ngauss=$ngauss_expected got $(length(self))"))
    end

    _fill_simple!(self, pars, FVALS_EXP, PVALS_EXP)
end


function fill_gauss!(self::GMix, pars::SimplePars)
    const ngauss_expected=1

    if pars.model != GAUSS
        throw(error("expected model==GAUSS got $pars.model"))
    end
    if length(self) != ngauss_expected
        throw(error("expected ngauss=$ngauss_expected got $(length(self))"))
    end

    _fill_simple!(self, pars, FVALS_GAUSS, PVALS_GAUSS)
end



function _fill_simple!(self::GMix,
                       pars::SimplePars,
                       fvals::Vector{MFloat},
                       pvals::Vector{MFloat})
    """
    fill the gaussian mixture with the input parameters.
    """

    # error checking
    npars=length(pars)
    if npars != 6
        throw(error("simple models should have 6 pars, got $npars"))
    end

    ng,nf,np = length(self), length(fvals), length(pvals)
    if nf != ng || nf != np
        throw(error("gmix has len $ng but fvals and pvals have len $nf,$np"))
    end


    y    = pars[1]
    x    = pars[2]
    g1   = pars[3]
    g2   = pars[4]
    T    = pars[5]
    flux = pars[6]

    sh = Shape(g1, g2)

    for i=1:length(self)

        T_i = T*fvals[i]
        flux_i=flux*pvals[i]

        Thalf=T_i/2.0

        iyy=Thalf*(1-sh.e1)
        ixy=Thalf*sh.e2
        ixx=Thalf*(1+sh.e1)

        cen = Point2D(x=x, y=y)
        cov = Cov2D(iyy=iyy, ixy=ixy, ixx=ixx)

        self[i] = Gauss2D(flux_i, cen, cov)

    end
end




# support string versions for when we read model
# specs from a config file
const GMIX_SIMPLE_NPARS = Dict(
    GAUSS=>6,
    EXP=>6,
    DEV=>6,
    TURB=>6,
    "gauss"=>6,
    "exp"=>6,
    "dev"=>6,
    "turb"=>6
    )

const GMIX_SIMPLE_NGAUSS = Dict(
    GAUSS=>1,
    EXP=>6,
    DEV=>10,
    TURB=>3,
    "gauss"=>1,
    "exp"=>6,
    "dev"=>10,
    "turb"=>3
    )

const SIMPLE_MODEL_NAMES = Dict(
    GAUSS     => "gauss",
    EXP       => "exp",
    DEV       => "dev",
    TURB      => "turb",
    "gauss"   => "gauss",
    "exp"     => "exp",
    "dev"     => "dev",
    "turb"    => "turb"
    )

const SIMPLE_MODEL_NUMS = Dict(
    "gauss"   => GAUSS,
    "exp"     => EXP,
    "dev"     => DEV,
    "turb"    => TURB,
    GAUSS     => GAUSS,
    EXP       => EXP,
    DEV       => DEV,
    TURB      => TURB
    )



# from Hogg & Lang, normalized
const FVALS_EXP = MFloat[0.002467115141477932, 
                         0.018147435573256168, 
                         0.07944063151366336, 
                         0.27137669897479122, 
                         0.79782256866993773, 
                         2.1623306025075739]
const PVALS_EXP = MFloat[0.00061601229677880041, 
                         0.0079461395724623237, 
                         0.053280454055540001, 
                         0.21797364640726541, 
                         0.45496740582554868, 
                         0.26521634184240478]

const FVALS_GAUSS = MFloat[1.0]
const PVALS_GAUSS = MFloat[1.0]


#
# rendering the mixture into an image
#

function make_image(self::GMix, dims)
    """
    Get a new image with a rendering of the mixture
    """
    im=zeros(MFloat,dims)
    draw_image!(self, im)
    im
end

function draw_image!(self::GMix, image::Array{MFloat,2})
    """
    The mixture is *added* to the input image

    note column-major memory layout
    """

    ny,nx = size(image)
    for ix=1:nx
        x = convert(MFloat, ix)

        for iy=1:ny
            y = convert(MFloat, iy)

            image[iy,ix] += gmix_eval(self, x=x, y=y)
        end
    end
end


function get_loglike(self::GMix,
                     image::Array{MFloat,2},
                     ivar::MFloat)
    """
    Calculate the likelihood of the mixture

    note column-major memory layout
    """

                 
    loglike::MFloat = 0.0
    s2n_numer::MFloat = 0.0
    s2n_denom::MFloat = 0.0

    ny,nx = size(image)

    for ix=1:nx
        x = convert(MFloat, ix)
        for iy=1:ny
            y = convert(MFloat, iy)

            model_val = gmix_eval(self, x=x, y=y)
            pixval = image[iy,ix]

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

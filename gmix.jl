"""
2-d Gaussian mixtures for images

GMix
SimplePars
"""
module gmix

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
using gauss2
using point2
using cov2

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

"""
pars=SimplePars(model, pars)

Constructor for SimplePars

parameters
----------
- model::`string` or `number`  e.g. "exp" or gmix.EXP
- pars::`sequence` Should be convertable to a Vector{Mfloat}.

    The layout is [ceny, cenx, g1, g2, T, flux].  The y comes first for
        consistency with 2-dimensionl arrays, which are accessed as image[y, x]

"""
type SimplePars
    model::GMixModel
    model_name::String
    pars::Vector{MFloat}

    function SimplePars(model, parsin)

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

"""
gm=GMix(num)

A gaussian mixture.  Wraps a vector of Gauss2.  You
can construct specific models using gmix.make_simple

parameters
----------
- num::`Int` Number of gaussians in mixture
"""
type GMix
    data::Vector{Gauss2}

    function GMix(num::Int) 
        tmp = ones(Gauss2, num)
        new(tmp)
    end
end

#
# creating models
#

"""
gm=make_simple(spars)

make a gmix model

parameters
----------
- pars::`SimplePars` A full SimplePars instance

returns
-------
- gm::`GMix`  The gaussian mixture


Alternatively, call as make_simple(model, pars)

- model::`model indicator` e.g. gmix.EXP or "exp"
- pars::`sequence` Should be convertable to a Vector{Mfloat}.  The
 layout is [ceny, cenx, g1, g2, T, flux].  The y comes first for
consistency with 2-dimensionl arrays, which are accessed
as image[y, x]
"""

function make_simple(pars::SimplePars)
    ngauss=GMIX_SIMPLE_NGAUSS[pars.model]
    self=GMix(ngauss)

    fill_simple!(self, pars)

    return self
end

function make_simple(model, pars)
    gp = SimplePars(model, pars)
    return make_simple(gp)
end

"""
fill_simple!(self, pars)

Fill an existing GMix object from parameters

parameters
----------
- self::`GMix` The gaussian mixture to fill
- pars::`SimplePars` The full SimplePars object
"""
function fill_simple!(self::GMix, pars::SimplePars)
    if pars.model==EXP
        _fill_exp!(self, pars)
    elseif pars.model==GAUSS
        _fill_gauss!(self, pars)
    else
        throw(error("fill_simple does not yet support $(pars.model)"))
    end
end



# indexing and iteration
# getindex can also be used like gmix[ind]
getindex(self::GMix, ind::Int) = self.data[ind]
setindex!(self::GMix, g::Gauss2, ind::Int) = self.data[ind] = g
length(self::GMix) = length(self.data)
start(self::GMix) = start(self.data)
next(self::GMix, i) = next(self.data,i)
done(self::GMix, i) = done(self.data,i)

# stringify and showing
function string(self::GMix)
    s=[]
    for i=1:length(self)
        push!(s, string(self[i]) )
    end

    join(s, " ")
end
function string(io::Base.IO, self::GMix)
    print(io, string(self))
end
show(self::GMix) = show(STDOUT,self)

# evaluation and properties
"""
val=gmix_eval(gmix, x=, y=, max_chi2=)

Evaluate the gaussian mixture at the specified location

parameters
----------
- self::`GMix` A gaussian mixture
- x::`MFloat` **keyword** The x position at which to evaluate
- y::`MFloat` **keyword** The y position at which to evaluate
- max_chi2::`MFloat`  **keyword, optional** The maximum chi^2 at which to 
evaluate.  Zero is returned for larger values.

returns
------------
- val::`MFloat` The value


Alternatively, can be called as val=gmix_eval(gmix,pt) instead of x=, y=

- pt::`Point2` The position at which to evaluate, as a `Point2` object
"""
function gmix_eval(self::GMix; x::MFloat=0.0, y::MFloat=0.0, max_chi2::MFloat=9.999e9)
    val::MFloat = 0.0

    for i=1:length(self)

        xd = x-self[i].cen.x
        yd = y-self[i].cen.y

        chi2 = (      self[i].cov.dxx * yd^2
                +     self[i].cov.dyy * xd^2
                - 2.0*self[i].cov.dxy * xd*yd )

        if chi2 < max_chi2
            val += self[i].pnorm*exp( -0.5*chi2 )
        end

    end

    val
end

function gmix_eval(self::GMix, pt::Point2, max_chi2::MFloat=9.999e9)
    gmix_eval(self, x=pt.x, y=pt.y, max_chi2=max_chi2)
end


"""
cen=get_cen(gmix)

Get the mean center of the gaussian mixture

parameters
----------
- gmix::`GMix` The gaussian mixture

returns
--------
- cen::`Point2` The new center as a Point2 object
"""
function get_cen(self::GMix)
    x::MFloat = 0
    y::MFloat = 0
    psum::MFloat = 0

    for i=1:length(self)
        psum += self[i].p
        x += self[i].p*self[i].cen.x
        y += self[i].p*self[i].cen.y
    end

    x /= psum
    y /= psum

    Point2(x=x, y=y)
end

"""
set_cen!(gmix, point)

Set the mean center of the gaussian mixture

parameters
----------
- gmix::`GMix` The gaussian mixture
- point::`Point2` The new center as a Point2
"""
function set_cen!(self::GMix, pt::Point2)
    pt_cur = get_cen(self)

    x_shift = pt.x-pt_cur.x
    y_shift = pt.y-pt_cur.y

    pshift = Point2(x=x_shift, y=y_shift)

    for i=1:length(self)
        newcen = self[i].cen + pshift

        self[i] = Gauss2(self[i].p, newcen, self[i].cov)
    end
end

"""
T=get_T(gmix)

Get T = Ixx + Iyy. use the parallel axis theorem for off-center components

parameters
----------
- gmix::`GMix` The gaussian mixture

returns
--------
- T::`MFloat` T=Ixx + Iyy
"""
function get_T(self::GMix)
    psum::MFloat = 0

    ixx::MFloat = 0
    iyy::MFloat = 0

    cen = get_cen(self)

    for i=1:length(self)
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

"""
p=get_psum(gmix)

Get the total flux "psum"

parameters
----------
- gmix::`GMix` The gaussian mixture

returns
--------
- p::`MFloat` The total flux
"""
function get_psum(self::GMix)

    psum::MFloat = 0

    for i=1:length(self)
        psum += self[i].p
    end
    psum
end

"""
set_psum!(gmix)

Set the total flux "psum"

parameters
----------
- gmix::`GMix` The gaussian mixture
- psum::`MFloat` The flux to set

returns
--------
- flux::`MFloat` The new flux
"""
function set_psum!(self::GMix, psum::MFloat)
    psum_cur = get_psum(self)
    rat = psum/psum_cur

    for i=1:length(self)
        self[i] = Gauss2(self[i].p*rat, self[i].cen, self[i].cov)
    end

    psum
end


#
# convolutions
#
"""
gm_new = convolve(gm1,psf)

Convolve a gaussian mixture with a psf to get a new mixture

parameters
----------
- obj::`GMix` The gaussian mixture
- psf::`GMix` The gaussian mixture representing a psf

returns
-------
- new_gmix::`GMix` The convolved GMix object
"""
function convolve(obj::GMix, psf::GMix)

    ntot = length(obj) * length(psf)

    self = GMix(ntot)

    convolve_inplace!(self, obj, psf)

    self
end

"""
convolve_inplace!(self,gm1,psf)

Convolve a gaussian mixture with a psf and store in the 
provided GMix

parameters
----------
- self::`GMix` The gaussian mixture to fill. Must have len=length(obj)*length(psf)
- obj::`GMix` The gaussian mixture
- psf::`GMix` The gaussian mixture representing a psf
"""

function convolve_inplace!(self::GMix, obj::GMix, psf::GMix)
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

            self[iself] = Gauss2(p, cen, cov)
            iself += 1
        end
    end
end

#
# rendering the mixture into an image
#

"""
im=make_image(gmix, dims)

Make a new image, filled with a rendering of the mixture

parameters
----------
- self::`GMix` The gaussian mixture
- dims::`sequence` dimensions of the image

returns
-------
- image::`Array{MFloat,2}`  Tne new image
"""
function make_image(self::GMix, dims)
    im=zeros(MFloat,dims)
    draw_image!(self, im)
    im
end

"""
draw_image!(gmix, image)

Draw the gaussian mixture into the image.  The mixture is *added* to the input
image

parameters
----------
- self::`GMix` The gaussian mixture
- image::`Array{MFloat,2}` The image
"""
function draw_image!(self::GMix, image::Array{MFloat,2})
    ny,nx = size(image)
    for ix=1:nx
        x = convert(MFloat, ix)

        for iy=1:ny
            y = convert(MFloat, iy)

            image[iy,ix] += gmix_eval(self, x=x, y=y)
        end
    end
end

"""
loglike,s2n_num,s2n_denom = get_loglike(gmix, image, weight, max_chi2=)

Calculate the likelihood of the mixture for the input image

parameters
----------
- self::`GMix` The gaussian mixture
- image::`Array{MFloat,2}` The image
- weight::`Array{MFloat,2}` The weight map
- max_chi2::`Array{MFloat,2}` **keyword**  The maximum chi^2
  at which to evaluate the exponential for each gaussian. Default
      is given as gmix.MAX_CHI2

returns
-------
(loglike,s2n_num,s2n_denom)
where

- loglike::`MFloat` The log likelihood
- s2n_num::`MFloat` sum(image*model*weight)
- s2n_denom::`MFloat` sum(model^2*weight)

the s/n can be calculated with s2n_num/sqrt(s2n_denom)
"""

function get_loglike(self::GMix,
                     image::Array{MFloat,2},
                     weight::Array{MFloat,2};
                     max_chi2::MFloat = MAX_CHI2)
                 
    loglike::MFloat = 0.0
    s2n_numer::MFloat = 0.0
    s2n_denom::MFloat = 0.0

    ny,nx = size(image)

    for ix=1:nx
        x = convert(MFloat, ix)
        for iy=1:ny
            y = convert(MFloat, iy)

            ivar = weight[iy,ix]

            if ivar > 0.0
                model_val = gmix_eval(self, x=x, y=y, max_chi2=max_chi2)
                pixval = image[iy,ix]

                diff = model_val-pixval
                loglike += diff*diff*ivar
                s2n_numer += pixval*model_val*ivar
                s2n_denom += model_val*model_val*ivar
            end
        end
    end

    loglike *= (-0.5)
    return loglike, s2n_numer, s2n_denom
end


#
# utility functions and data
#

function _fill_exp!(self::GMix, pars::SimplePars)
    const ngauss_expected=6

    if pars.model != EXP
        throw(error("expected model==EXP got $pars.model"))
    end
    if length(self) != ngauss_expected
        throw(error("expected ngauss=$ngauss_expected got $(length(self))"))
    end

    _fill_simple!(self, pars, FVALS_EXP, PVALS_EXP)
end


function _fill_gauss!(self::GMix, pars::SimplePars)
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

        cen = Point2(x=x, y=y)
        cov = Cov2(iyy=iyy, ixy=ixy, ixx=ixx)

        self[i] = Gauss2(flux_i, cen, cov)

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



function test_make_simple(model; T=16.0, g1=0.1, g2=0.3, flux=100.0)

    sigma=sqrt(T/2.0)
    dim=ceil(Int, 2.0*5.0*sigma)

    dims = (dim,dim)

    c1 = dim/2.0 + 0.01*rand()
    c2 = dim/2.0 + 0.01*rand()
    pars = [c1, c2, g1, g2, T, flux]

    gm = make_simple(model, pars)

    gm
end

function test_make_image(model; T=16.0, g1=0.1, g2=0.3, flux=100.0, show=false)

    gm = test_make_simple(model, T=T, g1=g1, g2=g2, flux=flux)

    sigma=sqrt(T/2.0)
    dim=ceil(Int, 2.0*5.0*sigma)

    dims = (dim,dim)

    image = make_image(gm, dims)

    image
end

end # module

module gmix

# will add methods to these
import Base.length, Base.fill!, Base.string
import Base.show, Base.start, Base.next, Base.done

import mfloat.MFloat
using shape

export GMix,
       GMixModel,
       Gauss2D,
       GMIX_GAUSS,
       GMIX_FULL,
       GMIX_COELLIP,
       GMIX_TURB,
       GMIX_EXP,
       GMIX_DEV,
       GMIX_BD


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
function new_model(pars::GMixPars)
    """
    The GMixPars fully describe the mixture
    """
    self=simple_zeros(pars.model)

    fill!(self, pars)

    return self
end
function new_model(model::GMixModel, pars)
    """
    the model and pars describe a GMixPars and in turn a GMix
    """
    gp = GMixPars(model, pars)
    return new_model(gp)
end




# from a parameters object
# currently only for simple
function fill!(self::GMix, pars::GMixPars)
    if pars.model==GMIX_EXP
        fill_exp6!(self, pars)
    elseif pars.model==GMIX_GAUSS
        fill_gauss!(self, pars)
    else
        throw(error("Bad GMixModel $(pars.model)"))
    end
    return self
end

# indexing and iteration
getindex(self::GMix, ind::Int) = self.data[ind]
length(self::GMix) = length(self.data)
start(self::GMix) = start(self.data)
next(self::GMix, i) = next(self.data,i)
done(self::GMix, i) = done(self.data,i)

function show(self::GMix; stream=STDOUT)
    for g in self.data
        show(g, stream=stream)
    end
end

function eval(self::GMix, x::MFloat, y::MFloat; max_chi2::MFloat = 100.0)
    val::MFloat = 0.0

    for g in self
        val += gauss2d.eval(g, x, y, max_chi2=max_chi2)
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

            gauss2d.set!(self[iself],
                         og.p*pg.p/psum,
                         x,y,ixx,ixy,iyy)
            iself += 1
        end
    end
end



# no enums in Julia yet
const GMIX_GAUSS=1
const GMIX_FULL=2
const GMIX_COELLIP=3
const GMIX_TURB=4
const GMIX_EXP=5
const GMIX_DEV=6
const GMIX_BD=7

const GMIX_SIMPLE_NPARS = {GMIX_GAUSS=>6,
                           GMIX_EXP=>6,
                           GMIX_DEV=>6,
                           GMIX_BD=>8,
                           GMIX_TURB=>8}
const GMIX_SIMPLE_NGAUSS = {GMIX_GAUSS=>1,
                            GMIX_EXP=>6,
                            GMIX_DEV=>10,
                            GMIX_BD=>16,
                            GMIX_TURB=>3}

const gmix_names = {"GMIX_GAUSS"=>GMIX_GAUSS,
                    "GMIX_FULL"=>GMIX_FULL,
                    "GMIX_COELLIP"=>GMIX_COELLIP,
                    "GMIX_EXP"=>GMIX_EXP,
                    "GMIX_DEV"=>GMIX_DEV,
                    "GMIX_BD"=>GMIX_BD,
                    "GMIX_TURB"=>GMIX_TURB}


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

    if pars.model != GMIX_EXP
        throw(error("expected model==GMIX_EXP got $pars.model"))
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

    if pars.model != GMIX_GAUSS
        throw(error("expected model==GMIX_GAUSS got $pars.model"))
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

        gauss2d.set!(gauss,
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

function render(self::GMix, dims::(Int,Int);
                nsub=1, max_chi2::MFloat = 100.0)
    """
    Get a new image with a rendering of the mixture
    """
    im=zeros(MFloat,dims)
    render!(self, im, nsub=nsub, max_chi2=max_chi2)
    return im
end

function render!(self::GMix, image::Array{MFloat,2}; 
                 nsub=1, max_chi2::MFloat = 100.0)
    """
    The mixture is *added* to the input image
    """

    if nsub < 1
        nsub=1
    end

    onebynsub2 = 1./(nsub*nsub)
    stepsize = 1./nsub
    offset = (nsub-1)*stepsize/2.

    nx,ny = size(image)
    for x=1:nx
        for y=1:ny

            tval=0.0
            tx = x-offset
            for ixsub=1:nsub
                ty = y-offset
                for iysub=1:nsub

                    tval += eval(self, tx, ty, max_chi2=max_chi2)

                    ty += stepsize
                end
                tx += stepsize
            end

            tval *= onebynsub2

            image[x,y] += tval

        end
    end
end


#
# Gauss2D
#

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

        set!(self,p,y,x,iyy,ixy,ixx)

        return self
    end

end

function string(s::Gauss2D)
    "p: $(s.p) x: $(s.x) y: $(s.y) ixx: $(s.ixx) ixy: $(s.ixy) iyy: $(s.iyy)"
end

# stdout
function show(s::Gauss2D; stream=STDOUT)
    println(stream,string(s))
end
function set!(self::Gauss2D,
              p::MFloat,
              x::MFloat,
              y::MFloat,
              ixx::MFloat,
              ixy::MFloat,
              iyy::MFloat)

    self.det = iyy*ixx - ixy*ixy;

    if self.det <= 0
        throw(DomainError()) 
    end

    self.p   = p
    self.x = x
    self.y = y
    self.ixx = ixx
    self.ixy = ixy
    self.iyy = iyy

    self.dxx = self.ixx/self.det
    self.dxy = self.ixy/self.det
    self.dyy = self.iyy/self.det
    self.norm = 1./(2*pi*sqrt(self.det))

    self.pnorm = p*self.norm

    return self
end

function eval(self::Gauss2D, x::MFloat, y::MFloat; max_chi2::MFloat = 100.0)
    u = y-self.y
    v = x-self.x

    chi2 = self.dxx*u*u + self.dyy*v*v - 2.0*self.dxy*u*v

    val=0.0
    if chi2 < max_chi2
        val = self.pnorm*exp( -0.5*chi2 )
    end

    return val
end



#
# GMixPars represent parameters that can be used to construct
# a GMix object
#

typealias GMixModel Int

# note GMIX_FULL not supported here
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
        # note for GMIX_FULL this would not make sense, as they can
        # have different shapes
        tshape = Shape(tpars[3], tpars[4])

        new(model, tpars, tshape)
    end

end

getindex(self::GMixPars, ind::Int) = self.pars[ind]
length(self::GMixPars) = length(self.pars)
function show(self::GMixPars; stream=STDOUT)
    println(stream,"model: $(self.model)")
    print(stream,"pars:\n$(self.pars)")
    println(stream,"shape: $(self.shape.g1) $(self.shape.g2)")
end




end # module

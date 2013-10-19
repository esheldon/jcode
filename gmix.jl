module gmix

# will add methods to these
import Base.length, Base.fill!
import Base.show, Base.start, Base.next, Base.done

import mfloat.MFloat
using gauss2d
using shape

export GMix,
       GMixModel,
       gmix_eval,
       gmix_get_cen,
       gmix_set_cen!,
       gmix_get_T,
       gmix_get_psum,
       gmix_set_psum!,
       gmix_simple_zeros,
       gmix_new_model,
       gmix_convolve,
       gmix_convolve_fill!,
       GMIX_GAUSS,
       GMIX_FULL,
       GMIX_COELLIP,
       GMIX_TURB,
       GMIX_EXP,
       GMIX_DEV,
       GMIX_BD


typealias GMixModel Int

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

function gmix_eval(self::GMix, x::MFloat, y::MFloat; max_chi2::MFloat = 100.0)
    val::MFloat = 0.0

    for g in self
        val += gauss2d_eval(g, x, y, max_chi2=max_chi2)
    end

    return val
end

function gmix_get_cen(self::GMix)
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

function gmix_set_cen!(self::GMix, x::MFloat, y::MFloat)
    x_cur, y_cur = gmix_get_cen(self)

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
function gmix_get_T(self::GMix)
    T::MFloat = 0
    psum::MFloat = 0

    for g in self
        psum += g.p
        T += g.p*(g.ixx + g.iyy)
    end

    T /= psum
    return T
end
function gmix_get_psum(self::GMix)
    psum::MFloat = 0

    for g in self
        psum += g.p
    end
    return psum
end
function gmix_set_psum!(self::GMix, psum::MFloat)
    psum_cur = gmix_get_psum(self)
    rat = psum/psum_cur

    for g in self
        g.p *= rat
    end

    return psum
end

function gmix_convolve(obj::GMix, psf::GMix)
    """
    Convolve a gaussian mixture with a psf to get a new mixture
    """
    ntot = length(obj) + length(psf)
    self = GMix(ntot)

    gmix_convolve_fill!(self, obj, psf)

    return self
end
function gmix_convolve_fill!(self::GMix, obj::GMix, psf::GMix)
    ntot = length(obj) + length(psf)

    sz=length(self)
    if ntot != sz
        throw(error("convolve expected $ntot gauss but have $sz"))
    end

    psf_xcen,psf_ycen = gmix_get_cen(psf)

    psum::MFloat = 0
    for g in psf
        psum += g.p
    end

    iself=0
    for og in obj
        for pg in psf

            ixx = og.ixx + pg.ixx
            ixy = og.ixy + pg.ixy
            iyy = og.iyy + pg.iyy

            x = og.x + (pg.x-psf_xcen)
            y = og.y + (pg.y-psf_ycen)

            gauss2d_set(self[iself],
                        og.p*pg.p/psum,
                        x,y,ixx,ixy,iyy)
            iself += 1
        end
    end
end

getindex(self::GMixPars, ind::Int) = self.pars[ind]
length(self::GMixPars) = length(self.pars)
function show(self::GMixPars; stream=STDOUT)
    println(stream,"model: $(self.model)")
    print(stream,"pars:\n$(self.pars)")
    println(stream,"shape: $(self.shape.g1) $(self.shape.g2)")
end


# all zeros
function gmix_simple_zeros(model::GMixModel) 
    ngauss=GMIX_SIMPLE_NGAUSS[model]
    return GMix(ngauss)
end



# from a parameters object
# currently only for simple
function fill!(self::GMix, pars::GMixPars)
    if pars.model==GMIX_EXP
        fill_exp6!(self, pars)
    else
        throw(error("Bad GMixModel $(pars.model)"))
    end
    return self
end

function gmix_new_model(pars::GMixPars)
    self=gmix_simple_zeros(pars.model)

    fill!(self, pars)

    return self
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

        gauss2d_set!(gauss,
                     counts_i,
                     x,
                     y,
                     (T_i/2.)*(1+pars.shape.e1), # ixx
                     (T_i/2.)*pars.shape.e2,     # ixy
                     (T_i/2.)*(1-pars.shape.e1)) # iyy

    end
    return self

end





end # module

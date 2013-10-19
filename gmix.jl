module gmix

# will extend length
import Base.length, Base.fill!
import mfloat.MFloat

export GMix,
       length,
       GMixModel,
       gmix_simple_zeros,
       GMIX_GAUSS,
       GMIX_FULL,
       GMIX_COELLIP,
       GMIX_TURB,
       GMIX_EXP,
       GMIX_DEV,
       GMIX_BD

using gauss2d
using shape

typealias GMixModel Int

type GMix
    data::Vector{Gauss2D}

    GMix() = new(Array(Gauss2D,0))
    function GMix(num::Int) 
        tmp = Array(Gauss2D, num)

        if num > 0
            val = Gauss2D()
            for i=1:num
                tmp[i] = val
            end
        end

        new(tmp)
    end
end

getindex(self::GMix, ind::Int) = self.data[ind]
length(self::GMix) = length(self.data)

# all zeros
function gmix_new_simple(model::GMixModel) 
    ngauss=GMIX_SIMPLE_NGAUSS[model]
    return GMix(ngauss)
end


type GMixPars
    model::GMixModel

    data::Vector{MFloat}

    shape::Shape
end

# from a parameters object
# currently only for simple
function fill!(self::GMix, pars::GMixPars)
    if pars.model==GMIX_EXP
        fill_exp6!(self, pars)
    else
        throw(error("Bad GMixModel $pars.model"))
    end
end

function gmix_new_model(pars::GMixPars)
    self=gmix_new_simple(pars.model)

    fill!(self, pars)
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

    gmix_resize(self, ngauss_expected)

    fill_simple!(self, pars, fvals_exp6, pvals_exp6)
end


function fill_simple!(self::GMix,
                      pars::GMixPars,
                      fvals::Vector{MFloat},
                      pvals::Vector{MFloat})

    x      = pars.data[1]
    y      = pars.data[2]
    T      = pars.data[5]
    counts = pars.data[6]

    for i=1:length(self)
        gauss=self[i]

        T_i = T*Fvals[i]
        counts_i=counts*pvals[i]

        set!(gauss,
             counts_i,
             x,
             y
             (T_i/2.)*(1+pars.shape.e1), # ixx
             (T_i/2.)*pars.shape.e2,     # ixy
             (T_i/2.)*(1-pars.shape.e1)) # iyy

    end
    return self

end





end # module

"""
julia clone of the python package emcee. Only the
affine invariante sampler is implemented

example
-------

    using emcee
    sampler=Sampler(nwalkers, ndim, lnprob_func, x, y, yerr)

    nwalkers=100
    pstart=zeros(nwalkers, ndim)

    # .. fill in your starting points

    burnin=1000
    nstep=1000

    plast,lnplast = sample!(sampler, pstart, burnin)
    p,lnp = sample!(sampler, plast, nstep, lnprob=lnplast)

    chain=flatchain(sampler)
    stats=getstats(chain)
"""
module emcee

import Base.string, Base.show

export Sampler, sample!, flatchain, getstats

type Sampler
    nwalkers::Int        # number of walkers
    ndim::Int            # number of dimensions
    niter::Int           # number of iterations
    lnprobfn::Function   # function to get ln(prob)

    args                 # optional extra arguments for lnprobfn

    a::Float64           # affine stretch parameter optional

    chain::Any                # holds the chain data
    lnprob::Any               # holds lnprob at each step
    naccepted::Any            # number accepted at each step

    function Sampler(nwalkers::Int,
                     ndim::Int,
                     lnprobfn::Function,
                     args...;
                     a=2.0)
        if (nwalkers % 2) != 0
            throw(error("number of walkers must be even, got $nwalkers"))
        end
        if nwalkers <= 2*ndim
            err="""
            The number of walkers needs to be more than twice the
            dimension of your parameter space, got nwalkers=$nwalkers and 
            ndim=$ndim
            """
            throw(error(err))
        end

        nochain=None
        noiter=0
        nolnprob=None
        noaccepted=None

        new(nwalkers,
            ndim,
            noiter,
            lnprobfn,
            args,
            convert(Float64,a),
            nochain,
            nolnprob,
            noaccepted)
    end
end

string(self::Sampler) = """
nwalkers: $(self.nwalkers)
ndim:     $(self.ndim)
a:        $(self.a)
"""
show(io::Base.IO, self::Sampler) = print(io,string(self))
show(self::Sampler) = show(STDOUT, self)

"""
p,lnp = sample!(sampler, pstart, niter, lnprob=)

Sample the posterior.  The internal chain is modified.

parameters
----------
- self::`Sampler` The sampler.
- pstart::Array{Float64,2} The starting position, shape (nwalkers,ndim)
- niter::Int
    Number of iterations to perform
lnprob::`Float64` **keyword, optional** The lnprob values for the input pstart, if not sent will be calculated

returns
-------
(p, lnprob)
- p::`Array{Float64,2}` The last position for each walker
- lnprob:`Vector{Float64}` The last lnprob values for each walker. Can be sent to this function to continue the chain.
"""
function sample!(self::Sampler,
                 pstart::Array{Float64,2}, # (nwalkers,ndim)
                 niter::Int;
                 lnprob=None)
    p = deepcopy(pstart) # current set of parameters
    check_inputs(self, p)
    halfk=fld(self.nwalkers,2)

    self.niter=niter
    self.chain = zeros(self.nwalkers, self.niter, self.ndim)
    self.lnprob = zeros(self.nwalkers, self.niter)
    self.naccepted = zeros(self.nwalkers)
    
    if lnprob == None
        lnprob = get_lnprob_walkers(self, p) # current set of lnprobs
    end

    first = collect(1:halfk)
    second = collect( (halfk+1):self.nwalkers )

    for iter=1:niter
        for (S1, S2) in [(first, second), (second, first)]
            q, newlnp, acc = propose_stretch(self,p[S1,:],p[S2,:],lnprob[S1])

            n=length(S1)
            for i=1:n
                if acc[i]
                    S1i=S1[i]
                    lnprob[S1i] = newlnp[i]
                    for j=1:self.ndim
                        p[S1i,j] = q[i,j]
                    end
                    self.naccepted[S1i] += 1
                end
            end
        end

        self.chain[:,iter,:] = p
        self.lnprob[:,iter] = lnprob
    end

    return p, lnprob
end


"""
chain=flatchain(sampler)

Flatten the chain

parameters
----------
- sampler::`Sampler` The emcee sampler

output
------
- chain::`Array{Float64,2}` The chain, flattened across walkers
"""
function flatchain(self::Sampler)
    dims=size(self.chain)
    flatdims = (dims[1] * dims[2], dims[3])
    return reshape(self.chain, flatdims)
end

"""
mean, cov = getstats(chain)

parameters
----------
- chain::`Array{Float64,2}` The flattened chain

returns
-------
(meanpars,cov)

- mean::`Vector{Float64}` [ndim] array of means
- cov::`Array{Float64,2}` [ndim,ndim] covariance array
"""

function getstats(fchain::Array{Float64,2})

    nstep,ndim=size(fchain)

    meanpars = zeros(ndim)
    cov = zeros( (ndim, ndim) )

    for dim=1:ndim
        meanpars[dim] = mean(fchain[:,dim])
    end

    for i=1:ndim
        idiff = fchain[:,i]-meanpars[i]
        for j=1:ndim
            if i==j
                jdiff = idiff
            else
                jdiff = fchain[:,j]-meanpars[j]
            end

            cov[i,j] = sum(idiff .* jdiff)/(nstep-1.)

            if i != j
                cov[j,i] = cov[i,j]
            end
        end
    end

    return meanpars, cov
end

"""
Get lnprob for each walker in the pars array
"""
function get_lnprob_walkers(self::Sampler,
                            pars::Array{Float64,2})
    nwalkers=size(pars,1)
    lnprob = zeros(nwalkers)
    for i=1:nwalkers
        tpars = vec( pars[i,:] )
        lnprob[i] = self.lnprobfn(tpars, self.args...)
    end

    return lnprob
end


"""
The Goodman and Weare proposal function and acceptance
criteria
"""
function propose_stretch(self::Sampler,
                         p1::Array{Float64,2},
                         p2::Array{Float64,2},
                         lnp1::Vector{Float64})
    nw1 = size(p1,1) # number of walkers
    nw2 = size(p2,1)

    q=zeros( (nw1,self.ndim) )
    newlnprob = zeros(nw1)
    accept = zeros(Bool, nw1)

    for i1=1:nw1
        z = ((self.a - 1.) * rand() + 1)^2 / self.a
        i2 = rand(1:nw2)

        for j=1:self.ndim
            q[i1,j] = p2[i2,j] - z * (p2[i2,j] - p1[i1,j])
        end

        vq = vec(q[i1,:] )
        newlnprob[i1] = self.lnprobfn(vq, self.args...)

        lnpdiff = (self.ndim - 1.) * log(z) + newlnprob[i1] - lnp1[i1]
        accept[i1] = lnpdiff > log(rand())
    end

    return q, newlnprob, accept

end


function check_inputs(self::Sampler, p0::Array{Float64})
    sz=size(p0)
    if length(sz) != 2
        err="""
        p0 must have two dimensions (ndim,nwalkers), got $sz
        """
        throw(error(err))
    end
    if sz[1] != self.nwalkers
        err="""
        p0 must have first dimension of size nwalkers=$(self.nwalkers)
        but got $sz[1]
        """
        throw(error(err))
    end

    if sz[2] != self.ndim
        err="""
        p0 must have second dimension of size ndim=$(self.ndim)
        but got $sz[2]
        """
        throw(error(err))
    end

end


#
# Some tests
#

function line_lnprob_func(pars::Vector{Float64},
                          x::Vector{Float64},
                          y::Vector{Float64},
                          ivar::Float64)
    lnprob::Float64 = 0.0
    tmp::Float64 = 0.0
    off = pars[1]
    slp = pars[2]

    npoints = length(x)
    for i=1:npoints
        tmp = off + slp*x[i] - y[i]
        tmp *= tmp
        lnprob += tmp
    end

    lnprob *= (-0.5*ivar)

    return lnprob
end


function test_line(; ntrial=1)

    nwalkers=20
    burnin=400
    nsteps=100
    ndim=2

    offset=1.0
    slope=2.0
    yerr=0.3
    ivar = 1./(yerr*yerr)
    npoints=25

    x = collect( linspace(0.0, 5.0, npoints) )

    y0 = offset + slope*x

    pstart=zeros( nwalkers, ndim )

    println("true offset: $offset")
    println("true slope:  $slope")
    o_offset=0.0
    o_slope=0.0
    o_offset_ivarsum = 0.0
    o_slope_ivarsum = 0.0

    last_sampler=null
    for i=1:ntrial

        y = y0 + yerr*randn(npoints)
        #args=(x,y,ivar)

        sampler=Sampler(nwalkers, ndim, line_lnprob_func, x, y, ivar)

        # adding offset to guesses
        pstart[:,1] = 0.2 + offset + 0.1*randn(nwalkers)
        pstart[:,2] = 0.1 + slope + 0.1*randn(nwalkers)

        p,lnp=sample!(sampler, pstart, burnin)
        p,lnp=sample!(sampler, p, nsteps, lnprob=lnp)

        fchain=flatchain(sampler)
        means,cov = getstats(fchain)

        offset_meas=means[1]
        offset_err=sqrt(cov[1,1])
        slope_meas=means[2]
        slope_err=sqrt(cov[2,2])

        if ntrial==1
            println("offset: $offset_meas +/- $offset_err")
            println("slope:  $slope_meas +/- $slope_err")
        end

        if ntrial > 1
            o_offset += offset_meas
            o_offset_ivarsum += (1.0/offset_err)^2
            o_slope += slope_meas
            o_slope_ivarsum += (1.0/slope_err)^2
        end

        last_sampler=sampler
    end

    if ntrial > 1
        o_offset /= ntrial
        o_slope /= ntrial

        o_offset_err = sqrt(1.0/o_offset_ivarsum)
        o_slope_err = sqrt(1.0/o_slope_ivarsum)

        println("overall offset: $o_offset +/- $o_offset_err")
        println("overall slope:  $o_slope +/- $o_slope_err")
    end

    last_sampler
end

end

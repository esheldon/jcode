"""
TODO 
    - add reset?
    - add parallel?
        - can split the loops in propose_stretch,
        fill out q first, then do pmap to calculate
        the newlnprob, then another loop to do lndiff
        and accept
     
"""
module emcee

export Sampler, sample!, flatchain, getstats

type Sampler
    nwalkers::Int        # number of walkers
    ndim::Int            # number of dimensions
    niter::Int           # number of iterations
    lnprobfn::Function   # function to get ln(prob)
    args::(Any...)       # extra arguments

    a::Float64           # affine stretch parameter optional

    chain                # holds the chain data
    lnprob               # holds lnprob at each step
    naccepted            # number accepted at each step

    function Sampler(nwalkers::Int,
                     ndim::Int,
                     lnprobfn::Function,
                     args::(Any...);
                     a::Float64=2.0)
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
        new(nwalkers, ndim, noiter, lnprobfn, args, a,
            nochain, nolnprob, noaccepted)
    end
end

function sample!(self::Sampler,
                 pstart::Array{Float64,2}, # (ndim,nwalkers)
                 niter::Int;
                 lnprob=None)
    """
    sample the posterior

    parameters
    ----------
    self::Sampler
        The sampler. The internal chain, lnprop,niter,naccepted
        values are modified.
    pstart::Array{Float64,2}
        The starting position, shape (ndim, nwalkers)
    niter::Int
        Number of iterations to perform
    lnprob: optional 
        The lnprob values for the input pstart, if not
        sent will be calculated

    output
    ------
    p, lnprob
        The last positions and lnprob values from the
        chain.  Can be sent to this function to continue
        the chain.
    """
    p = deepcopy(pstart) # current set of parameters
    check_inputs(self, p)
    halfk=fld(self.nwalkers,2)

    self.niter=niter
    self.chain = zeros(self.ndim, self.niter, self.nwalkers)
    self.lnprob = zeros(self.niter, self.nwalkers)
    self.naccepted = zeros(self.nwalkers)
    
    if lnprob == None
        lnprob = get_lnprob_walkers(self, p) # current set of lnprobs
    end

    first = [1:halfk]
    second = [(halfk+1):self.nwalkers]

    for iter=1:niter
        for (S1, S2) in [(first, second), (second, first)]
            q, newlnp, acc = propose_stretch(self,p[:,S1],p[:,S2],lnprob[S1])

            if any(acc)
                lnprob[S1[acc]] = newlnp[acc]
                p[:,S1[acc]] = q[:,acc]
                self.naccepted[S1[acc]] += 1
            end
        end

        self.chain[:,iter,:] = reshape(p, (self.ndim,1,self.nwalkers))
        self.lnprob[iter,:] = lnprob
    end

    return p, lnprob
end


function propose_stretch(self::Sampler,
                         p1::Array{Float64,2},
                         p2::Array{Float64,2},
                         lnp1::Vector{Float64})
    """
    The Goodman and Weare proposal function
    """
    nw1 = size(p1,2) # number of walkers
    nw2 = size(p2,2)

    q=zeros( (self.ndim, nw1) )
    newlnprob = zeros(nw1)
    accept = zeros(Bool, nw1)

    for i1=1:nw1
        z = ((self.a - 1.) * rand() + 1)^2 / self.a
        i2 = rand(1:nw2)

        q[:,i1] = p2[:,i2] - z * (p2[:,i2] - p1[:,i1])
        newlnprob[i1] = self.lnprobfn(q[:,i1], self.args...)

        lnpdiff = (self.ndim - 1.) * log(z) + newlnprob[i1] - lnp1[i1]
        accept[i1] = lnpdiff > log(rand())
    end

    return q, newlnprob, accept

end


function get_lnprob_walkers(self::Sampler,
                            pars::Array{Float64,2})
    """
    Get lnprob for each walker in the pars array
    """
    nwalkers=size(pars,2)
    lnprob = zeros(nwalkers)
    for i=1:nwalkers
        tpars = pars[:, i] # makes a copy
        lnprob[i] = self.lnprobfn(tpars, self.args...)
    end

    return lnprob
end

function flatchain(self::Sampler)
    """
    Flatten the chain

    parameters
    ----------
    sampler:
        The sampler

    output
    ------
    The flattened chain of shape (ndim, niter*nwalkers)
    """
    dims=size(self.chain)
    flatdims = (dims[1], dims[2]*dims[3])
    return reshape(self.chain, flatdims)
end

function getstats(fchain::Array{Float64,2})
    """
    parameters
    ----------
    fchain: Array{Float64,4}
        The flattened chain
    
    returns
    -------
    meanpars,cov
    """

    (ndim,nstep)=size(fchain)

    meanpars = zeros(ndim)
    cov = zeros( (ndim, ndim) )

    for dim=1:ndim
        meanpars[dim] = mean(fchain[dim,:])
    end

    for i=1:ndim
        idiff = fchain[i,:]-meanpars[i]
        for j=1:ndim
            if i==j
                jdiff = idiff
            else
                jdiff = fchain[j,:]-meanpars[j]
            end

            cov[i,j] = sum(idiff .* jdiff)/(nstep-1.)

            if i != j
                cov[j,i] = cov[i,j]
            end
        end
    end

    return meanpars, cov
end


function check_inputs(self::Sampler, p0::Array{Float64})
    sz=size(p0)
    if length(sz) != 2
        err="""
        p0 must have two dimensions (ndim,nwalkers), got $sz
        """
        throw(error(err))
    end
    if sz[1] != self.ndim
        err="""
        p0 must have first dimension of size ndim=$(self.ndim)
        but got $sz[1]
        """
        throw(error(err))
    end
    if sz[2] != self.nwalkers
        err="""
        p0 must have second dimension of size nwalkers=$(self.nwalkers)
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
    npoints=25*25

    x = linspace(0.0, 5.0, npoints)

    y0 = offset + slope*x

    pstart=zeros( (ndim, nwalkers) )

    println("true offset: $offset")
    println("true slope:  $slope")
    o_offset=0.0
    o_slope=0.0
    o_offset_ivarsum = 0.0
    o_slope_ivarsum = 0.0
    for i=1:ntrial

        y = y0 + yerr*randn(npoints)
        args=(x,y,ivar)

        sampler=Sampler(nwalkers, ndim, line_lnprob_func, args)

        # adding offset to guesses
        pstart[1,:] = 0.2 + offset + 0.1*randn(nwalkers)
        pstart[2,:] = 0.1 + slope + 0.1*randn(nwalkers)

        p,lnp=sample!(sampler, pstart, burnin)
        p,lnp=sample!(sampler, p, nsteps, lnprob=lnp)

        fchain=flatchain(sampler)
        means,cov = getstats(fchain)

        offset_meas=means[1]
        offset_err=sqrt(cov[1,1])
        slope_meas=means[2]
        slope_err=sqrt(cov[2,2])

        println("offset: $offset_meas +/- $offset_err")
        println("slope:  $slope_meas +/- $slope_err")

        if ntrial > 1
            o_offset += offset_meas
            o_offset_ivarsum += (1.0/offset_err)^2
            o_slope += slope_meas
            o_slope_ivarsum += (1.0/slope_err)^2
        end
    end

    if ntrial > 1
        o_offset /= ntrial
        o_slope /= ntrial

        o_offset_err = sqrt(1.0/o_offset_ivarsum)
        o_slope_err = sqrt(1.0/o_slope_ivarsum)

        println("overall offset: $o_offset +/- $o_offset_err")
        println("overall slope:  $o_slope +/- $o_slope_err")
    end
end

end

module emcee

export Sampler, sample!, flatchain, getstats

type Sampler
    nwalkers::Int        # number of walkers
    ndim::Int            # number of dimensions
    niter::Int           # number of iterations
    lnprobfn::Function   # function to get ln(prob)
    a::Float64           # affine stretch parameter
    args                 # extra arguments
    chain                # holds the chain data
    lnprob               # holds lnprob at each step
    naccepted            # number accepted at each step

    function Sampler(nwalkers::Int,
                     ndim::Int,
                     lnprobfn::Function;
                     a::Float64=2.0,
                     args=None)
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
        new(nwalkers, ndim, noiter, lnprobfn, a, args,
            nochain, nolnprob, noaccepted)
    end
end

function sample!(self::Sampler,
                 pstart::Array{Float64,2}, # (ndim,nwalkers)
                 niter::Int;
                 lnprob=None)

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

function flatchain(self::Sampler)
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
function propose_stretch(self::Sampler,
                         p1::Array{Float64,2},
                         p2::Array{Float64,2},
                         lnp1::Vector{Float64})
    
    nw1 = size(p1,2) # number of walkers
    nw2 = size(p2,2)

    q=zeros( (self.ndim, nw1) )
    newlnprob = zeros(nw1)
    accept = zeros(Bool, nw1)

    for i1=1:nw1
        z = ((self.a - 1.) * rand() + 1)^2 / self.a
        i2 = rand(1:nw2)

        q[:,i1] = p2[:,i2] - z * (p2[:,i2] - p1[:,i1])
        newlnprob[i1] = self.lnprobfn(q[:,i1])

        lnpdiff = (self.ndim - 1.) * log(z) + newlnprob[i1] - lnp1[i1]
        accept[i1] = lnpdiff > log(rand())
    end

    return q, newlnprob, accept

end
function propose_stretch_old(self::Sampler,
                         p1::Array{Float64,2},
                         p2::Array{Float64,2},
                         lnp1::Vector{Float64})
    
    nw1 = size(p1,2) # number of walkers
    nw2 = size(p2,2)

    # Generate the vectors of random numbers that will produce the
    # proposal.
    z = ((self.a - 1.) * rand(nw1) + 1) .^2 / self.a
    zz = reshape(z, (1,nw1))

    rint = Array(Int, nw1)
    rand!(1:nw2, rint)

    # Calculate the proposed positions and the log-probability there.
    p2r = p2[:,rint]
    q = p2r .- zz .* (p2r .- p1)
    newlnprob = get_lnprob_walkers(self, q)

    # Decide whether or not the proposals should be accepted.
    lnpdiff = (self.ndim - 1.) * log(z) .+ newlnprob .- lnp1
    accept = lnpdiff > log(rand(length(lnpdiff)))

    return q, newlnprob, accept

end

function get_lnprob_walkers(self::Sampler,
                            pars::Array{Float64,2})
    """

    Get lnprob for each walker in the pars in array (npars, nwalkers) shape.
    Can make this parallel in the future

    """
    nwalkers=size(pars,2)
    lnprob = zeros(nwalkers)
    for i=1:nwalkers
        tpars = pars[:, i] # makes a copy
        lnprob[i] = self.lnprobfn(tpars)
    end

    return lnprob
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

function test_line()

    nwalkers=100
    burnin=400
    nsteps=200
    ndim=2

    offset=1.0
    slope=2.0
    yerr=0.01
    ivar = 1./(yerr*yerr)
    npoints=20

    x = linspace(0.0, 5.0, npoints)

    # closures
    line_func = (pars)->pars[1] + pars[2] * x
    lnprob_func = (pars)-> sum( ivar*(line_func(pars) .- y ).^2 )

    y = line_func([offset,slope])
    y += yerr*randn(npoints)

    sampler=Sampler(nwalkers, ndim, lnprob_func)

    pstart=zeros( (ndim, nwalkers) )
    pstart[1,:] = offset + 0.1*randn(nwalkers)
    pstart[2,:] = slope + 0.1*randn(nwalkers)

    println("pstart")
    print(max(pstart))
    println()
    #for i=1:nwalkers
    #    println(pstart[1,i]," ",pstart[2,i])
    #end

    p,lnp=sample!(sampler, pstart, burnin)
    println("after burnin")
    for i=1:nwalkers
        println(p[1,i]," ",p[2,i]," ",lnp[i])
    end
    #sample!(sampler, p, nsteps, lnprob=lnp)
    p,lnp=sample!(sampler, p, nsteps)

    println("after steps")
    for i=1:nwalkers
        println(p[1,i]," ",p[2,i]," ",lnp[i])
    end
    fchain=flatchain(sampler)
    means,cov = getstats(fchain)

    offset_meas=means[1]
    offset_err=sqrt(cov[1,1])
    slope_meas=means[2]
    slope_err=sqrt(cov[2,2])

    println("true offset: $offset")
    println("true slope:  $slope")
    println("offset: $offset_meas +/- $offset_err")
    println("slope:  $slope_meas +/- $slope_err")
end

end

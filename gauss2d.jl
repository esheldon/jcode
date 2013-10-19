module gauss2d

export Gauss2D, gauss2d_set!

using mfloat

type Gauss2D
    p::MFloat
    col::MFloat
    row::MFloat

    icc::MFloat
    icr::MFloat
    irr::MFloat

    # derived quantities
    det::MFloat

    dcc::MFloat
    dcr::MFloat
    drr::MFloat # irr/det

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
                     col::MFloat,
                     row::MFloat,
                     icc::MFloat,
                     icr::MFloat,
                     irr::MFloat)
        self=Gauss2D()

        gauss2d_set!(self,p,row,col,irr,icr,icc)

        return self
    end

end

function gauss2d_set!(self::Gauss2D,
                      p::MFloat,
                      row::MFloat,
                      col::MFloat,
                      irr::MFloat,
                      icr::MFloat,
                      icc::MFloat)

    self.det = irr*icc - icr*icr;

    if self.det <= 0
        throw(DomainError()) 
    end

    self.p   = p
    self.row = row
    self.col = col
    self.irr = irr
    self.icr = icr
    self.icc = icc

    self.drr = self.irr/self.det
    self.dcr = self.icr/self.det
    self.dcc = self.icc/self.det
    self.norm = 1./(2*pi*sqrt(self.det))

    self.pnorm = p*self.norm

    return self
end



end

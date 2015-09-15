module jacobian

import Base.string, Base.show, Base.IO
import mfloat.MFloat

using point2

export Jacobian, DiagonalJacobian, UnitJacobian


"""
p = Jacobian(cen, jmatrix)

Object representing a jacobian for WCS transformation

parameters
----------
- cen::`Point2` The reference position
- jmatrix::`Array{MFloat,2}` Jacobian matrix

    [ [dudrow, dudcol,
      dvdrow, dvdcol] ]

fields
------
- cen::`Point2` The reference center
- jmatrix::`Jacobian` The jacobian matrix
- det::`MFloat` The jacobian
- scale::`Mfloat` The pixel scale, sqrt(det)
"""

immutable Jacobian
    cen::Point2
    jmatrix::Array{MFloat,2}
    det::MFloat
    scale::MFloat

    function Jacobian(cen::Point2, jmatrix_in)
        
        jmatrix = convert(Array{MFloat,2}, jmatrix_in)
        nrows, ncols = size(jmatrix)
        if nrows != 2 || ncols != 2
            throw(error("jacobian matrix should have shape 2,2, got $nrows,$ncols"))
        end

        det = abs( jmatrix[1,1]*jmatrix[2,2] - jmatrix[1,2]*jmatrix[2,1] )
        scale = sqrt(det)

        new(cen, jmatrix, det, scale)

    end

end

function string(self::Jacobian) 
    """
    cen:
        $(self.cen)
    jmatrix:
     $(self.jmatrix)
    """
end
function show(io::Base.IO, self::Jacobian) 
    print(io, string(self))
end
show(self::Jacobian) = show(STDOUT, self)

"""
j = jacobian.diagonal_jacobian(cen, scale)

Object representing a jacobian for WCS transformation, with
    identity diagonal jacobian matrix.

parameters
----------
- cen::`Point2` The reference position
- scale::`number` The pixel scale
"""
function DiagonalJacobian(cen::Point2, scale)
    jmatrix = zeros(2,2)
    jmatrix[1,1]=scale
    jmatrix[2,2]=scale
    Jacobian(cen, jmatrix)
end


"""
j = jacobian.unit_jacobian(cen)

Object representing a jacobian for WCS transformation, with
    identity diagonal jacobian matrix.

parameters
----------
- cen::`Point2` The reference position
"""
function UnitJacobian(cen::Point2)
    DiagonalJacobian(cen, 1.0)
end


function test()
     pt1=Point2(x=1.5, y=2.3)
     pt2=Point2(x=1.0, y=1.0)

     println("pt1:",pt1)
     println("pt2:",pt2)

     println("-pt1:",-pt1)
     println("+pt1:",+pt1)

     println("pt1 + pt2:",pt1 + pt2)
     println("pt1 - pt2:",pt1 - pt2)

end

end # module

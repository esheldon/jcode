module point2

import Base.string, Base.show, Base.IO
import mfloat.MFloat

export Point2

#
# Point2
#

"""
p = Point2(y=, x=)

Object representing a point in two dimensions

parameters
----------
- y::`MFloat` **keyword** the y position
- x::`MFloat` **keyword** the x position

operations
----------

The + and - operators work as expected, so you can translate points, e.g.

    -pt        Negate the point
    +pt        Just gets a copy of the point
    pt1 + pt2  Add two points to get a new point
    pt1 - pt2  Subtract two points to get a new point
"""
immutable Point2
    x::MFloat
    y::MFloat

    function Point2(;
                     x::MFloat=0.0,
                     y::MFloat=0.0)

        new(x,y)
    end

end
Base.(:-)(self::Point2) = Point2(y=-self.y, x=-self.x)
Base.(:+)(self::Point2) = self
Base.(:+)(pt1::Point2, pt2::Point2) = Point2(x=pt1.x+pt2.x, y=pt1.y+pt2.y)
Base.(:-)(pt1::Point2, pt2::Point2) = Point2(x=pt1.x-pt2.x, y=pt1.y-pt2.y)

function string(self::Point2) 
    "x: $(self.x) y: $(self.y)"
end
function show(io::Base.IO, self::Point2) 
    print(io, string(self))
end
show(self::Point2) = show(STDOUT, self)


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

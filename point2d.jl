module point2d

import Base.string, Base.show, Base.IO
import mfloat.MFloat

export Point2D

#
# Point2D
#

type Point2D
    x::MFloat
    y::MFloat

    #Point2D() = new(0.0,0.0)

    function Point2D(;
                     x::MFloat=0.0,
                     y::MFloat=0.0)

        new(x,y)
    end

end

Base.(:+)(pt1::Point2D, pt2::Point2D) = Point2D(x=pt1.x+pt2.x, y=pt1.y+pt2.y)
Base.(:-)(pt1::Point2D, pt2::Point2D) = Point2D(x=pt1.x-pt2.x, y=pt1.y-pt2.y)

function string(self::Point2D) 
    "x: $(self.x) y: $(self.y)"
end
function show(io::Base.IO, self::Point2D) 
    print(io, string(self))
end
show(self::Point2D) = show(STDOUT, self)


function test()
     pt1=Point2D(x=1.5, y=2.3)
     pt2=Point2D(x=1.0, y=1.0)

     println("pt1:",pt1)
     println("pt2:",pt2)

     println("pt1 + pt2:",pt1 + pt2)
     println("pt1 - pt2:",pt1 - pt2)

end

end # module

module point2

import Base.string, Base.show, Base.IO
import mfloat.MFloat

export Point2

#
# Point2
#

immutable Point2
    x::MFloat
    y::MFloat

    function Point2(;
                     x::MFloat=0.0,
                     y::MFloat=0.0)

        new(x,y)
    end

end

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

     println("pt1 + pt2:",pt1 + pt2)
     println("pt1 - pt2:",pt1 - pt2)

end

end # module

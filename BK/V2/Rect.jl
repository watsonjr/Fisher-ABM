module Rect
#Implementation of a 2D axis-aligned rectangle
# Ported from COS 226 Java code from Princeton University
# For numerical purposes: rectangle assumed to be within R2 with float arguments
# NOTE: methods comparing rectangles need to be passed the objects in question.
export Rectangle, width, height, intersects, distanceTo, distanceSquaredTo, contains


type Rectangle
        xmin::Float64
        xmax::Float64
        ymin::Float64
        ymax::Float64
end

function width(r::Rectangle)
        return r.xmax - r.xmin;
end

function height(r::Rectangle)
        return r.ymax - r.ymin;
end

#Check if this Rectangle intersects that Rectangle. Returns boolean. 
function intersects(r::Rectangle,that::Rectangle)
        return r.xmax >= that.xmin && r.ymax >= that.ymin && that.xmax >= r.xmin && that.ymax >= r.ymin;
end

#Return distance from p to closest point on this Rectangle
function distanceTo(r::Rectangle,p::Array{Float64,2})
        return sqrt(distanceSquaredTo(r,p));
end

#Return squared distance from p to closest point on this rectangle
function distanceSquaredTo(r::Rectangle,p::Array{Float64,2})

        if p[1] < r.xmin
                dx = p[1] - r.xmin;
        elseif p[1] > r.xmin
                dx = p[1] - r.xmax;
        else
        end
        if p[2] < r.ymin
                dy = p[2] - r.ymin;
        elseif p[2] < r.ymax
                dx = p[2] - r.ymax;
        else
        end

	return dx*dx + dy*dy; 
end

# Does this rectangle contain p? 
function contains(r::Rectangle,p::Array{Float64,2})
        return p[1] >= r.xmin && (p[1] <= r.xmax) && (p[2] >= r.ymin) && (p[2] <= r.ymax);
end

end

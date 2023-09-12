draw_default_config = Dict(
    :output => "luxor-drawing-$(Dates.format(Dates.now(), "HHMMSS_s")).png",
    :scale => 100.0,
    :line_width => 2.0,
    :line_color => [0, 0, 0],
    :line_color_positive_feedback => [0.9, 0, 0],
    :line_color_negative_feedback => [0, 0, 0.9],
    :line_style => "dashed",
    :face_color => [0.25, 0.25, 0.25, 0.5],
    :face_color_coherent => [0.0, 0.0, 0.8, 0.5],
    :face_color_incoherent => [0.8, 0.0, 0.0, 0.5],
    :relative_canvas_size => 3.5,
    :relative_node_radius => 0.25,
    :relative_path_offset => 0.05,
    :relative_space_between_arrows => 0.1,
    :relative_arc_offset => 0.375,
    :relative_arc_radius => 0.2,
    :relative_arrowhead_width => 0.15,
    :arc_start_angle => 5π / 8,
    :arc_end_angle => 11π / 16,
    :arrowhead_angle => π / 6
)
# Parameters named relative_ are defined relative to the radius of the unit circle where nodes lie


"""
   get_polygon_vertex_coordinates(n::Int)

   Returns vertex coordinates of an n-gon inscribed in the unit circle. The first node is placed in the 12 o'clock direction and the nodes are indexed clockwise. Note that the definition of axes follows that of svg's.

# Arguments (Required)
- `n::Int`: Number of vertices

# Returns
- `x::AbstractArray`: Vector of x coordinates
- `y::AbstractArray`: Vector of y coordinates
"""
function get_polygon_vertex_coordinates(n::Int)
    angle = [-1 * π / 2 + 2π * (i - 1) / n for i in 1:n]
    # Starting from upper center, run clockwise (following svg angle definition)

    x = [real(exp(i * im)) for i in angle]
    y = [imag(exp(i * im)) for i in angle]
    # Also follows svg definition of y axis

    return x, y
end


"""
   draw_node(type::Int, x::Float64, y::Float64, draw_config::Dict)

   Draws a node.

# Arguments (Required)
- `type::Int`: Type of node. 1 if the node gets coherent inputs and -1 if it gets incoherent inputs
- `x::Float64`: x coordinate of node
- `y::Float64`: y coordinate of node

# Arguments (Optional)
- `draw_config::Dict`: Specification of drawing styles
"""
function draw_node(type::Int, x::Float64, y::Float64; draw_config::Dict = draw_default_config)
    # Unpack drawing settings
    scale = draw_config[:scale]
    relative_node_radius = draw_config[:relative_node_radius]
    line_color = draw_config[:line_color]
    line_width = draw_config[:line_width]
    if type == 1
        # Coherent input
        face_color = draw_config[:face_color_coherent]
    elseif type == -1
        # Incoherent input
        face_color = draw_config[:face_color_incoherent]
    else
        face_color = draw_config[:face_color]
    end

    # Fill face
    setcolor(face_color...)
    circle(Point(scale * x, scale * y), relative_node_radius * scale, action=:fill)

    # Draw border
    setcolor(line_color...)
    setline(line_width)
    circle(Point(scale * x, scale * y), relative_node_radius * scale, action=:stroke)
end


"""
   draw_edge(type::Int, x1::Float64, y1::Float64, draw_config::Dict)

   Draws an edge (self-loop). Note that this function is overloaded.

# Arguments (Required)
- `type::Int`: Type of interaction. 1 for positive interaction and -1 for negative interaction
- `x::Float64`: x coordinate of node
- `y::Float64`: y coordinate of node

# Arguments (Optional)
- `draw_config::Dict`: Specification of drawing styles
"""
function draw_edge(type::Int, x1::Float64, y1::Float64; draw_config::Dict = draw_default_config)
    # Unpack drawing settings
    scale = draw_config[:scale]
    line_width = draw_config[:line_width]
    relative_arc_offset = draw_config[:relative_arc_offset]
    relative_arc_radius = draw_config[:relative_arc_radius]
    arc_start_angle = draw_config[:arc_start_angle]
    arc_end_angle = draw_config[:arc_end_angle]
    relative_arrowhead_width = draw_config[:relative_arrowhead_width]
    arrowhead_angle = draw_config[:arrowhead_angle]

    # Get the position and the radius of the arc
    centerpos = Point(scale * x1 * (1 + relative_arc_offset), scale * y1 * (1 + relative_arc_offset))
    radius = scale * relative_arc_radius

    # Get the angular span of the arc
    startangle = atan(y1, x1) - arc_start_angle
    endangle = atan(y1, x1) + arc_end_angle

    arrowheadlength = tan(π / 2 - arrowhead_angle) * relative_arrowhead_width * scale / 2

    if type == 1
        # Positive input
        arrow(centerpos, radius, startangle, endangle, arrowheadlength=arrowheadlength, arrowheadangle=arrowhead_angle, linewidth=line_width)
    elseif type == -1
        # Negative input
        arrow(centerpos, radius, startangle, endangle, arrowheadlength=0, linewidth=line_width)
        arrow(centerpos, radius, startangle - eps(Float64), endangle, linewidth=line_width, arrowheadfunction=negative_arrowhead(scale * relative_arrowhead_width))
    else
        # Unknown input
        arrow(centerpos, radius, startangle, endangle, arrowheadlength=0, linewidth=line_width)
    end
end


"""
   draw_edge(type::Int, x1::Float64, y1::Float64, x2::Float64, y2::Float64, offset::Bool, draw_config::Dict)

   Draws an edge between two distinct nodes. Note that this function is overloaded.

# Arguments (Required)
- `type::Int`: Type of interaction. 1 for positive interaction and -1 for negative interaction
- `x1::Float64`: x coordinate of first node
- `y1::Float64`: y coordinate of first node
- `x2::Float64`: x coordinate of second node
- `y2::Float64`: y coordinate of second node

# Arguments (Optional)
- `offset::Bool`: If true, the edge will be drawn slightly away from the line that joins the exact centers of nodes. If there are both forward and backward interactions between nodes, consider using this option to avoid edges being overlapped.
- `draw_config::Dict`: Specification of drawing styles
"""
function draw_edge(type::Int, x1::Float64, y1::Float64, x2::Float64, y2::Float64; offset::Bool = false, draw_config::Dict = draw_default_config)
    # Unpack drawing settings
    scale = draw_config[:scale]
    relative_node_radius = draw_config[:relative_node_radius]
    line_width = draw_config[:line_width]
    relative_path_offset = draw_config[:relative_path_offset]
    relative_space_between_arrows = draw_config[:relative_space_between_arrows]
    relative_arrowhead_width = draw_config[:relative_arrowhead_width]
    arrowhead_angle = draw_config[:arrowhead_angle]

    # Angle of the vector connecting given nodes, definition of angle follows svg's
    theta = atan(y2 - y1, x2 - x1)

    dx = (relative_node_radius + relative_path_offset) * cos(theta)
    dy = (relative_node_radius + relative_path_offset) * sin(theta)

    offx = relative_space_between_arrows / 2 * sin(theta) * offset
    offy = -1 * relative_space_between_arrows / 2 * cos(theta) * offset

    start_x = x1 + dx + offx
    start_y = y1 + dy + offy

    end_x = x2 - dx + offx
    end_y = y2 - dy + offy

    startpoint = Point(scale * start_x, scale * start_y)
    endpoint = Point(scale * end_x, scale * end_y)

    arrowheadlength = tan(π / 2 - arrowhead_angle) * relative_arrowhead_width * scale / 2

    if type == 1
        # Positive input
        arrow(startpoint, endpoint, arrowheadlength=arrowheadlength, arrowheadangle=arrowhead_angle, linewidth=line_width)
    elseif type == -1
        # Negative input
        arrow(startpoint, endpoint, arrowheadlength=0, linewidth=line_width)
        arrow(startpoint, endpoint, arrowheadlength=arrowheadlength, arrowheadangle=arrowhead_angle, linewidth=line_width, arrowheadfunction=negative_arrowhead(scale * relative_arrowhead_width))
    else
        # Unknown input
        arrow(startpoint, endpoint, arrowheadlength=0, linewidth=line_width)
    end
end


"""
   negative_arrowhead(width)

   Returns a function that can be passed to the `arrowheadfunction` keyword argument of the `Luxor.arrow` function to generate flat arrowheads.

# Arguments (Required)
- `width::Float64`: Width of the flat arrowhead

# Returns
- `Function`: A function to be passed to the `arrowheadfunction` keyword argument
"""
function negative_arrowhead(width::Float64)
    return ((originalendpoint, newendpoint, shaftangle) ->
        @layer begin
            translate(newendpoint)
            rotate(shaftangle + π / 2)
            ngon(O, width * 0.5, 2, 0, :stroke)
        end
    )
end


"""
   draw_cycle(type::Int, xs::AbstractArray, ys::AbstractArray, draw_config::Dict)

   Draws a Bezier curve around given points. This is to highlight cycles in interaction networks.

# Arguments (Required)
- `type::Int`: Type of cycle. 1 for positive feedback loops and -1 for negative feedback loops
- `xs:AbstractArray`: Array of x coordinates of nodes in a cycle
- `ys:AbstractArray`: Array of y coordinates of nodes in a cycle

# Arguments (Optional)
- `draw_config::Dict`: Specification of drawing styles
"""
function draw_cycle(type::Int, xs::AbstractArray, ys::AbstractArray; draw_config::Dict = draw_default_config)
    # Unpack drawing settings
    scale = draw_config[:scale]
    line_width = draw_config[:line_width]
    line_style = draw_config[:line_style]
    relative_arc_offset = draw_config[:relative_arc_offset]
    relative_arc_radius = draw_config[:relative_arc_radius]
    relative_space_between_arrows = draw_config[:relative_space_between_arrows]

    cycle_length = size(xs)[1]

    if type == 1
        line_color = draw_config[:line_color_positive_feedback]
    elseif type == -1
        line_color = draw_config[:line_color_negative_feedback]
    else
        line_color = draw_config[:line_color]
    end

    setcolor(line_color...)
    setline(line_width)
    setdash(line_style)

    if cycle_length == 1
        centerpos = Point(scale * xs[1] * (1 + relative_arc_offset), scale * ys[1] * (1 + relative_arc_offset))
        radius = scale * (relative_arc_radius + 0.05)

        arc(centerpos, radius, 0, 2π, action=:stroke)
    elseif cycle_length == 2
        theta = atan(ys[2] - ys[1], xs[2] - xs[1])

        offx = relative_space_between_arrows * sin(theta)
        offy = -1 * relative_space_between_arrows * cos(theta)

        arc_midpoint_1 = Point(scale * ((xs[1] + xs[2]) / 2 + offx), scale * ((ys[1] + ys[2]) / 2 + offy))
        arc_midpoint_2 = Point(scale * ((xs[1] + xs[2]) / 2 - offx), scale * ((ys[1] + ys[2]) / 2 - offy))

        p1 = Point(scale * xs[1], scale * ys[1])
        p2 = Point(scale * xs[2], scale * ys[2])

        path = makebezierpath([p1, arc_midpoint_1, p2, arc_midpoint_2])
        drawbezierpath(path, action=:stroke, close=true)
    elseif cycle_length > 2
        points = [Point(scale * xs[i], scale * ys[i]) for i in 1:cycle_length]

        path = makebezierpath(points)
        drawbezierpath(path, action=:stroke, close=true)
    end
end


"""
   draw_connectivity(connectivity::AbstractMatrix; coherent_nodes::AbstractArray = [], incoherent_nodes::AbstractArray = [], positive_cycles::AbstractArray = [], negative_cycles::AbstractArray = [], kwargs...)

   Draws a connectivity and optionally highlight node and cycle types.

# Arguments (Required)
- `connectivity::AbstractMatrix`: Connectivity matrix of a network

# Arguments (Optional)
- `coherent_nodes::AbstractArray`: Vector containing indices of nodes receiving coherent inputs
- `incoherent_nodes::AbstractArray`: Vector containing indices of nodes receiving incoherent inputs
- `positive_cycles::AbstractArray`: Vector containing positive feedback cycles. Each cycle should be given as a vector of node indices
- `negative_cycles::AbstractArray`: Vector containing negative feedback cycles
"""
function draw_connectivity(connectivity::AbstractMatrix; coherent_nodes::AbstractArray = [], incoherent_nodes::AbstractArray = [], positive_cycles::AbstractArray = [], negative_cycles::AbstractArray = [], kwargs...)
    n = size(connectivity)[1]

    # Input sanity checks
    nodes_with_coherency = vcat(coherent_nodes, incoherent_nodes)
    if !isempty(nodes_with_coherency) && findmax(nodes_with_coherency)[1] > n
        return error("Some nodes in the optional arguments are not represented in the connectivity matrix.")
    end

    nodes_in_any_cycles = vcat(positive_cycles..., negative_cycles...)
    if !isempty(nodes_in_any_cycles) && findmax(nodes_in_any_cycles)[1] > n
        return error("Some nodes in the optional arguments are not represented in the connectivity matrix.")
    end
    ambiguous_nodes = intersect(coherent_nodes, incoherent_nodes)
    if !isempty(ambiguous_nodes)
        return error("Coherency for following nodes is ambiguous: $(join(ambiguous_nodes, ',')).")
    end
    ambiguous_cycles = intersect(positive_cycles, negative_cycles)
    if !isempty(ambiguous_cycles)
        return error("Type for following cycles is ambiguous: $(join(ambiguous_cycles, ',')).")
    end

    draw_config = merge(draw_default_config, Dict(kwargs))

    output = draw_config[:output]
    ext = lowercase(splitext(output)[end])

    scaled_size = draw_config[:relative_canvas_size] * draw_config[:scale]

    xs, ys = get_polygon_vertex_coordinates(n)

    # Drawing routine
    function draw()
        # Draw nodes
        for i in 1:n
            if !(i in nodes_with_coherency)
                draw_node(0, xs[i], ys[i], draw_config=draw_config)
            end
        end
        for i in coherent_nodes
            draw_node(1, xs[i], ys[i], draw_config=draw_config)
        end
        for i in incoherent_nodes
            draw_node(-1, xs[i], ys[i], draw_config=draw_config)
        end

        # Draw self-loops
        for i in 1:n
            if connectivity[i, i] != 0
                draw_edge(connectivity[i, i], xs[i], ys[i], draw_config=draw_config)
            end
        end

        # Draw interactions
        for i in 1:n
            for j in i + 1:n
                forward = connectivity[j, i]
                backward = connectivity[i, j]

                forward != 0 && backward != 0 ? sep = true : sep = false
                # If both forward and backward paths exist, make the arrows to be separated from each other

                if forward != 0
                    draw_edge(forward, xs[i], ys[i], xs[j], ys[j], offset=sep, draw_config=draw_config)
                end

                if backward != 0
                    draw_edge(backward, xs[j], ys[j], xs[i], ys[i], offset=sep, draw_config=draw_config)
                end
            end
        end

        # Highlight cycles
        for c in positive_cycles
            draw_cycle(1, xs[c], ys[c], draw_config=draw_config)
        end
        for c in negative_cycles
            draw_cycle(-1, xs[c], ys[c], draw_config=draw_config)
        end
    end

    if ext == ".png"
        @png begin
            draw()
        end scaled_size scaled_size output
    elseif ext == ".svg"
        @svg begin
            draw()
        end scaled_size scaled_size output
    else
        return error("Support for the extension $ext is not implemented.")
    end
end

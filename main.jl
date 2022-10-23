#using GLMakie
#GLMakie.activate!()
using Plots

# function heatmap_contour_and_contourf()
#     x, y, z = peaks()
#     fig = Figure(resolution=(1200, 400))
#     axs = [Axis(fig[1, i]; aspect=DataAspect()) for i = 1:3]
#     hm = heatmap!(axs[1], x, y, z)
#     contour!(axs[2], x, y, z; levels=20)
#     contourf!(axs[3], x, y, z)
#     Colorbar(fig[1, 4], hm, height=Relative(0.5))
#     fig
# end
# JDS.heatmap_contour_and_contourf()

function insert(mem, xvec,yvec,zvec)
    for z in zvec
        for y in yvec
            for x in xvec
                mem[x,y,z] += 1
            end
        end
    end
    for x in xvec
        for y in yvec
            for z in zvec
                mem[z,y,x] += 16
            end
        end
    end
end 


N = 1000
P = 10
global errors = 0

CYCLES = 10000
println("preparing $CYCLES vectors")
w = @elapsed begin
    global CYCLES 
    mem = zeros(UInt8,N,N,N);
    cache = zeros(UInt8,N,N);
    # prepare test vectors
    x_input = zeros(Int64,CYCLES,P)
    y_input = zeros(Int64,CYCLES,P)
    z_input = zeros(Int64,CYCLES,P)

    for i in 1:CYCLES
        # ceil will make sure that the smallest number is 1 
        x_input[i,:] = sort(ceil.(Int,rand(P)*N))
        y_input[i,:] = sort(ceil.(Int,rand(P)*N))
        z_input[i,:] = sort(ceil.(Int,rand(P)*N))
    end 

    # check that the vectors dont have identical numbers 

    while true 
        redo = false
        for i in 1:CYCLES
            for k in 1:P-1
                if x_input[i,k] == x_input[i,k+1]
                    x_input[i,k+1] = ceil(Int,rand()*N)
                    x_input[i,:] = sort(x_input[i,:])
                    redo = true
                    break
                end
                if y_input[i,k] == y_input[i,k+1]
                    y_input[i,k+1] = ceil(Int,rand()*N)
                    y_input[i,:] = sort(y_input[i,:])
                    redo = true
                    break
                end
                if z_input[i,k] == z_input[i,k+1]
                    z_input[i,k+1] = ceil(Int,rand()*N)
                    z_input[i,:] = sort(z_input[i,:])
                    redo = true
                    break
                end
            end
        end
        redo || break
    end
end
speed = round(CYCLES/w, digits=2)
println("preparing $CYCLES vectors took $w seconds: $speed per sec")



w = @elapsed  for i in 1:CYCLES
    x = x_input[i,:]
    y = y_input[i,:]
    z = z_input[i,:]
    insert(mem,x,y,z)
end
speed = round(CYCLES/w, digits=2)
println("inserting $CYCLES vectors took $w seconds: $speed per sec")

function binarize(x,r,N,P)
    sorted = sort(x)
    plot(sorted)
    rankedmax = sorted[ N - P + 1]
    if rankedmax == 0
		rankedmax = 1
    end

	j=1
	for i = 1:N
        if x[i] >= rankedmax
            if (j<= P) 
                r[j] = i
            else
                println("Warning: j>P - rankedmax is $rankedmax,  retrieved vector is probably wrong")
            end
            j+=1
        end
    end
end

function binarize_by_certainty(x,r,N,P,c)
    sorted = sort(x)
    plot(sorted)
    max_val = maximum(sorted)
    certainty_thresh = max_val * c
    # rankedmax = sorted[ N - P + 1]
    # if rankedmax == 0
	# 	rankedmax = 1
    # end

	j=1
	for i = 1:N
        if x[i] >= certainty_thresh
            if (j<= P) 
                r[j] = i
            else
                println("Warning: j>P - certainty_thresh is $certainty_thresh,  retrieved vector is probably wrong")
            end
            j+=1
        end
    end
end

function query(mem,xvec,yvec,zvec)
    temp_vec_N = zeros(Int16,N)
    ret_vec = zeros(Int16,P)
    certainty = zeros(Int16,P)
    # notice that z sequencial access is better done using 
    # the x "like" locations - just looking at the high nibble
    # this allows us to speed up memory access since Julia is storing
    # sequece in column order see https://discourse.numenta.org/t/triadic-memory-a-fundamental-algorithm-for-cognitive-computing/9763/73
    # for more details..
    if ismissing(zvec)
        for x in xvec
            @simd  for y in yvec
                @views temp_vec_N += (mem[:,y,x] .>> 4)               
            end
        end
    elseif  ismissing(xvec)
        for y in yvec
            @simd  for z in zvec
                @views temp_vec_N += (mem[:,y,z] .& 0x0f)
            end
        end
    elseif ismissing(yvec)
        for x in xvec
            @simd for z in zvec
                @views temp_vec_N += (mem[x,:,z] .& 0x0f)                
            end
        end
    end
    
    binarize(temp_vec_N,ret_vec, N,P)
    return ret_vec
end 

# this query1 function operates sligtly different from the original query 
# function. Insteap of binarize for the first P values in the sorted summation vectors
# it takes the first few values with > certainty% which is 50% by default
function query1(mem,xvec,yvec,zvec)
    temp_vec_N = zeros(Int16,N)
    ret_vec = zeros(Int16,P)
    certainty = zeros(Int16,P)
    # notice that z sequencial access is better done using 
    # the x "like" locations - just looking at the high nibble
    # this allows us to speed up memory access since Julia is storing
    # sequece in column order see https://discourse.numenta.org/t/triadic-memory-a-fundamental-algorithm-for-cognitive-computing/9763/73
    # for more details..
    if ismissing(zvec)
        for x in xvec
            @simd  for y in yvec
                @views temp_vec_N += (mem[:,y,x] .>> 4)               
            end
        end
    elseif  ismissing(xvec)
        for y in yvec
            @simd  for z in zvec
                @views temp_vec_N += (mem[:,y,z] .& 0x0f)
            end
        end
    elseif ismissing(yvec)
        for x in xvec
            @simd for z in zvec
                @views temp_vec_N += (mem[x,:,z] .& 0x0f)                
            end
        end
    end
    
    # use a default 50%
    binarize_by_certainty(temp_vec_N,ret_vec, N,P,0.5)
    return filter(!iszero, ret_vec)
end 

# Query Z given x,y
####################
errors = 0
@time w = @elapsed for i in 1:CYCLES
    x = x_input[i,:]
    y = y_input[i,:]
    z = z_input[i,:]
    qz = query(mem,x,y,missing)
     if  qz != z
        global errors += 1
    end
end
speed = round(CYCLES/w, digits=2)
println("querying z $CYCLES vectors took $w seconds: $speed per sec, $errors errors")


# Query Y  given x,z
####################
errors = 0
@time w = @elapsed for i in 1:CYCLES
    x = x_input[i,:]
    y = y_input[i,:]
    z = z_input[i,:]
    qy = query(mem,x,missing,z)
     if  qy != y
        global errors += 1
    end
end
speed = round(CYCLES/w, digits=2)
println("querying y $CYCLES vectors took $w seconds: $speed per sec, $errors errors")

# Query X given z,y
####################
errors = 0
@time w = @elapsed for i in 1:CYCLES
    x = x_input[i,:]
    y = y_input[i,:]
    z = z_input[i,:]
    qx = query(mem,missing,y,z)
     if  qx != x
         global errors += 1
    end
end
speed = round(CYCLES/w, digits=2)
println("querying x $CYCLES vectors took $w seconds: $speed per sec, $errors errors")

# another example of flexible sdr size insertion and queries
insert(mem,[1,2,3],[4,5,6],[7,8,9])
query1(mem,[1,2,3,5,9],missing,[1,7,8,9])
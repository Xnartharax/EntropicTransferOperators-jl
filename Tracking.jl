function greedyscan(frames; tol=0.01, minlength=20)
    ls = []
    for (k, (frame1, frame2)) in enumerate(zip(frames[1:end-1], frames[2:end]))
        candidates = Set(1:length(frame2))
        done = Set()
        ends = lineends(ls)
        for (l, e) in ends
            x = frame1[e[2]]
            j = pickcandidate!(candidates, j -> abs(frame2[j] - x), tol=tol)
            push!(l, (k+1, j))
            push!(done, e)
        end
        for (i, x) in enumerate(frame1) 
            if i  ∉ done
                j = pickcandidate!(candidates, j -> abs(frame2[j] - x), tol=tol)
                push!(ls, [(k, i), (k+1, j)])
            end
        end
    end
    chopped = chop.(ls)
    pruned = pruneminlength(chopped, minlength)
    detangle.(pruned)
end

function greedyscan2(frames; tol=0.01, minlength=20)
    ls = []
    for (k, (frame1, frame2)) in enumerate(zip(frames[1:end-1], frames[2:end]))
        candidates = Set(1:length(frame2))
        done = Set()
        ends = lineends(ls)
        for (l, e) in ends
            x1 = e[3]
            x0 = l[end-1][3]
            j = pickcandidate!(candidates, j -> abs(2*x1 -x0 - frame2[j]), tol=tol)
            if j != -1
                push!(l, (k+1, j, frame2[j]))
            else
                push!(l, (k+1, j, -1))
            end
            push!(done, e)
        end
        for (i, x) in enumerate(frame1) 
            if i  ∉ done
                j = pickcandidate!(candidates, j -> abs(frame2[j] - x), tol=tol)
                if j != -1
                    push!(ls, [(k, i, x), (k+1, j, frame2[j])])
                end
            end
        end
    end
    chopped = chop.(ls)
    pruned = pruneminlength(chopped, minlength)
    detangle.(pruned)
end

function greedyscan3(frames, vecs; tol=0.01, minlength=20, α =1)
    ls = []
    for (k, (frame1, frame2)) in enumerate(zip(frames[1:end-1], frames[2:end]))
        candidates = Set(1:length(frame2))
        done = Set()
        ends = lineends(ls)
        for (l, e) in ends
            x1 = e[3]
            x0 = l[end-1][3]
            j = pickcandidate!(candidates, j -> abs(2*x1 -x0 - frame2[j]) + α*norm(vecs[k+1][e[2]]- vecs[k+1][j]), tol=tol)
            if j != -1
                push!(l, (k+1, j, frame2[j]))
            else
                push!(l, (k+1, j, -1))
            end
            push!(done, e[2])
        end
        for (i, x) in enumerate(frame1) 
            if i  ∉ done
                j = pickcandidate!(candidates, j -> abs(frame2[j] - x), tol=tol)
                if j != -1
                    push!(ls, [(k, i, x), (k+1, j, frame2[j])])
                end
            end
        end
    end
    chopped = chop.(ls)
    pruned = pruneminlength(chopped, minlength)
    detangle.(pruned)
end

lineends(ls) = [(l, l[end]) for l in ls if l[end][2] != -1]

chop(l) = [e for e in l if e[2] !=-1]

detangle(l) = [e[1] for e in l], [e[2] for e in l]

pruneminlength(ls, minlength) = [l for l in ls if length(l) > minlength]
pruneminlength2(ls, minlength) = [l for l in ls if length(l[1]) > minlength]
function pickcandidate!(candidates, dist; tol=0.01)
    if length(candidates) == 0
        return -1
    end 
    d, i = findmin(dist, collect(candidates))
    picked = collect(candidates)[i]
    if d < tol
        delete!(candidates, picked)
        return picked
    end
    return -1
end

function longesthorizontal(l, λs; tol=0.01)
    xs=0.8:tol/5:1
    horizontals = [horizontal(l, λs, x; tol=tol) for x in xs]
    _, i = findmax(length, horizontals)
    return detangle(horizontals[i])
end

function horizontal(l, λs, x; tol=0.01)
    hor = []
    for (i,j) in zip(l[1], l[2])
        if abs(λs[i][j] -x) <= tol
            push!(hor, (i, j))
        end
    end
    return hor
end

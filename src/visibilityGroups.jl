module visibilityGroups

using attitudeFunctions
using LinearAlgebra

export findVisGroup, findAllVisGroups

function findVisGroup(sat,scen,attitude)

    n = scen.obsNo
    m = sat.facetNo

    (sunVec,obsVecs) = toBodyFrame(attitude,scen.sunVec,scen.obsVecs)

    visGroup = Array{Bool,2}(undef,n+1,m)

    for i = 1:n+1
        for j = 1:m
            if size(obsVecs) == (3,n)
                if i == n+1
                    if dot(sunVec,view(sat.nvecs,:,j)) > 0
                        visGroup[i,j] = true
                    else
                        visGroup[i,j] = false
                    end
                elseif dot(view(obsVecs,:,i),view(sat.nvecs,:,j)) > 0
                    visGroup[i,j] = true
                else
                    visGroup[i,j] = false
                end
            elseif size(obsVecs) == (n,)
                if i == n+1
                    if dot(sunVec,sat.nvecs[j]) > 0
                        visGroup[i,j] = true
                    else
                        visGroup[i,j] = false
                    end
                elseif dot(obsVecs[i],sat.nvecs[j]) > 0
                    visGroup[i,j] = true
                else
                    visGroup[i,j] = false
                end
            end
        end
    end
    return visGroup
end

function findAllVisGroups(sat,scen)

    n = scen.obsNo
    m = sat.facetNo

    vgnmax = 2^(n+m)





end



end # module

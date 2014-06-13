module KDTree

#A KDTree is a recursive data structure/geometric algorithm that
#impliements nearest-neighbor search in N log N time.

#Port of COS 226 KDTree implementation following the structure in tree.jl.
#https://github.com/JuliaLang/DataStructures.jl/blob/master/src/tree.jl

#Keys are array representations of (x,y) coordinates.

using Rect
import Base: haskey, getindex,setindex!, delete!



abstract Tree{K,V}

type EmptyKDTree{K,V} <: KDTree{K,V}
end

type TreeNode{K,V} <: KDTree{K,V}
        key::K
        val::V
        rect
        lb::KDTree{K,V}
        rt::KDTree{K,V}
end

type KDTree{K,V}
       root::Tree{K,V}

       KDTree() = new(EmptyKDTree{K,V}())
end

haskey(t::EmptyKDTree{K,V}()) = false
haskey(t::KDTree,key) = haskey(t.root,key,true)

#t.key is the current node we are examining (Node x in the original Java implementation)
#key is the key we are ultimately looking for (key in original Java implementation)
function haskey(t::Treenode, key, orientation::Bool)
        if isequal(t.key,key) #objectwise comparison since keys are likely to be 2-element arrays
                return true
        end

        #Orientation alternates at each level and governs whether we go to left/right tree or
        #top/bottom tree
        if orientation
                cmp = t.key[1] - key[1];
        else
                cmp = t.key[2] - key[2];
        end

        if cmp < 0
                return haskey(t.key.lb, key, !orientation);
        else
                return haskey(t.key.rt, key, !orientation);
        end
end

getindex(t::EmptyKDTree, k) = throw(KeyError(k))
getindex(t::KDTree,k,true) = t.root[k] #Start searching from the root with "true" orientation

function getindex(t::TreeNode, key,orientation::Bool)
        if isequal(t.key,key)
                return t.data
        end
        # Compare the keys
        if orientation
                cmp = t.key[1] - key[1];
        else
                cmp = t.key[2] - key[2];
        end
        # Return left-bottom key or right-top depending on key compare
        if cmp < 0
                return t.key.lb[key];
        else
                return t.key.rt[key];
        end
end

setindex!{K,V}(t::EmptyKDTree{K,V},v,k) = TreeNode{K,V}(k, v, t, t)
#Helper function to start the recursive calls
setindex!(t::KDTree,v,k, orientation) = (t.root = setindex!(t.root,v,k, true); t)

#Set index k to value v
function setindex!(t::TreeNode, v, k, orientation::Bool)
        #Update value if key is already in KDTree
        if isequals(t.key,k)
                t.data = v
        end

        #Compare keys
        if orientation
                cmp = t.key[1] - k[1];
        else
                cmp = t.key[2] - k[2];
        end

        if cmp < 0
                t.lb = setindex!(t.lb, v, k, !orientation)
        else
                t.rt = setindex!(t.rt, v, k, !orientation)
        end
        return t #I'm assuming it returns t
end

#Delete key k and its value v from the KDTree (via Hibbard delete?)
delete!(t::EmptyKDTree, k) = throw(KeyError(k))
#Seed the first call with the true orientation
delete!(t::KDTree, k) = (t.root = delete!(t.root, k, true); t)


function delete!(t::TreeNode, k, orientation::Bool)

        if orientation
                cmp = t.key[1] - k[1];
        else
                cmp = t.key[2] - k[2];
        end

        if isequal(t.key, k)
                if isa(t.rt, EmptyKDTree)
                        t = t.lb;
                elseif isa(t.lt, EmptyKDTree)
                        t = t.rt;
                else
                        r = t.rt;
                        t = t.lb
                        treeinsert!(t, r)
                end
        elseif cmp < 0
                t.lb = delete!(t.lb, k, !orientation);
        else
                t.rt = delete!(t.rt, k, !orientation)
        end
        t
end


treeinsert!(t::EmptyKDTree, r::TreeNode) = r;

#Keeping the naming conventions consistent with DataStructures.jl's Tree
#implementation. t corresponds to the current node being examined;
#r is the node being inserted.
function treeinsert!(t::TreeNode, r::TreeNode, orientation::Bool)

        if orientation
                cmp = r.key[1] - t.key[1];
        else
                cmp = r.key[2] - t.key[2];
        end

        if cmp < 0 #if r.key < t.key
                t.lb = treeinsert!(t.lb,r,!orientation);
        else
                t.rt = treeinsert!(t.rt,r,!orientation);
        end
        t

end


end

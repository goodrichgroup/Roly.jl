module ReverseSearch

import SciMLBase
using Base.Threads: nthreads, @spawn

export RSSystem, RSState, RSIterator, reversesearch

include("./reverse_search.jl")
end
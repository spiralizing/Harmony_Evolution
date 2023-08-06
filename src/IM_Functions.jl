#Information-based measures
#
"""
    get_tprobs_prev(f_eij, prev_eij, beta)

    Returns a dictionary of probabilites (normalized frequency) from a dictionary
    f_eij and the dictionary of counts prev_eij (previous pool), using a Laplacian
    smoothing within the transition space n_states or PT_exj.
"""
function get_tprobs_prev(f_eij, prev_eij, beta)
    f_ei = Dict{Any, Int64}()
    for eij in keys(prev_eij)
        ei = split(eij,"+")[1]
        f_ei[ei] = get(f_ei,ei,0) + prev_eij[eij]
    end
    p_ij = Dict{Any,Float64}()
 
    for eij in keys(f_eij)
        e1 = split(eij,"+")[1]
        if haskey(prev_eij, eij)
            prob_ij = (prev_eij[eij] + beta) / (f_ei[e1] + beta*24) #case when the transition exists
        elseif haskey(f_ei, e1) #transition doesnt exists but first token does
            prob_ij = beta / (f_ei[e1] + beta*24)
        else
            prob_ij = beta / (beta*24) #completely new pair of elements
        end
        p_ij[eij] = get(p_ij,eij,prob_ij)
    end

    return sort(p_ij)
end
###
"""
    get_tprobs(f_eij)

    Returns probabilities (normalized frequencies) from a dictionary of transition
    counts, the transitions are normalized by the transitions of the first state,
    constructing an stochastic matrix without a matrix.
"""
function get_tprobs(f_eij)
    f_ei = Dict{Any,Int64}()
    p_ei = Dict{Any,Float64}()
    n_ei = []
    for eij in keys(f_eij)
        ei = split(eij,"+")[1]
        push!(n_ei,ei)
        f_ei[ei] = get(f_ei,ei,0) + f_eij[eij]
        p_ei[ei] = get(f_ei,ei,0) + f_eij[eij]
    end
    p_e = Dict{Any,Float64}()
    n_e = sum(values(f_ei))
    for e in keys(f_ei)
        prob_e = f_ei[e]/n_e
        p_e[e] = get(p_e, e, prob_e)
    end
   # m_p = Float64[]
   # for i in 2:length(n_ei)
   #     push!(m_p, p_e[n_ei[i]])
   # end
    m_p = Float64[]
    p_ij = Dict{Any,Float64}()
    for eij in keys(f_eij)
        e1 = split(eij,"+")[1]
        prob_ij = f_eij[eij] / f_ei[e1]
        p_ij[eij] = get(p_ij, eij, prob_ij)
    end
    p_ij = sort(p_ij)
    for eij in keys(p_ij)
        e1 = split(eij,"+")[1]
        push!(m_p, p_e[e1])
    end
    return m_p, sort(p_ij)
end
###--
"""
    KLD_Piece_Pool(piece, prev_pool, PT_exj)

    Returns the Kullback-Leibler divergence between a piece (counts of transitions)
    and a pool of previous pieces (counts of transitions).
"""
function KLD_Piece_Pool(piece, prev_pool, PT_exj)
    #P_piece = sort(get_tprobs_prev(piece,piece,PT_exj))
    P_piece = get_tprobs(piece)[2]
    Q_prev = get_tprobs_prev(piece, prev_pool, PT_exj)
    #KL divergence, K(p|q) = Σ p(x)log(p(x)/q(x)P
    P = collect(values(P_piece)) #distribution of the transitions in the piece
    Q = collect(values(Q_prev)) #probabilites of the transitions for the piece in pool
    #T = collect(values(piece))./n_e #times the transition is made.
    KLD = mapreduce((x,y) -> x* log10(x/y), +, P, Q)
    return KLD 
end
##--
"""
    KLD_rate(piece,prev_pool,beta)

    Computes the Kullback-Leibler Divergence rate between two stochastic matrices:
    KL(P||Q) = Σ μ p(x)log(p(x)/q(x)), where μ is the asymptotic distribution of marginals.
    Beta is a parameter for additive smoothing of the distribution Q.
"""
function KLD_rate(piece, prev_pool, beta)
    μ_i, P_ij = get_tprobs(piece) #asymptotic distribution of marginals, transition probabilities
    Q_prev = get_tprobs_prev(piece, prev_pool, beta)
    P = collect(values(P_ij))[2:end] #distribution of the transitions in the piece
    Q = collect(values(Q_prev))[2:end] #probabilites of the transitions for the piece in pool
    #KL divergence rate, K(p|q) = Σ μ p(x)log(p(x)/q(x)P
    return mapreduce((x,y,z) -> z * x* log10(x/y), +, P, Q, μ_i)
end
##--
function KLD_Piece_Piece(piece, PT_exj)
    Q_prev = sort(get_tprobs_prev(piece,piece,PT_exj))
    P_piece = sort(get_tprobs(piece))
    n_e = sum(collect(values(piece)))
    #KL divergence, K(p|q) = Σ p(x)log(p(x)/q(x)
    P = collect(values(P_piece)) #distribution of the transitions in the piece
    Q = collect(values(Q_prev)) #probabilites of the transitions for the piece in pool
    T = collect(values(piece)) #times the transition is made.
    KLD = mapreduce((x,y,z) -> z*(x* log10(x/y)), +, P, Q, T)
    return KLD
end


###--
"""
    IC_Piece_Pool(piece, prev_pool, PT_exj)

    Returns the information content of a piece (dictionary of pair transitions)
    with the probabilities computed by the MLE in the pool of previous pieces,
    with a laplacian smoothing over the whole transition space.
"""
function IC_Piece_Pool(piece, prev_pool, PT_exj; normed=true)
    Q_prev = sort(get_tprobs_prev(piece, prev_pool, PT_exj))
    T = collect(values(piece)) #times the transition is made.
    Q = collect(values(Q_prev))
    #information content -Σp(x)
    IC = mapreduce((x,y)-> - y * log10(x), +,Q,T)
    if normed
        return IC / sum(T)
    else
        return IC
    end
end

###--
"""
    get_tspace(Xi)

    Returns the number of different transitions x->j in the probability space
    for each x in the alphabet (different states).
"""
function get_tspace(Xi)
    PT_exj = Dict{Any,Int64}() #number of possible transitions in the whole space.
    a_ei = map(x-> split(x,"+")[1],collect(keys(Xi)))
    for i = 1:length(a_ei)
        PT_exj[a_ei[i]] = get(PT_exj, a_ei[i],0)+1
    end
    return PT_exj
end

###---
"""
    tran_freq_to_mat(f_eij)

    Returns two outputs, first is a transition matrix for the pairs of transitions i -> j from a 
    dictionary f_eij, and the second output is the list of the tokens (or elements) of the sequence.
"""
function tran_freq_to_mat(f_eij)
    f_ei = Any[]
    f_ej = Any[]
    eij_c = Int64[]
    for eij in keys(f_eij)
        ei,ej = split(eij,"+")
        push!(f_ei, ei); push!(f_ej, ej)
        push!(eij_c,f_eij[eij])
    end

    ##--
    tok = sort(unique(f_ej))
    t_mat = zeros(Int64,length(tok),length(tok))
    ##--
    for i = 1:length(tok)
        ti = findall(x-> x==tok[i], f_ei)
        if isempty(ti); continue; end
        for j = 1:length(ti)
            tj = findfirst(x-> x==f_ej[ti[j]], tok)
            t_mat[i,tj] = eij_c[ti[j]]
        end
    end
    return t_mat, tok
end
###############################CODEWORDS
#NOVELTY FUNCTIONS


function get_novelty_measure(piece, Om, N_cws)
    ν_ξ = 0 #the value for the novelty of piece ξ
    ξ = piece
    Ω = collect(keys(Om)) #conventional pool
    Α = collect(setdiff(keys(ξ), keys(Om))) #Novelty pool, elements in the piece that are not in previous works.
    γ = permutedims(hcat(map(x->split(x, "+"), Ω)...)) #re-ordering the pools into 2-dim array.
    α = permutedims(hcat(map(x->split(x, "+"), Α)...))

    s_m = sum(values(ξ)) #size of the piece
    for k in keys(ξ)
        γ_1, γ_2 = split(k,"+") #separate chord 1 and chord 2
        γ_Ω = findall(x->x==γ_1, γ[:,1]) #pairs that start with γ_1 in the Codewordspace
        if isempty(Α)
            s_γ = 0
        else
            γ_a = findall(x->x==γ_1, α[:,1]) #pairs that start with γ_1 in the novelty pool
            s_γ = length(γ_a) #number of pairs that start with γ_1 in novelty pool.
        end
        for e in γ_Ω  #looping over elements that start with γ_1 in the conventional pool
            s_γ += Om[Ω[e]]  #number of pairs in the conventional pool
        end
        s_γ += N_cws
        ix = findall(x->x==k,Ω) #finding the element k in the conventional pool
        if isempty(ix)
            ν_ξ += ξ[k] * log10(s_γ)  #this is the case when the element is a novelty
        else
            ν_ξ += ξ[k] * log10(s_γ / (Om[Ω[ix[1]]] + 1) ) #when the element is in the conventional pool
        end
    end
    return ν_ξ / s_m #dividing by the size of the piece
end

function get_novelty_measure(piece, Om, N_cws; nov_w=1)
    ν_ξ = 0 #the value for the novelty of piece ξ
    ξ = piece
    Ω = collect(keys(Om)) #conventional pool
    Α = collect(setdiff(keys(ξ), keys(Om))) #Novelty pool, elements in the piece that are not in previous works.
    γ = permutedims(hcat(map(x->split(x, "+"), Ω)...)) #re-ordering the pools into 2-dim array.
    α = permutedims(hcat(map(x->split(x, "+"), Α)...))

    s_m = sum(values(ξ)) #size of the piece
    for k in keys(ξ)
        γ_1, γ_2 = split(k,"+") #separate chord 1 and chord 2
        γ_Ω = findall(x->x==γ_1, γ[:,1]) #pairs that start with γ_1 in the Codewordspace
        if isempty(Α)
            s_γ = 0
        else
            γ_a = findall(x->x==γ_1, α[:,1]) #pairs that start with γ_1 in the novelty pool
            s_γ = length(γ_a) #number of pairs that start with γ_1 in novelty pool.
        end
        for e in γ_Ω  #looping over elements that start with γ_1 in the conventional pool
            s_γ += Om[Ω[e]] #number of pairs in the conventional pool
        end
        s_γ += N_cws * nov_w #number of different keys the current key can go to.
        ix = findall(x->x==k,Ω) #finding the element k in the conventional pool
        if isempty(ix)
            ν_ξ += ξ[k]  * log10(s_γ / nov_w)  #this is the case when the element is a novelty
        else
            ν_ξ += ξ[k] * log10((s_γ+nov_w) / (Om[Ω[ix[1]]] + nov_w) ) #when the element is in the conventional pool
        end
    end
    return ν_ξ / s_m #dividing by the size of the piece
end

function indice(tmax::Int64) #esta funcion solo regresa el arreglo que lleva el indice(numero) de las notas (o el eje x)
    ind = Array{Float64}(undef,tmax)
    for i = 1:tmax
        ind[i] = convert(Float64, i)
    end
    return(ind)
end

function filter_undef!(voces)
    tt = map(x-> isassigned(x), voces)
    deleteat!(voces,findall(x->x==false,tt))
end

#Se definen las funciones
function gaps_notas!(s, q) #La funcion toma de entrada el arreglo de notas y corrige los pequenios gaps que hay entre notas que terminan y que empiezan
    for v in s
        for i = 1:(size(v)[1]-1)
            m = abs(v[i,2] - v[i+1,1])
            #printlrintln(Voz[i,2],'\t',Voz[i+1,1],'\t',m)
            if m < q
                v[i,2] = v[i+1,1]
                #println(Voz[i,2],'\t', Voz[i+1,1])
            end
        end
      end
end
#########################################################################################################################################################################
function gaps_silencios!(s,q) # La funcion toma de entrada el arreglo de notas y corrige los pequenios gaps que puede haber entre los silencios y donde empiezan las notas
    sm = q
    for v in s
        for i = 1:(size(v)[1]-1) #corrige los gaps que hay entre silencios y notas
            m = mod(v[i,2] - v[i+1,1],sm)
            if m > 0
                v[i,2] = v[i,2] + (sm - m)
            end
        end
        x = length(v[:,1])
        n = mod(v[x,2] - v[x,1], sm) #aqui se corrige la ultima nota del arrelgo
        if n > 0
            v[x,2] = v[x,2] + (sm - n)
        end
    end
end
#################################################################################################################################################
function min_voces(n_vs, dv) #la funcion regresa el minimo valor de duracion
    ns = size(n_vs)[1]
    mins = Array{Float64}(undef,ns)
    for i = 1:ns
        dif = n_vs[i][:,2] - n_vs[i][:,1]
        difn = filter(x->x >= dv,dif)
        mins[i] = minimum(difn)
    end
    return minimum(mins)
end
################################################################################################################################################
function max_tempo(Voces, ns) #la funcion encuentra el tiempo maximo de termino de notas, es decir donde termina la pieza
    Tfinal = Array{Float64}(undef,ns)
    for i = 1:ns
        Tfinal[i] = maximum(Voces[i][:,2])
    end
    return maximum(Tfinal)
end
################################################################################################################################
#La siguiente funcion corrige los numeros fraccionarios que existen en las Voces
function rounding!(voces, q)
    nv = size(voces)[1]
    for i =1:nv
        voces[i][:,1] = map(x-> ceil(x/q), voces[i][:,1])
        voces[i][:,2] = map(x-> ceil(x/q), voces[i][:,2])
    end
end
#######################################################################################################################################
function serie_notas(Voz, tmax) #construye la serie de tiempo teniendo como entrada el arreglo de inicio - final - voz
    serie = zeros(tmax)
    ini = Voz[:,1]
    dura = Voz[:,2] - Voz[:,1] #este es el arreglo de duraciones
    for i = 1:(length(dura)) #aqui se asignan las notas
        j = 0
        x = convert(Int64, ini[i]) + 1
        while dura[i] > 0
            serie[x+j] = Voz[i,3]
            j += 1
            dura[i] -= 1
        end
    end
    return(serie)
end

###################################################################################
#New functions to construct time series and pitch sequences
function get_onoff_notes(s)
    ini = s[s[:,3].==" Note_on_c", :] #getting the the starting places
    fin = s[s[:,3].==" Note_off_c",:] #ending
    dn = unique(ini[:,5]) #different notes.
    ndn = length(dn) #number of different notes
    posin_dn = Array{Vector}(undef, ndn) #initialize arrays for getting information of each note
    posfn_dn = Array{Vector}(undef, ndn)
    in_fi = Array{Matrix}(undef,ndn)
    for i = 1:ndn       #for each different note it gets when are they played
        tmp = ini[ini[:,5].==dn[i],:][:,2]
        tmp2 = fin[fin[:,5].==dn[i],:][:,2]
        l1 = length(tmp); l2 = length(tmp2)
        if l1 != l2
            nt = min(l1,l2)
            note = [dn[i] for j=1:nt]
            in_fi[i] = [tmp[1:nt] tmp2[1:nt] note]
        else
            note = [dn[i] for j = 1:length(tmp)]
            in_fi[i] = [tmp tmp2 note]
        end
    end
    nmat = vcat(in_fi...) #getting all notes in the same array
    return nmat[sortperm(nmat[:,1]),:]  #returns the array sorted by appereance.
end
################################################################################
function get_new_voice(v)
    new_v = []
    ids = []
    iids = [i for i=1:size(v)[1]]
    c = 1
    while c <= size(v)[1]-1
        if v[c+1,1] <= v[c,1]
            #println(v[c,:],'\t',v[c+1,:])
            push!(new_v, v[c+1,:])
            push!(ids,c+1)
            c+=1
        else
            #push!(ids,c)
            c+=1
        end
        #c += 1
        #println(c)
    end
    deleteat!(iids, ids)
    old_v = v[iids,:]
    new_v = convert(Array{Any,2},transpose(reshape(vcat(new_v...),3,:)))
    return [old_v, new_v]
end
function find_more_voices(v)
    v_1 =[]
    nvs = get_new_voice(v)
    push!(v_1,nvs[1])
    lv = nvs[2]
    br = false
    while br == false #this is to check if there are more voices
        ns = get_new_voice(lv)
        push!(v_1, ns[1])
        br = isempty(ns[2])
        if br; break; end
        #println(br)
        lv = ns[2]
    end
    return v_1
end
################################################################################
function get_onon_notes(s)
    ini = s[s[:,6].!=0, :] #getting the the starting places
    fin = s[s[:,6].==0,:] #ending
    dn = unique(ini[:,5]) #different notes.
    ndn = length(dn) #number of different notes
    posin_dn = Array{Vector}(undef, ndn) #initialize arrays for getting information of each note
    posfn_dn = Array{Vector}(undef, ndn)
    in_fi = Array{Matrix}(undef,ndn)
    for i = 1:ndn       #for each different note it gets when are they played
        tmp = ini[ini[:,5].==dn[i],:][:,2]
        tmp2 = unique(fin[fin[:,5].==dn[i],:][:,2])#added the unique to be sure they only finish "once", weird stuff was happening because of this.
        l1 = length(tmp); l2 = length(tmp2)
        if l1 != l2
            nt = min(l1,l2)
            note = [dn[i] for j=1:nt]
            in_fi[i] = [tmp[1:nt] tmp2[1:nt] note]
        else
            note = [dn[i] for j = 1:length(tmp)]
            in_fi[i] = [tmp tmp2 note]
        end
    end
    nmat = vcat(in_fi...) #getting all notes in the same array
    return nmat[sortperm(nmat[:,1]),:]  #returns the array sorted by appereance.
end
################################################################################


#next function gets the multivariate time series of a pieces from a csv file
function get_time_series(s)
    sd = 16
    #mq / sd
    #sd is the smallest reference duration, it would divide the quarter of a note.
    mq = s[1,6] #quarter of a note
    nv = s[findlast(x-> x==" Note_on_c",s[:,3]),1] #- 1 #estimates how many voices are in the midi
    voces = Array{Matrix}(undef,nv) #initialze an array
    #next function is to get the time in miliseconds when the pitch starts and ends.
    if nv ==1
        #nv = 1
        if isempty(findall(x-> x==" Note_off_c",s[:,3]))
            b = s[s[:,1].==1,:] #takes the events of the voice i
            mat = b[b[:,3].==" Note_on_c", :]
            mat = mat[findall(x->x!="", mat[:,5]),:]
            voces = get_onon_notes(mat) #construct an array of of information of initial time, finish time, pitch and intensity.
        else
            mat = s[s[:,1].==1,:]
            mat = mat[findall(x->x!="", mat[:,5]),:]
            voces = get_onoff_notes(mat)
        end
    else
        voces = Array{Matrix}(undef,nv) #initialze an array
        if isempty(findall(x-> x==" Note_off_c",s[:,3]))
            #checks if the midi has events of note_off
            for i = 1:(nv) #if does not, it construct the series in this way
                b = s[s[:,1].==i,:] #takes the events of the voice i
                if size(b[b[:,3].==" Note_on_c", :])[1] == 0; continue; end #checks if there are notes in the channel
                mat = b[b[:,3].==" Note_on_c", :]
                mat = mat[findall(x->x!="", mat[:,5]),:]
                voces[i] =  get_onon_notes(mat)#construct an array of information of initial time, finish time, pitch and intensity.
            end
        else
            for i = 1:nv #if has note_off events , it does this way
                ini = findfirst(x-> x==" Note_on_c", s[s[:,1].==i,3]) #initial time
                fin = findlast(x-> x==" Note_off_c", s[s[:,1].==i,3]) #finish time
                if ini == 0 || fin == 0; continue; end
                mat = s[s[:,1].==i,:]
                mat = mat[findall(x->x!="", mat[:,5]),:]
                voces[i] = get_onoff_notes(mat) #construct an array of information of initial time, finish time, pitch and intensity
            end
        end
    end
    if nv > 1
        filter_undef!(voces)
        filter!(x->length(x)>0,voces)
        n_vs = []
        nv = length(voces)
        for i=1:nv
            m_v = find_more_voices(voces[i])
            for j =1:length(m_v)
                push!(n_vs, m_v[j])
            end
        end
    else
        try
            n_vs = []
            m_v = find_more_voices(voces[1])
            for j =1:length(m_v)
                push!(n_vs, m_v[j])
            end
        catch
            n_vs = []
            m_v = find_more_voices(voces)
            for j =1:length(m_v)
                push!(n_vs, m_v[j])
            end
        end
    end
    filter!(x->length(x)>0,n_vs)
    nv = size(n_vs)[1] #the real number of voices
    q = round(Int, mq / (mq / min_voces(n_vs, mq / sd)))
    #q = mq/ round(Int, mq / min_voces(n_vs, mq / sd)) #defines the unit of time  q
    gaps_notas!(n_vs,q) #fixes some gaps of miliseconds between the notes
    gaps_silencios!(n_vs,q) #fixes some gaps of milisecons of rests
    rounding!(n_vs, q) #rounds up some possible decimals in the time
    tmax = round(Int,max_tempo(n_vs,nv)) #gets the total lenght of the piece
    series = zeros(tmax,nv+1)  #this is the initialization of the output array
    series[:,1] = indice(tmax) #this is just the index
    for i=2:(nv+1) #constructs the output time series
        series[:,i] = serie_notas(n_vs[i-1],tmax)
    end
    return series
end


################################################################################
function chunk_voices(s, n_chunks)
    mq = s[1,6] #quarter of a note
    nv = s[findlast(x-> x==" Note_on_c",s[:,3]),1] #- 1 #estimates how many voices are in the midi
    voces = Array{Matrix}(undef,nv) #initialze an array
    #next function is to get the time in miliseconds when the pitch starts and ends.
    if nv ==1
        #nv = 1
        if isempty(findall(x-> x==" Note_off_c",s[:,3]))
            b = s[s[:,1].==1,:] #takes the events of the voice i
            mat = b[b[:,3].==" Note_on_c", :]
            mat = mat[findall(x->x!="", mat[:,5]),:]
            voces = get_onon_notes(mat) #construct an array of of information of initial time, finish time, pitch and intensity.
        else
            mat = s[s[:,1].==1,:]
            mat = mat[findall(x->x!="", mat[:,5]),:]
            voces = get_onoff_notes(mat)
        end
    else
        voces = Array{Matrix}(undef,nv) #initialze an array
        if isempty(findall(x-> x==" Note_off_c",s[:,3]))
            #checks if the midi has events of note_off
            for i = 1:(nv) #if does not, it construct the series in this way
                b = s[s[:,1].==i,:] #takes the events of the voice i
                if size(b[b[:,3].==" Note_on_c", :])[1] == 0; continue; end #checks if there are notes in the channel
                mat = b[b[:,3].==" Note_on_c", :]
                mat = mat[findall(x->x!="", mat[:,5]),:]
                voces[i] =  get_onon_notes(mat)#construct an array of information of initial time, finish time, pitch and intensity.
            end
        else
            for i = 1:nv #if has note_off events , it does this way
                ini = findfirst(x-> x==" Note_on_c", s[s[:,1].==i,3]) #initial time
                fin = findlast(x-> x==" Note_off_c", s[s[:,1].==i,3]) #finish time
                if ini == 0 || fin == 0; continue; end
                mat = s[s[:,1].==i,:]
                mat = mat[findall(x->x!="", mat[:,5]),:]
                voces[i] = get_onoff_notes(mat) #construct an array of information of initial time, finish time, pitch and intensity
            end
        end
    end
    if nv > 1
        filter_undef!(voces)
        filter!(x->length(x)>0,voces)
        n_vs = []
        nv = length(voces)
        for i=1:nv
            m_v = find_more_voices(voces[i])
            for j =1:length(m_v)
                push!(n_vs, m_v[j])
            end
        end
    else
        try
            n_vs = []
            m_v = find_more_voices(voces[1])
            for j =1:length(m_v)
                push!(n_vs, m_v[j])
            end
        catch
            n_vs = []
            m_v = find_more_voices(voces)
            for j =1:length(m_v)
                push!(n_vs, m_v[j])
            end
        end
    end
    filter!(x->length(x)>0,n_vs)
    nv = size(n_vs)[1]
    t_max = max_tempo(n_vs, nv)
    window_size, w_res = divrem(t_max, n_chunks)
    chunked_voices = Array{Array{Matrix}}(undef,nv) #initialze an array
    for v = 1:nv
        chunked_voice = []
        for i = 1:n_chunks
            st =  (i-1) * window_size + 1
            en = st + window_size - 1
            in_notes = findall(x-> st< x < en, n_vs[v][:,1])
            push!(chunked_voice, n_vs[v][in_notes,:])
        end
        chunked_voices[v] = chunked_voice
    end
    return chunked_voices, mq
end

function chunk_voices_ws(s, window_size)
    mq = s[1,6] #quarter of a note
    nv = s[findlast(x-> x==" Note_on_c",s[:,3]),1] #- 1 #estimates how many voices are in the midi
    voces = Array{Matrix}(undef,nv) #initialze an array
    #next function is to get the time in miliseconds when the pitch starts and ends.
    if nv ==1
        #nv = 1
        if isempty(findall(x-> x==" Note_off_c",s[:,3]))
            b = s[s[:,1].==1,:] #takes the events of the voice i
            mat = b[b[:,3].==" Note_on_c", :]
            mat = mat[findall(x->x!="", mat[:,5]),:]
            voces = get_onon_notes(mat) #construct an array of of information of initial time, finish time, pitch and intensity.
        else
            mat = s[s[:,1].==1,:]
            mat = mat[findall(x->x!="", mat[:,5]),:]
            voces = get_onoff_notes(mat)
        end
    else
        voces = Array{Matrix}(undef,nv) #initialze an array
        if isempty(findall(x-> x==" Note_off_c",s[:,3]))
            #checks if the midi has events of note_off
            for i = 1:(nv) #if does not, it construct the series in this way
                b = s[s[:,1].==i,:] #takes the events of the voice i
                if size(b[b[:,3].==" Note_on_c", :])[1] == 0; continue; end #checks if there are notes in the channel
                mat = b[b[:,3].==" Note_on_c", :]
                mat = mat[findall(x->x!="", mat[:,5]),:]
                voces[i] =  get_onon_notes(mat)#construct an array of information of initial time, finish time, pitch and intensity.
            end
        else
            for i = 1:nv #if has note_off events , it does this way
                ini = findfirst(x-> x==" Note_on_c", s[s[:,1].==i,3]) #initial time
                fin = findlast(x-> x==" Note_off_c", s[s[:,1].==i,3]) #finish time
                if ini == 0 || fin == 0; continue; end
                mat = s[s[:,1].==i,:]
                mat = mat[findall(x->x!="", mat[:,5]),:]
                voces[i] = get_onoff_notes(mat) #construct an array of information of initial time, finish time, pitch and intensity
            end
        end
    end
    if nv > 1
        filter_undef!(voces)
        filter!(x->length(x)>0,voces)
        n_vs = []
        nv = length(voces)
        for i=1:nv
            m_v = find_more_voices(voces[i])
            for j =1:length(m_v)
                push!(n_vs, m_v[j])
            end
        end
    else
        try
            n_vs = []
            m_v = find_more_voices(voces[1])
            for j =1:length(m_v)
                push!(n_vs, m_v[j])
            end
        catch
            n_vs = []
            m_v = find_more_voices(voces)
            for j =1:length(m_v)
                push!(n_vs, m_v[j])
            end
        end
    end
    filter!(x->length(x)>0,n_vs)
    nv = size(n_vs)[1]
    t_max = max_tempo(n_vs, nv)
    #window_size, w_res = divrem(t_max, n_chunks)
    n_chunks, w_res = divrem(t_max, window_size)
    chunked_voices = Array{Array{Matrix}}(undef,nv) #initialze an array
    for v = 1:nv
        chunked_voice = []
        for i = 1:n_chunks
            st =  (i-1) * window_size + 1
            en = st + window_size - 1
            in_notes = findall(x-> st< x < en, n_vs[v][:,1])
            push!(chunked_voice, n_vs[v][in_notes,:])
        end
        chunked_voices[v] = chunked_voice
    end
    return chunked_voices
end
function chunk_voices_measure(s)
    frac = s[findall(x-> x == " Time_signature", s[:,3]),2:5]  #this is the measure of the piece in fractional way a / b
    if length(frac) > 4
        measure = map((x,y) -> x / 2^y ,frac[:,3], frac[:,4])
        meas_time = frac[:,1]
    else
        measure = frac[3] / 2^frac[4]
    end

    q_t = s[1,6] #the time in ms for a quarter of note 1/4
    w_s = measure * q_t / 0.25 #this is the miliseconds of a measure
    nv = s[findlast(x-> x==" Note_on_c",s[:,3]),1] #- 1 #estimates how many voices are in the midi
    voces = Array{Matrix}(undef,nv) #initialze an array
    #next function is to get the time in miliseconds when the pitch starts and ends.
    if nv ==1
        #nv = 1
        if isempty(findall(x-> x==" Note_off_c",s[:,3]))
            b = s[s[:,1].==1,:] #takes the events of the voice i
            mat = b[b[:,3].==" Note_on_c", :]
            mat = mat[findall(x->x!="", mat[:,5]),:]
            voces = get_onon_notes(mat) #construct an array of of information of initial time, finish time, pitch and intensity.
        else
            mat = s[s[:,1].==1,:]
            mat = mat[findall(x->x!="", mat[:,5]),:]
            voces = get_onoff_notes(mat)
        end
    else
        voces = Array{Matrix}(undef,nv) #initialze an array
        if isempty(findall(x-> x==" Note_off_c",s[:,3]))
            #checks if the midi has events of note_off
            for i = 1:(nv) #if does not, it construct the series in this way
                b = s[s[:,1].==i,:] #takes the events of the voice i
                if size(b[b[:,3].==" Note_on_c", :])[1] == 0; continue; end #checks if there are notes in the channel
                mat = b[b[:,3].==" Note_on_c", :]
                mat = mat[findall(x->x!="", mat[:,5]),:]
                voces[i] =  get_onon_notes(mat)#construct an array of information of initial time, finish time, pitch and intensity.
            end
        else
            for i = 1:nv #if has note_off events , it does this way
                ini = findfirst(x-> x==" Note_on_c", s[s[:,1].==i,3]) #initial time
                fin = findlast(x-> x==" Note_off_c", s[s[:,1].==i,3]) #finish time
                if ini == 0 || fin == 0; continue; end
                mat = s[s[:,1].==i,:]
                mat = mat[findall(x->x!="", mat[:,5]),:]
                voces[i] = get_onoff_notes(mat) #construct an array of information of initial time, finish time, pitch and intensity
            end
        end
    end
    # if nv > 1
    #     filter_undef!(voces)
    #     filter!(x->length(x)>0,voces)
    #     n_vs = []
    #     nv = length(voces)
    #     for i=1:nv
    #         m_v = find_more_voices(voces[i])
    #         for j =1:length(m_v)
    #             push!(n_vs, m_v[j])
    #         end
    #     end
    # else
    #     try
    #         n_vs = []
    #         m_v = find_more_voices(voces[1])
    #         for j =1:length(m_v)
    #             push!(n_vs, m_v[j])
    #         end
    #     catch
    #         n_vs = []
    #         m_v = find_more_voices(voces)
    #         for j =1:length(m_v)
    #             push!(n_vs, m_v[j])
    #         end
    #     end
    # end
    if nv != 1
        n_vs = filter_undef!(voces)
        filter!(x->length(x)>0,n_vs)
    else
        n_vs = []
        push!(n_vs,voces)
    end
    #filter!(x->length(x)>0,n_vs)
    nv = size(n_vs)[1]
    t_max = max_tempo(n_vs, nv)
    #window_size, w_res = divrem(t_max, n_chunks)
    number_measures = length(measure)
    if number_measures > 1
        chunked_voices = Array{Array{Matrix}}(undef,nv) #initialze an array
        for v = 1:nv
            chunked_voice = []
            for ct = 1:number_measures
                chunked_meas = []
                if ct < number_measures
                    t_final = meas_time[ct+1] - meas_time[ct]
                else
                    t_final = t_max - meas_time[ct]
                end
                n_chunks, w_res = divrem(t_final,w_s[ct])
                for i = 1:n_chunks
                    st =  meas_time[ct] + (i-1) * w_s[ct]
                    en = st + w_s[ct] - 1
                    in_notes = findall(x-> st<= x < en, n_vs[v][:,1])
                    push!(chunked_meas, n_vs[v][in_notes,:])
                end
                push!(chunked_voice, chunked_meas)
            end
            #println(chunked_voice)
            chunked_voices[v] = vcat(chunked_voice...)
        end
    else
        n_chunks, w_res = divrem(t_max, w_s)
        chunked_voices = Array{Array{Matrix}}(undef,nv) #initialze an array
        for v = 1:nv
            chunked_voice = []
            for i = 1:n_chunks
                st =  (i-1) * w_s + 1
                en = st + w_s - 1
                in_notes = findall(x-> st<= x < en, n_vs[v][:,1])
                push!(chunked_voice, n_vs[v][in_notes,:])
            end
            chunked_voices[v] = chunked_voice
        end
    end
    return chunked_voices
end

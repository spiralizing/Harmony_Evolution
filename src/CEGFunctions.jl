#CENTER OF EFFECT ---- KEY FINDING ALGORITHM
#This set of functions corresponds to the spiral representations
#writing the names of the keys, ordered by fifths
#WARNING!!!!
#SOME FUNCTIONS NEED THE PACKAGE PyCall AND NEED THE LIBRARY music21 IN THE PYTHON Python_path SEE BELOW
using PyCall
m21 = pyimport("music21") 
numpy = pyimport("numpy")
math = pyimport("math")
####
###VARIABLES#####
pitch_names = ["C","G","D","A","E","B","Gb/F#","Db","Ab","Eb","Bb","F"]
nMajor_keys = ["C", "G", "D", "A","E", "B/Cb","Gb/F#","Db/C#","Ab","Eb", "Bb","F"]
nminor_keys = ["c","g","d","a","e","b","f#","c#","g#/ab","eb/d#","bb","f"]
nall_keys = [] #All keys
for (a,b) in zip(nminor_keys,nMajor_keys)
    push!(nall_keys,a); push!(nall_keys,b)
end
cf_notes = [0, 7, 2, 9, 4, 11, 6, 1, 8, 3, 10, 5]
midi_notes = ["C","C#","D","Eb","E","F","F#","G","G#","A","Bb","B"]
#the real value of the notes in the circle of fifths, chromatic scale representation from C = 0 to B = 11
# r = 1
# h = sqrt(2/15)
# w = [0.536, 0.274, 0.19] #weights used for the center of effect algorithm
# u = v = o = w
# a = 0.75
# b = 0.75
h_octav = 4.3818
"""
    get_pitch(k; r=1, h=sqrt(2/15))

    Returns the location (x,y,z) of a pitch "k"
    in the spiral representation, pitches are in module 12
    the variable "r" is the radius of the cylinder
    and "h" is the vertical distance between thirds.

"""
function get_pitch(k; r=1, h=sqrt(2/15))
    p = [r*sin(k*pi/2), r*cos(k*pi/2), k*h]
    return p
end
"""
    get_Major_chord(k; w=[0.536, 0.274, 0.19], r=1, h=sqrt(2/15))

    Returns the location (x,y,z) of a Major chord given by a linear combination
    of a triad constructed from a fundamental note "k", the fifth (k+1) and the third (k+4)
    in the spiral representation. The vector "w" contains the weights used for each term.
"""
function get_Major_chord(k; w=[0.536, 0.274, 0.19], r=1, h=sqrt(2/15))
    CM = w[1]*get_pitch(k, r=r, h=h) + w[2]*get_pitch(k+1,r=r,h=h) + w[3]*get_pitch(k+4,r=r,h=h)
    return CM
end

"""
    get_minor_chord(k; u=[0.536, 0.274, 0.19], r=1, h=sqrt(2/15))

    Returns the location (x,y,z) of a Minor chord given by a linear combination
    of a triad constructed from a fundamental note "k", the fifth (k+1) and the minor third (k-3)
    in the spiral representation. The vector "u" contains the weights used for each term.
"""
function get_minor_chord(k; u=[0.536, 0.274, 0.19], r=1, h=sqrt(2/15))
    Cm = u[1]*get_pitch(k, r=r, h=h) + u[2]*get_pitch(k+1,r=r,h=h) + u[3]*get_pitch(k-3,r=r,h=h)
    return Cm
end

"""
    get_Major_key(k; o=[0.536, 0.274, 0.19], w=[0.536, 0.274, 0.19], r=1, h=sqrt(2/15))

    Returns the location (x,y,z) of a Major Key given by a linear combination
    of chord triads constructed from the major chords of a fundamental "tonic" note "k",
    the dominant/fifth (k+1) and the subdominant/minor third (k-3) in the spiral representation.
    The vector "o" contains the weights used for each chord and "w" are the weights on each note in the chord.
"""
function get_Major_key(k; o=[0.536, 0.274, 0.19], w=[0.536, 0.274, 0.19], r=1, h=sqrt(2/15))
    TM = o[1]*get_Major_chord(k, w=w, r=r, h=h) + o[2]*get_Major_chord(k+1, w=w, r=r, h=h) + o[3]*get_Major_chord(k-1, w=w, r=r, h=h)
    return TM
end

"""
    get_minor_key(k; v=[0.536, 0.274, 0.19], w=[0.536, 0.274, 0.19],
    u=[0.536, 0.274, 0.19], a=0.75, b=0.75, r=1, h=sqrt(2/15))

    Returns the location (x,y,z) of a Minor Key given by a linear combination
    of chord triads constructed from the major chords of a fundamental note "k",
    the possible "dominant" triad (k+1) and the possible "subdominant" triad (k-3) in the spiral representation.
    The vector "v" contains the weights used for each chord and "w" are the weights on each note in the major chords.
    "u" are the weights for the minor chords. The dominant and subdominant minor/major chords are weighted by "a" and "b".
"""
function get_minor_key(k; v=[0.536, 0.274, 0.19], w=[0.536, 0.274, 0.19], u=[0.536, 0.274, 0.19], a=0.75, b=0.75, r=1, h=sqrt(2/15))
    Tm = v[1]*get_minor_chord(k, u=u, r=r, h=h) + v[2]*(a*get_Major_chord(k+1, w=w, r=r, h=h) + (1-a)*get_minor_chord(k+1, u=u, r=r, h=h)) + v[3]*(b*get_minor_chord(k-1, u=u, r=r, h=h) + (1-b)*get_Major_chord(k-1, w=w, r=r, h=h))
    return Tm
end

"""
    get_cfpitch_mod12(pitch_seq)

    Converts a pitch sequence of MIDI values (0-127) into a pitch sequence in mod12 ordered by fifths
    to be consistent with the cylindrical representation, starting in C=0, G=1, ..., etc.
"""
function get_cfpitch_mod12(pitch_seq)
    m12v = convert(Array{Int64,1},map(x->mod(x,12),pitch_seq))
    return map(y->findfirst(x->x==y,cf_notes)-1,m12v)
end

"""
    get_cfpitch(pitch_seq)

    Converts a pitch sequence of MIDI values (0-127) into a pitch sequence ordered by fifths
    to be consistent with the cylindrical representation, starting in C=0, G=1, ..., etc.
"""
function get_cfpitch(pitch_seq)
    return map(y->findfirst(x->x==y,all_cf_notes)-1,pitch_seq)
end

function isminor(s::Any) 
    return occursin(r"^[a-z,#,/,\s]+$",s)
end
function isminor(k_seq::Array{Any,1})
    gkey = get_rank_freq(k_seq)[1,1]
    return occursin(r"^[a-z,#,/,\s]+$",gkey)
end
#functional harmony sequence; sequence of the relative key from a given global key in roman numerals.
"""
    funhar_seq(kseq, fun_key)
    Returns a sequence of Roman numerals given a reference for the Tonic Key.
"""
function funhar_seq(kseq, fun_key)
    fh_kseq = []
    if isminor(fun_key)
        nk = findfirst(x-> x==fun_key, nminor_keys)
    else
        nk = findfirst(x-> x==fun_key, nMajor_keys)
    end
    minork = circshift(nminor_keys, 13-nk) 
    majork = circshift(nMajor_keys, 13-nk)
    for k = 1:length(kseq)
        if isminor(kseq[k])
            pos = findfirst(x-> x==kseq[k],minork) - 1
            push!(fh_kseq, Minor_RN[pos])
        else
            pos = findfirst(x-> x== kseq[k],majork) - 1
            push!(fh_kseq, Major_RN[pos])
        end
    end
    return fh_kseq
end
##--

#handwritting the functional harmony notations "key or chord relative to the fundamental"
#Major roman numerals, for the coe notation.
Major_RN = Dict(0 => "I", 
    7 =>"#I/bII",
    2 => "II",
    9 => "#II/bIII",
    4 => "III",
    #"III#/IVb" = 11
    11 => "IV",
    6 => "#IV/bV",
    1 => "V",
    8 => "#V/bVI",
    3 => "VI",
    10 => "#VI/bVII",
    5 => "VII",
   )
##--
Minor_RN = Dict(0 => "i", 
    7 =>"#i/bii",
    2 => "ii",
    9 => "#ii/biii",
    4 => "iii",
    #"III#/IVb" = 11
    11 => "iv",
    6 => "#iv/bv",
    1 => "v",
    8 => "#v/bvi",
    3 => "vi",
    10 => "#vi/bvii",
    5 => "vii",
   )
#defining all locations for the array for the "circle of fifths" going up and down, taking C as k=0
notes = [i for i=0:127]
pitches = map(x-> get_pitch(x), notes)
major_chords = map(x-> get_Major_chord(x), notes)
minor_chords = map(x-> get_minor_chord(x), notes)
major_keys = map(x-> get_Major_key(x), notes)
minor_keys = map(x-> get_minor_key(x), notes)
pos_all_keys = vcat(major_keys, minor_keys)
#Defining all circle of fifth notes and names of major and minor keys
midi_note_names = []
all_cf_notes = []
all_Major_keys = []
all_minor_keys = []
all_pitch_names = []
push!(all_cf_notes, cf_notes)
push!(all_Major_keys, nMajor_keys)
push!(all_minor_keys, nminor_keys)
push!(midi_note_names, midi_notes)
push!(all_pitch_names, pitch_names)

for i = 1:9
    push!(all_cf_notes, cf_notes.+(12*i))
    push!(all_Major_keys, nMajor_keys)
    push!(all_minor_keys, nminor_keys)
    push!(midi_note_names, midi_notes)
    push!(all_pitch_names, pitch_names)
end
push!(all_cf_notes, cf_notes[1:8].+(120))
push!(all_Major_keys, nMajor_keys[1:8])
push!(all_minor_keys, nminor_keys[1:8])
push!(midi_note_names, midi_notes[1:8])
push!(all_pitch_names, pitch_names[1:8])


all_cf_notes = vcat(all_cf_notes...)
all_Major_keys = vcat(all_Major_keys...)
all_minor_keys = vcat(all_minor_keys...)
all_keys = vcat(all_Major_keys,all_minor_keys)
midi_note_names = vcat(midi_note_names...)
all_pitch_names = vcat(all_pitch_names...)

"""
    function get_coe_path(notes_chunk)

    Returns the path of the center of effect for a given set of notes.
"""

function get_coe_path(notes_chunk; r=1, h=sqrt(2/15), mod_12=false,all_keys=all_keys, pos_all_keys=pos_all_keys, sbeat_w=[[1.],[1.]], lin_w=1)
    ptcs = notes_chunk[:,6] #Pitches
    durs = notes_chunk[:,5] #durations
    pbeat = notes_chunk[:,1] #beat where the note starts
    beat_w = ones(length(durs)) #array of the beat weights
    for b = 1:length(sbeat_w[1])
        loc_b = findall(x-> x==sbeat_w[1][b], pbeat) #finding all notes that start at beat sbeat_w[1][b]
        beat_w[loc_b] .= sbeat_w[2][b] #this is the weight.
    end
    notas, n_we = get_local_lin_w(ptcs, lin_w) #doing the linear weight in the pitches
    ii = vcat(map(x-> findall(y-> y==x, notas), ptcs)...)
    #println(ptcs)
    b_wei = n_we[ii] #getting the linear weight for every note i n the array of pitches
    if mod_12
        mod12_seq = get_cfpitch_mod12(ptcs)
        spi_ix = reorder_seq_closest(mod12_seq).+12
        #println("WARNING! \n Doing module 12 notes.")
    else
        spi_ix = get_cfpitch(ptcs)
        if length(spi_ix) > 1
            max_dif = 10
            while max_dif > 6
                nex_seq = reorder_seq_closest(spi_ix)
                max_dif = maximum(map(x-> abs(x),diff(nex_seq)))
                spi_ix = nex_seq
            end
        end
    end
    spi_p = map(x-> get_pitch(x, r=r, h=h), spi_ix) #getting the location (x,y,z) for each pitch
    t_ws = map((x,y,z)-> x*y*z, beat_w,durs, b_wei) #computing the total weights
    
    cv_i = []
    #cv_i = map((x,y)-> x*y, t_ws, spi_p) / sum(t_ws) #computing the location of the pitches with their relative weights
    for i in 1:length(spi_p)
        push!(cv_i, map((x,y)-> x*y, t_ws[1:i],spi_p[1:i])/sum(t_ws[1:i]))  
    #c_i = sum(cv_i) #finding the center of effect
    end
    return map(x-> sum(x),cv_i)
end

"""
    function get_piece_by_measure(data; csv=true, qdiv=32, abs_time=true, start_measure=1)

    Returns 
"""

function get_piece_by_measure(data; csv=true, qdiv=32, abs_time=true, start_measure=1)
    if csv
        return get_piece_by_measure_csv(data, qdiv=qdiv, abs_time=abs_time)
    else
        return get_piece_by_measure_m21(data, start_measure=start_measure)
    end
end

"""
    get_piece_by_measure_csv(s; qdiv=32, abs_time=true)

    Returns a vector of two dimensional arrays of the basic information of the piece for analyzing its properties,
    input "s" should be a CSV table from a CSV file converted from a MIDI file with the midicsv script available for free
    in: https://www.fourmilab.ch/webtools/midicsv/ .

    The first output consist of a vector which components represent each measure of the music piece,
    each measure contains four columns; Beat: number of beat where the note starts (beat or halfbeat)
    Measure: the fraction of the measure (n/d), Duration: the duration of the note in midiclocks, if abs_time=true
    the output returns the time in ms where the note started and ended before the duration of the note,
    Pitch: the pitch of the note represented in MIDI notation (0-127).

    The second output corresponds to the total number of notes in the piece and the number of notes that
    fall outside the threshold given (1/qdiv), most of the time if notes fall outside this threshold is because
    the MIDI was sequenced (recorded) instead of generated with a music score software.
"""
function get_piece_by_measure_csv(s; qdiv=32, abs_time=true)
    frac = s[findall(x-> x == " Time_signature", s[:,3]),2:5]  #this is the measure of the piece in fractional way a / b
    if length(frac) > 4
        bot = map(x ->  2^x ,frac[:,4])
        sym_measure = map((x,y) -> join([x y],"/"), frac[:,3], bot)
        measure = map((x,y) -> x / 2^y ,frac[:,3], frac[:,4])
        meas_time = frac[:,1]
        numer = frac[:,3]
    else
        measure = frac[3] / 2^frac[4]
        numer = []
        push!(numer,frac[3])
        sym_measure = join([frac[3] 2^frac[4]], "/")
    end

    q_t = s[1,6] #the time in ms for a quarter of note 1/4
    th = q_t / qdiv #64nd of tolerance with qdiv=32
    w_s = measure * q_t / 0.25 #this is the miliseconds of a measure
    smdiv = 32 #this is the threshold for the tolerance in the beats
    nv = s[findlast(x-> x==" Note_on_c",s[:,3]),1] #- 1 #estimates how many voices are in the midi
    voces = Array{Matrix}(undef,nv) #initialze an array
    subdiv = map((x,y)-> x/y, w_s, numer)
    notes_measure = []
    for i = 1:length(measure)
        push!(notes_measure,cumsum([subdiv[i] for j=1:numer[i]]))
    end
    n_off = findall(x-> x==" Note_off_c",s[:,3])
    n_on = findall(x-> x==" Note_on_c",s[:,3])
    #next lines are to get the time in miliseconds when the pitch starts and ends.
    if nv ==1
        #nv = 1
        if isempty(n_off) || length(n_off) < length(n_on)/2
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
        if isempty(n_off) || length(n_off) < length(n_on)/2
            #checks if the midi has events of note_off
            for i = 1:(nv) #if does not, it construct the series in this way
                b = s[s[:,1].==i,:] #takes the events of the voice i
                #if size(b[b[:,3].==" Note_on_c", :])[1] == 0; continue; end #checks if there are notes in the channel
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
    if nv != 1
        n_vs = filter_undef!(voces)
        filter!(x->length(x)>0,n_vs)
    else
        n_vs = []
        push!(n_vs,voces)
    end
    nv = size(n_vs)[1]
    t_max = max_tempo(n_vs, nv)
    number_measures = length(measure)
    all_v = vcat(n_vs...)
    norder = sortperm(all_v[:,1])
    all_v = all_v[norder,:]
    all_divis = []
    if number_measures > 1
        chunked_voice = []
        for ct = 1:number_measures
            th = (q_t / qdiv) / notes_measure[ct][1]
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
                in_notes = findall(x-> st<= x < en, all_v[:,1])
                st_times = all_v[in_notes,1]
                in_times = all_v[in_notes,1] .- st
                loc = zeros(length(in_notes))
                num_n = [sym_measure[ct] for i = 1:length(in_notes)]
                dur = all_v[in_notes,2] - all_v[in_notes,1]
                pitc = all_v[in_notes,3]
                for z = 1:length(in_notes)
                    divis = in_times[z] / notes_measure[ct][1]
                    dif, bt = modf(divis)
                    if dif <= th #check if the note starts in the tempo is found, within a threshold
                        loc[z] = bt + 1
                    elseif 1 - dif <= th #if the note starts in the next tempo
                        loc[z] = bt + 2
                    elseif dif  >= 0.5 - th && dif <= 0.5 + th  #check if the note starts on a half beat
                        loc[z] = bt + 1.5
                    else
                        loc[z] = 0 #the note doesn't start in the tempo, is in between.
                    end
                    push!(all_divis, dif*smdiv)
                    #println("Pitch: $(pitc[z])",'\t',"Starting time in piece (midi clock): $(st_times[z]), in the measure:$(in_times[z]) ",'\t',"Ratio Starting/LengthBeat: $(divis)",'\t',"Selected starting beat (in measure): $(loc[z])")
                end
                if abs_time
                    push!(chunked_meas, [loc num_n all_v[in_notes,1:2] dur all_v[in_notes,3]])
                else
                    push!(chunked_meas, [loc num_n all_v[in_notes,3] dur])
                end
            end
            #next lines are for the last measure or the last fraction left.

            if ct==number_measures && w_res > 0
                st = meas_time[ct] + n_chunks*w_s[ct]
                en = t_max
                in_notes = findall(x-> st<= x < en, all_v[:,1])
                st_times = all_v[in_notes,1]
                in_times = all_v[in_notes,1] .- st
                loc = zeros(length(in_notes))
                num_n = [sym_measure for i = 1:length(in_notes)]
                dur = all_v[in_notes,2] - all_v[in_notes,1]
                pitc = all_v[in_notes,3]
                for z = 1:length(in_notes)
                    #pl = findfirst(x-> x== in_times[z], sort([in_times[z];notes_measure[ct]]))
                    divis = in_times[z] / notes_measure[1][1]
                    dif, bt = modf(divis) #- (notes_measure[ct][pl] - notes_measure[ct][1])
                    #hbeat = notes_measure[ct][1]/2
                    if dif <= th #check if the note starts in the tempo is found, within a threshold
                        loc[z] = bt + 1
                    elseif 1 - dif <= th #if the note starts in the next tempo
                        loc[z] = bt + 2
                    elseif dif  >= 0.5 - th && dif <= 0.5 + th  #check if the note starts on a half beat
                        loc[z] = bt + 1.5
                    else
                        loc[z] = 0 #the note doesn't start in the tempo, is in between.
                    end
                    push!(all_divis, dif*smdiv)
                    #println("Pitch: $(pitc[z])",'\t',"Starting time in piece (midi clock): $(st_times[z]), in the measure:$(in_times[z]) ",'\t',"Ratio Starting/LengthBeat: $(divis)",'\t', "Ratio Decimals multiplied by 32: $(dif*smdiv)")#"Selected starting beat (in measure): $(loc[z])")
                end
                if abs_time
                    push!(chunked_meas, [loc num_n all_v[in_notes,1:2] dur all_v[in_notes,3]]) #returns also the values for the absolute time in ms when the note started and when it ended
                else
                    push!(chunked_meas, [loc num_n all_v[in_notes,3] dur])
                end
            end
            push!(chunked_voice, chunked_meas)
        end
    else
        n_chunks, w_res = divrem(t_max, w_s)
        th = (q_t / qdiv) / notes_measure[1][1]
        chunked_voice = []
        for i = 1:n_chunks
            st =  (i-1) * w_s
            en = st + w_s - 1
            in_notes = findall(x-> st<= x < en, all_v[:,1])
            st_times = all_v[in_notes,1]
            in_times = all_v[in_notes,1] .- st
            loc = zeros(length(in_notes))
            num_n = [sym_measure for i = 1:length(in_notes)]
            dur = all_v[in_notes,2] - all_v[in_notes,1]
            pitc = all_v[in_notes,3]
            for z = 1:length(in_notes)
                divis = in_times[z] / notes_measure[1][1]
                dif, bt = modf(divis)
                if dif <= th #check if the note starts in the tempo is found, within a threshold
                    loc[z] = bt + 1
                elseif 1 - dif <= th #if the note starts in the next tempo
                    loc[z] = bt + 2
                elseif dif  >= 0.5 - th && dif <= 0.5 + th  #check if the note starts on a half beat
                    loc[z] = bt + 1.5
                else
                    loc[z] = 0 #the note doesn't start in the tempo, is in between.
                end
                push!(all_divis, dif*smdiv)
                #println("Pitch: $(pitc[z])",'\t',"Starting time in piece (midi clock): $(st_times[z]), in the measure:$(in_times[z]) ",'\t',"Ratio Starting/LengthBeat: $(divis)",'\t',"Selected starting beat (in measure): $(loc[z])")
            end
            if abs_time
                push!(chunked_voice,[loc num_n all_v[in_notes,1:2] dur all_v[in_notes,3]])
            else
                push!(chunked_voice,[loc num_n all_v[in_notes,3] dur])
            end
                
        end
        if w_res > 0
            st = n_chunks*w_s
            if st != t_max
                en = t_max
                #println(st,"\t",en)
                in_notes = findall(x-> st<= x < en, all_v[:,1])
                st_times = all_v[in_notes,1]
                in_times = all_v[in_notes,1] .- st
                loc = zeros(length(in_notes))
                num_n = [sym_measure for i = 1:length(in_notes)]
                dur = all_v[in_notes,2] - all_v[in_notes,1]
                pitc = all_v[in_notes,3]
                for z = 1:length(in_notes)
                    #pl = findfirst(x-> x== in_times[z], sort([in_times[z];notes_measure[ct]]))
                    divis = in_times[z] / notes_measure[1][1]
                    dif, bt = modf(divis) #- (notes_measure[ct][pl] - notes_measure[ct][1])
                    #hbeat = notes_measure[ct][1]/2
                    if dif <= th #check if the note starts in the tempo is found, within a threshold
                        loc[z] = bt + 1
                    elseif 1 - dif <= th #if the note starts in the next tempo
                        loc[z] = bt + 2
                    elseif dif  >= 0.5 - th && dif <= 0.5 + th  #check if the note starts on a half beat
                        loc[z] = bt + 1.5
                    else
                        loc[z] = 0 #the note doesn't start in the tempo, is in between.
                    end
                    push!(all_divis, dif*smdiv)
                    #println("Pitch: $(pitc[z])",'\t',"Starting time in piece (midi clock): $(st_times[z]), in the measure:$(in_times[z]) ",'\t',"Ratio Starting/LengthBeat: $(divis)",'\t', "Ratio Decimals multiplied by 32: $(dif*smdiv)")#"Selected starting beat (in measure): $(loc[z])")
                end
                if abs_time
                    push!(chunked_voice,[loc num_n all_v[in_notes,1:2] dur all_v[in_notes,3]])
                else
                    push!(chunked_voice,[loc num_n all_v[in_notes,3] dur])
                end
            end
        end
    end
    #next lines are for the last measure or the last fraction left.
    residues = map(x-> modf(x)[1], all_divis)
    n_off = length(findall(x-> x > 10^-6, residues))
    res_frac = n_off / length(residues)
    if number_measures > 1
        out_measures = vcat(chunked_voice...)
    else
        out_measures = chunked_voice
    end
    if res_frac > 0.05
        println("WARNING!!! \n THE FRACTION OF NOTES FALLING OUTSIDE THE MIDI CLOCK IS: $(res_frac)")
    end
    return filter(x-> !isempty(x),out_measures), [length(residues) n_off]
end
"""
    function get_piece_by_measure_m21(m21_data; start_measure=1)

    Returns an array with the information of pitches and duration on each measure in the following format:

    | #measure | time signature | start_quarter | end_quarter | duration (in quarters) | pitch (0-127) |
"""
function get_piece_by_measure_m21(m21_data; start_measure=1)

    notes = []
    for parts in m21_data.parts #getting all the notes in all parts, with duration and pitch
        for note in parts.flat.notes
            try
                local start, duration
                measure = note.measureNumber + (start_measure - 1)
                start = float(note.offset)
                duration = float(note.quarterLength)
                pitch = note.pitch.ps
                #print(start,'\n')
                push!(notes, [measure start start + duration duration pitch])
            catch
                local start, duration
                measure = note.measureNumber + (start_measure - 1)
                try

                    start = float(note.offset.())
                catch
                    try
                        numer = note.offset.numerator
                        denom = note.offset.denominator
                        start = float(numer / denom)
                    catch
                        start = float(note.offset)
                    end
                end
                #print(note.quarterLength,'\n')
                try
                    #print("here 1 \n")
                    duration = float(note.quarterLength.())
                catch
                    try
                        #print("here 2 \n")
                        numer = note.quarterLength.numerator
                        denom = note.quarterLength.denominator
                        duration = float(numer / denom)
                    catch
                        #print("here 3 \n")
                        duration = float(note.quarterLength)

                    end
                end
                #print("duration before appending: ",duration, '\n')
                #print("start before appending: ", start, '\n')
                for chord_note in note.pitches
                    pitch = chord_note.ps
                    #print("aqui \n")
                    #print(start,'\t',duration,'\n')
                    push!(notes, [measure start start + duration duration pitch])

                end
            end
        end
    end
    all_notes = vcat(notes...)
    n_measures = Int(maximum(all_notes[:,1]))
    timesig = []
    lquarters = []
    for i in 0:n_measures #extracting time signatures
        try
            ts = m21_data.parts[1].measure(i).timeSignature
            push!(timesig,[i ts.ratioString])
            nquarters = ts.numerator / ts.denominator * 4
            push!(lquarters, [i nquarters])
        catch
        end
    end

    measures_piece = []
    if length(timesig) > 1
        for n_ts in 2:length(timesig)
            last_m = timesig[n_ts][end,1] #get the last measure where the time signature applies 
            first_m = timesig[n_ts-1][1,1]
            ix = findall(x-> first_m<=x<=last_m, all_notes[:,1])
            m_ts = [timesig[n_ts-1][1,2] for i in 1:length(ix)]
            push!(measures_piece, [all_notes[ix,1] m_ts all_notes[ix, 2:end]])
        end
        #now last time signature
        last_m = n_measures
        first_m = timesig[end][1,1]
        ix = findall(x-> first_m<=x<=last_m, all_notes[:,1])
        m_ts = [timesig[end][1,2] for i in 1:length(ix)]
        push!(measures_piece, [all_notes[ix,1] m_ts all_notes[ix, 2:end]])
    else
        ix = all_notes[:,1]
        m_ts = [timesig[1][1,2] for i in 1:length(ix)]
        push!(measures_piece, [all_notes[:,1] m_ts all_notes[:,2:end]])
    end
#=     if length(timesig) > 1 #if there is more than one time signature
        for n_ts in 2:length(timesig)
            last_m = timesig[n_ts][end,1] #get the last measure where the time signature applies 
            first_m = timesig[n_ts-1][1,1]
            ix = findall(x-> first_m<=x<=last_m, all_notes[:,1])
            m_ts = [timesig[n_ts-1][1,2] for i in 1:length(ix)]
            ini_time = map(x-> mod(x, lquarters[n_ts-1][1,2]+1), all_notes[ix,2])
            fin_time = map(x-> mod(x, lquarters[n_ts-1][1,2]+1), all_notes[ix,3])
            push!(measures_piece, [all_notes[ix,1] m_ts ini_time fin_time all_notes[ix,4:end]])
        end
        #now last time signature
        last_m = n_measures
        first_m = timesig[end][1,1]
        ix = findall(x-> first_m<=x<=last_m, all_notes[:,1])
        m_ts = [timesig[end][1,2] for i in 1:length(ix)]
        ini_time = map(x-> mod(x, lquarters[end][1,2]+1), all_notes[ix,2])
        fin_time = map(x-> mod(x, lquarters[end][1,2]+1), all_notes[ix,3])
        push!(measures_piece, [all_notes[ix,1] m_ts ini_time fin_time all_notes[ix,4:end]])

    else
        ix = all_notes[:,1]
        m_ts = [timesig[1][1,2] for i in 1:length(ix)]
        ini_time = map(x-> mod(x, lquarters[1][1,2]+1), all_notes[:,2])
        fin_time = map(x-> mod(x, lquarters[1][1,2]+1), all_notes[:,3])
        push!(measures_piece, [all_notes[:,1] m_ts ini_time fin_time all_notes[:,4:end]])
    end =#
    out_measures = vcat(measures_piece...)
    return out_measures[sortperm(out_measures[:,1]),:]
end
"""
    get_key_sequence(notes_chunk; r=1, h=sqrt(2/15), mod_12=false,
        all_keys=all_keys, pos_all_keys=pos_all_keys, sbeat_w=[[1.],[1.]], lin_w=1)

    Returns the center of effect "c_i(x,y,z)" of a given sequence of notes in the format
    | beat | measure | duration | pitch |

    The second output is a list of the most likely keys with their respective distance
    to the Center of Effect and ranked in ascending order.

    Third output returns a table of the values for the computation in the format:
    | beat | measure | duration | pitch | pitch name | CF pitch (spiral representation) | beat weight | linear weight | total weight |

    mod_12 argument (false by default), is for computing the CoE in module 12 representation
    having C=0 and B=11.

    lin_w argument is for the linear weight added to the height of the pitch (the lower the more important)
    the weighting is from 1-factor being lin_w the factor, default is lin_w=1 (no extra weighting)

    sbeat_w argument is the extra weight for the first beat in the sequence, e.g. the notes that fall in
    exactly the first beat. This works for any measure, the argument is an array of array: Array{Array{Float64,1},1}
    first array corresponds to the number of beat to be modified:
    If in a 4/4 measure you want to weight the 1st and 3rd notes the array would be [1.,3.] and second [w1,w2]
    being w1 and w2 the weights for each beat.
    Default value is no extra weight, sbeat_w = [[1.],[1.]].
"""
function get_key_sequence(notes_chunk; r=1, h=sqrt(2/15), mod_12=false,all_keys=all_keys, pos_all_keys=pos_all_keys, sbeat_w=[[1.],[1.]], lin_w=1)
    ptcs = notes_chunk[:,6] #Pitches
    durs = notes_chunk[:,5] #durations
    pbeat = notes_chunk[:,1] #beat where the note starts
    beat_w = ones(length(durs)) #array of the beat weights
    for b = 1:length(sbeat_w[1])
        loc_b = findall(x-> x==sbeat_w[1][b], pbeat) #finding all notes that start at beat sbeat_w[1][b]
        beat_w[loc_b] .= sbeat_w[2][b] #this is the weight.
    end
    notas, n_we = get_local_lin_w(ptcs, lin_w) #doing the linear weight in the pitches
    ii = vcat(map(x-> findall(y-> y==x, notas), ptcs)...)
    #println(ptcs)
    b_wei = n_we[ii] #getting the linear weight for every note i n the array of pitches
    if mod_12
        mod12_seq = get_cfpitch_mod12(ptcs)
        spi_ix = reorder_seq_closest(mod12_seq).+12
        #println("WARNING! \n Doing module 12 notes.")
    else
        spi_ix = get_cfpitch(ptcs)
        if length(spi_ix) > 1
            max_dif = 10
            while max_dif > 6
                nex_seq = reorder_seq_closest(spi_ix)
                max_dif = maximum(map(x-> abs(x),diff(nex_seq)))
                spi_ix = nex_seq
            end
        end
    end
    spi_p = map(x-> get_pitch(x, r=r, h=h), spi_ix) #getting the location (x,y,z) for each pitch
    t_ws = map((x,y,z)-> x*y*z, beat_w,durs, b_wei) #computing the total weights
    cv_i = map((x,y)-> x*y, t_ws, spi_p) / sum(t_ws) #computing the location of the pitches with their relative weights

    c_i = sum(cv_i) #finding the center of effect

    d_to_keys = map(x-> euclidean(c_i, x)^2, pos_all_keys) #computing the eclidean distance to all keys

    ranking = sortperm(d_to_keys) #ranking the distances from the closest to the farthest

    return c_i, [all_keys[ranking] d_to_keys[ranking]], [notes_chunk midi_note_names[ptcs.+1] spi_ix beat_w b_wei t_ws]
end

"""
    get_rand_key_sequence(notes_chunk; n_iter=1000,r=1, h=sqrt(2/15),
        mod_12=false,all_keys=all_keys, pos_all_keys=pos_all_keys, sbeat_w=[[1.],[1.]], lin_w=1)

    Returns a list of all the closest keys with their respective normalized frequency and the average distance:
        | Closest Key | Frequency | Average Distance |
    The algorithm for finding the closest key is the get_key_sequence which computes the center of effect of a set of notes
    since the center of effect is sensitive to the order in time of the notes, this function computes n_iter times the
    get_key_sequence function for a random shuffle version of the original set of notes, after that it returns the
    frequency distribution of the closest keys and the average distance to that key.
"""
function get_rand_key_sequence(notes_chunk; center_effect=false, n_iter=1000,r=1, h=sqrt(2/15), mod_12=false,all_keys=all_keys, pos_all_keys=pos_all_keys, sbeat_w=[[1.],[1.]], lin_w=1)
    n_pitches = size(notes_chunk,1)
    c_e = []
    if n_iter != 1000
    	n_iter = n_iter
    elseif n_pitches < 7
        n_iter = factorial(n_pitches)
    end
    key_list = Array{Any}(undef, n_iter,2)
    for t=1:n_iter
        ixes = [j for j=1:n_pitches]
        f_m = notes_chunk[shuffle(ixes),:]
        out = get_key_sequence(f_m, mod_12=mod_12,r=r,h=h,all_keys=all_keys,pos_all_keys=pos_all_keys,sbeat_w=sbeat_w,lin_w=lin_w)
        push!(c_e, round.(out[1], digits = 4))
        m_key = out[2][1,:]
        key_list[t,1] = m_key[1]; key_list[t,2] = m_key[2]
    end
    rf = get_rank_freq(key_list[:,1])
    n_keys = size(rf,1)
    t_count = sum(rf[:,2])
    k_vals = Array{Any}(undef, n_keys,3)
    for i = 1:n_keys
        m_v = mean(key_list[findall(x-> x==rf[i,1], key_list[:,1]),2])
        k_vals[i,1] = rf[i,1]; k_vals[i,2] = rf[i,2] / t_count; k_vals[i,3] = m_v
    end
    if center_effect
        return c_e, k_vals
    else
        return k_vals
    end
end

"""
    get_kseq_properties(key_seqs; p_distance=false)

    Returns properties computed from a sequence of keys: An array of tables | Key | Distance | if p_distance is true
    if p_distance is false the imput should have the probabilities instead of distances | Key | Probability |.
    The properties computed are:
        - Sequence of Keys (e.g. ["G","A","B"...])
        - Number of measure (an Int array with the number of measure)
        - Sequence of Key Uncertainty (Shannon entropy for the distribution of probabilities
          for each key in the key distances)
        - Average Uncertainty
        - Key entropy for the piece, computed for the distribution of keys in the piece.
        - Array of probability distributions of the keys;
          if is given distances the probabilites are computed as p(x) = exp(-x) and normalized.
"""
function get_kseq_properties(key_seqs; p_distance=false)
    if p_distance
        p_distros = get_distance_dists(key_seqs)
        dk_seq_ent = get_distance_seq_entropy(key_seqs)
    else
        p_distros, dk_seq_ent = get_kseq_pdent(key_seqs)
    end
    n_k_s = get_KS(key_seqs)

    chord_n = convert(Array{String,1},n_k_s)

    D = get_freq(n_k_s)
    p_keys = collect(values(D))./sum(collect(values(D)))
    ent_dseq = mapreduce(x-> -x*log2(x),+,p_keys)

    av_e = mean(filter(x-> !isnan(x),dk_seq_ent))
    return chord_n,  dk_seq_ent, av_e, ent_dseq
end
"""
    get_kseq_randomized(measures)

    Given a set of group of notes, it computes the most likely key for each set in a stochastic fashion.
    
    Returns a list of keys, list of uncertanty values for the key sequence, mean uncertainty and
    the Entropy for the key distribution in the piece (set of group of notes).

    For more details see function get_kseq_properties

"""
function get_kseq_randomized(measures)
    keys = []
    for i = 1:length(measures)
        push!(keys,get_rand_key_sequence(measures[i])[:,1:2])
    end
    return get_kseq_properties(keys)
end
##--
"""
    divide_measure(notes_chunk)

    Returns a measure divided by beat, it does not take the full lenght
    for each note if the note is longer than the beat.
    
"""
function divide_measure(notes_chunk)
    t_sig = notes_chunk[:,2]
    n_win = parse(Int64,split(t_sig[1],"/")[1])
    ptcs = notes_chunk[:,6] #Pitches
    durs = notes_chunk[:,5] #durations
    pbeat = notes_chunk[:,1] #beat where the note starts
    ini = notes_chunk[:,3]
    fin = notes_chunk[:,4]
    #we need to divide the measure into the number of notes (in the denominator)
    tam_win = (maximum(fin) - minimum(ini))/n_win #size of the window (per beat)
    sub_chunks = []

    r_ini = ini.-ini[1]
    r_fin = fin.-ini[1]

    for nw in 0:(n_win-1)
        sub_t = []
        ini_win = nw * tam_win 
        fin_win = nw * tam_win + tam_win
        ix_n = unique(vcat(findall(x-> x < fin_win && x>= ini_win, r_ini),findall(x-> x<= fin_win && x>ini_win, r_fin)))
        #println(ix_n," window $(nw+1)")
        if isempty(ix_n)
            ix = intersect(findall(x-> x>fin_win, r_fin),findall(x-> x<ini_win, r_ini))
            if isempty(ix) continue; end
            #println("$ix window $(nw+1)")
            for p in ix
                push!(sub_t, [pbeat[p] t_sig[p] ini[p] fin[p] tam_win ptcs[p]])
            end
        else
            for p in ix_n
                #println("window number $(nw+1)")
                #println(r_ini[p],"\t",r_fin[p], "\t", ini_win,"\t", fin_win)
                if r_fin[p] <= fin_win
                    if r_ini[p]< ini_win #starts before the beat
                        push!(sub_t, [pbeat[p] t_sig[p] ini[p] fin[p] r_fin[p] - ini_win ptcs[p]]) #case it finishes before or in the beat
                    else
                        push!(sub_t, [pbeat[p] t_sig[p] ini[p] fin[p] r_fin[p] - r_ini[p] ptcs[p]])
                    end
                else
                    if r_ini[p] < ini_win
                        push!(sub_t, [pbeat[p] t_sig[p] ini[p] fin[p] fin_win - ini_win ptcs[p]]) #case it is longer than the beat
                    else
                        push!(sub_t, [pbeat[p] t_sig[p] ini[p] fin[p] fin_win - r_ini[p] ptcs[p]])
                    end
                end
            end
        end
        push!(sub_chunks, vcat(sub_t...))
    end
    return sub_chunks
end


function get_distance_to_keys(c_i)
    d_to_keys = round.(map(x-> euclidean(c_i, x), pos_all_keys), digits = 4) #computing the eclidean distance to all keys

    ranking = sortperm(d_to_keys) #ranking the distances from the closest to the farthest

    return [all_keys[ranking] d_to_keys[ranking] pos_all_keys[ranking]]
end


function cluster_notes(ptcs)
    p_m12 = map(x-> mod(x,12),ptcs)
    spi_notes = get_cfpitch(p_m12)
    low_notes = findall(x-> x<=6, spi_notes)
    high_notes = findall(x-> x>6, spi_notes)
    if !isempty(high_notes) && !isempty(low_notes)
        if length(high_notes) <= length(low_notes)
            for i in 1:length(high_notes)
                spi_notes = shift_outlier(spi_notes, high_notes[i])
            end
        else
            for i in 1:length(low_notes)
                spi_notes = shift_outlier(spi_notes, low_notes[i])
            end
        end
    end
    dmean = Float64[]
    d_olier = Float64[]
    oliers = Int64[]
    new_spi = []
    conv = false
    while conv == false
        dt_mean = zeros(length(spi_notes))
        for i in 1:length(spi_notes)
            dt_mean[i] = abs(spi_notes[i] - mean(spi_notes[1:end .!= i]))
        end
        d, olier = findmax(dt_mean)
        push!(dmean, mean(dt_mean)); push!(oliers, olier);push!(d_olier, round(d, digits=4)); push!(new_spi, spi_notes)
        #println(mean(dt_mean),'\t', olier, '\t', p_cf)
        #println(mean(dt_mean),'\t', olier)
        if length(oliers) > 20 && length(unique(oliers[end-4:end])) <= 2
            conv = true
            break
        end
        spi_notes = shift_outlier(spi_notes, olier)
    end
    mmin = findmin(dmean[end-3:end])[2]
    dmin = findmin(d_olier[end-3:end])[2]
    if dmean[end-3:end][mmin] == d_olier[end-3:end][dmin]
        spi_out = new_spi[end-3:end][mmin]
    elseif oliers[end-3:end][mmin] == oliers[end-3:end][dmin]
        spi_out = new_spi[end-3:end][dmin]
    else
        spi_out = new_spi[end-3:end][mmin]
    end
    return spi_out
end

function shift_outlier(notes, olier)
    dif = notes[olier] - median(notes)
    notes_new = copy(notes)
    if dif > 0
        notes_new[olier] = notes_new[olier] - 12
    elseif dif < 0
        notes_new[olier] = notes_new[olier] + 12
    end
    return notes_new
end

function get_center_effect(chunk_notes; r=1, h=sqrt(2/15), mod_12=false,all_keys=all_keys, pos_all_keys=pos_all_keys, sbeat_w=[[1.],[1.]], lin_w=1)    
    ptcs = chunk_notes[:,6]
    durs = chunk_notes[:,5]
    pbeat = chunk_notes[:,1]
    beat_w = ones(length(durs)) #array of the beat weights
    for b = 1:length(sbeat_w[1])
        loc_b = findall(x-> x==sbeat_w[1][b], pbeat) #finding all notes that start at beat sbeat_w[1][b]
        beat_w[loc_b] .= sbeat_w[2][b] #this is the weight.
    end
    notas, n_we = get_local_lin_w(ptcs, lin_w) #doing the linear weight in the pitches
    ii = vcat(map(x-> findall(y-> y==x, notas), ptcs)...)
    #println(ptcs)
    b_wei = n_we[ii] #getting the linear weight for every note i n the array of pitches
    
    spi_ix = cluster_notes(ptcs) .+ 24

    spi_p = map(x-> get_pitch(x, r=r, h=h), spi_ix) #getting the location (x,y,z) for each pitch
    t_ws = map((x,y,z)-> x*y*z, beat_w,durs, b_wei) #computing the total weights
    cv_i = map((x,y)-> x*y, t_ws, spi_p) / sum(t_ws) #computing the location of the pitches with their relative weights
    c_i = sum(cv_i) #finding the center of effect
    return c_i
end

function get_distance_ces(ce1, ce2)
    z_dif = ce2[3] - ce1[3]
    while abs(z_dif) > h_octav / 2 #translating over z to be in the same octave (same CE region)
        if z_dif > 0
            ce2[3] = ce2[3] - h_octav
        else
            ce2[3] = ce2[3] + h_octav
        end
        z_dif = ce2[3] - ce1[3]
    end
    return round(euclidean(ce1,ce2), digits = 4)
end

function get_xml_df(piece_xml)
    piece = get_piece_by_measure(piece_xml, csv=false)
    df_piece = DataFrame(
        :Measure => convert(Array{Int64,1},piece[:,1]),
        :TimeSignature => piece[:,2],
        :StartQuarter => piece[:,3],
        :EndQuarter => piece[:,4],
        :Duration => piece[:,5],
        :Pitch => convert(Array{Int64,1},piece[:,6])
    )
    return df_piece
end
function get_csv_df(piece_csv)
    piece = get_piece_by_measure(piece_csv, csv=true)[1]
    num_mea = []
    for nm in 1:length(piece)
        push!(num_mea,[nm for i in 1:size(piece[nm],1)])
    end
    num_mea = vcat(num_mea...)
    piece = vcat(piece...)
    df_piece = DataFrame(
        :Measure => convert(Array{Int64,1},num_mea),
        :TimeSignature => piece[:,2],
        :StartTime => piece[:,3],
        :EndTime => piece[:,4],
        :Duration => piece[:,5],
        :Pitch => convert(Array{Int64,1},piece[:,6])
    )
    return df_piece
end

################################################################################
#####TOOLS


"""
    get_kseq_pdent(key_seqs)

    Returns the probabilities and the entropies of a key sequence (Array of tables):
    | Key | Probability | using Shannon definition in module 2.
"""
function get_kseq_pdent(key_seqs)
    pd = []; entro = []
    for i = 1:length(key_seqs)
        push!(pd, key_seqs[i][:,2])
        push!(entro, mapreduce(x-> -x*log2(x),+,key_seqs[i][:,2]))
    end
    return pd, entro
end
"""
    reorder_seq_closest(ex_seq)

    Returns a re-ordered sequence of pitches in the spiral representation order,
    the sequence satisfies the property that any step between two consecutive notes
    is not farther than 6, meaning that the steps are to the closest note not to the
    real one.
"""
function reorder_seq_closest(ex_seq)
    nex_seq = []
    push!(nex_seq,ex_seq[1])
    for i = 1:(length(ex_seq)-1)
        delta = ex_seq[i+1] - nex_seq[i]
        if abs(delta) > 6
            if delta < 0
                push!(nex_seq, nex_seq[i]+(12+delta))
            else
                push!(nex_seq, nex_seq[i]-(12-delta))
            end
        else
            push!(nex_seq, ex_seq[i+1])
        end
    end
    return nex_seq
end

"""
    get_local_lin_w(ptcs, factor)

    Returns a list of pitches and their respective weights from a sequence of pitches (ptcs)
    and a weighting factor (factor), the weighting gives more importance to lower notes
    within the 1-factor range.
"""
function get_local_lin_w(ptcs, factor)
    notas = sort(unique(ptcs))
    if length(notas) == 1
        n_we = []
        push!(n_we, 1)
    else
        p_min = notas[1]
        p_max = notas[end]
        tam = length(notas)
        n_we = zeros(tam)
        for i = 1:tam
            n_we[i] = 1 + (factor-1) * (notas[i] - p_max)/(p_min-p_max)
        end
    end
    return notas, n_we
end

"""
    get_KS(prob_keys)

    Returns a sequence of the most likely key, given a sequence of 2-dimensional arrays
    with the format |key | distance|
"""
function get_KS(prob_keys)
    k_s = []
    for i = 1:length(prob_keys)
        if prob_keys[i] == "No Key/Notes"
            #println(i)
            #push!(out_triads, [0,0,0])
        else
            push!(k_s, prob_keys[i][1,1])
        end
    end
    return k_s
end

"""
    get_icontent_key(dis)

    Returns the information content (IC) of the most likely key (most probable)
    from an array of distances "dis".
"""
function get_icontent_key(dis)
    prob_d = map(x-> exp(-x), dis)
    prob_d = prob_d ./ sum(prob_d)
    if sum(prob_d) == 0 || isnan(sum(prob_d))
        prob_d = [1/length(dis) for i = 1:length(dis)]
    end
    return log2(1/prob_d[1]) #returns the information content for the most likely key
end

"""
    get_distance_entropy(dis)

    Returns the Shannon Entropy log2 from an array of distances "dis"
    The probabilites p(x) are computed in the form of p(x) = exp(-x),
    where "x" is the distance.
"""
function get_distance_entropy(dis)
    prob_d = map(x-> exp(-x), dis)
    prob_d = prob_d ./ sum(prob_d)
    if sum(prob_d) == 0 || isnan(sum(prob_d))
        prob_d = [1/length(dis) for i = 1:length(dis)]
    end
    return mapreduce(x-> -x*log2(x),+,prob_d)
end
"""
    get_distros(dis)

    Returns a probability distribution (a normalized array) for a set of distances "dis".
    The probabilites p(x) are computed in the form of p(x) = exp(-x),
    where "x" is the distance.
"""
function get_distros(dis)
    prob_d = map(x-> exp(-x), dis)
    prob_d = prob_d ./ sum(prob_d)
    if sum(prob_d) == 0 || isnan(sum(prob_d))
        prob_d = [1/length(dis) for i = 1:length(dis)]
    end
    return prob_d
end

"""
    get_average_entropy(key_seq)

    Returns the average Uncertainty (Shannon Entropy) of a key sequence
    in the format |key | distance|.
    The probabilites p(x) are computed in the form of p(x) = exp(-x),
    where "x" is the distance.
"""
function get_average_entropy(key_seq)
    dis = []
    for i = 1:length(key_seq)
        if key_seq[i] == "No Key/Notes"
            continue
        end
        #println(i)
        push!(dis, convert(Array{Float64,1},key_seq[i][2][:,2]))
    end
    dist = map(x-> get_distance_entropy(x), dis)
    filter!(x-> !isnan(x),dist)
    return sum(dist) / length(dist)
end

"""
    get_distance_seq_entropy(key_seq)

    Returns the Uncertainties (Shannon Entropy) of a key sequence
    in the format |key | distance|, uncertainty is computed for each key.
    The probabilites p(x) are computed in the form of p(x) = exp(-x),
    where "x" is the distance.
"""
function get_distance_seq_entropy(key_seq)
    dis = []
    for i = 1:length(key_seq)
        if key_seq[i] == "No Key/Notes"
            continue
        end
        #println(i)
        push!(dis, convert(Array{Float64,1},key_seq[i][:,2]))
    end
    dist = map(x-> get_distance_entropy(x), dis)
    #filter!(x-> !isnan(x),dist)
    return dist
end

"""
    get_distance_dists(key_seq)

    Returns an array of Array{Float64,1} of the distances for all the keys to be able
    to compute probabilites and information measures with them. Input is a sequence of
    keys with their respective distance (array of 2-dimensional arrays |key | distance |)
"""
function get_distance_dists(key_seq)
    dis = []
    for i = 1:length(key_seq)
        if key_seq[i] == "No Key/Notes"
            continue
        end
        #println(i)
        push!(dis, get_distros(convert(Array{Float64,1},key_seq[i][:,2])))
    end
    return dis
end

"""
    get_freq(s)

    Returns a dictionary with the "key" and the "frequency" from a list of elements "s"
    frequency is not normalized, returns how many times each symbol appears in the list.
"""
function get_freq(s)
    T = Dict{Any,Int64}()
    for i = 1:length(s)
        T[s[i]] = get(T, s[i],0) + 1
    end
    return T
end

"""
    get_rank_freq(series)

    Returns a 2-Dimensional Array of symbols and their respective frequencies
    from a series of data. The output is ordered by rank (most to least frequent).
"""
function get_rank_freq(series)
    tam = length(series)
    M = Dict{Any,Int64}()
    for i = 1:tam
        M[series[i]] = get(M, series[i], 0) + 1
    end
    dist = sort(collect(M), by = tuple -> last(tuple), rev=true)
    rf = Array{Any}(undef, length(dist),2)
    for i = 1:length(dist)
        rf[i,1] = dist[i][1]; rf[i,2] = dist[i][2]
    end
    return rf
end

"""
    get_hamming_distance(s1,s2)

    Returns the Hamming distance between two sequences (Arrays of any type),
    only works for arrays of the same lenght.
"""
function get_hamming_distance(s1,s2)
    if length(s1) != length(s2)
        println("STRINGS ARE NOT SAME LENGHT!!!!")
        return NaN
    else
        h_d = sum([s1[i] != s2[i] for i = 1:length(s1)])
    end
    return h_d
end
################################################################################

################################################################################
#OLD CODE
function duration_series(voice)
    st = 1
    n_p = []
    t_c = []
    n_c = 0
    i = 1
    while i < length(voice)

        if voice[st]==voice[i]
         n_c+=1
         i+=1
        else
         push!(n_p, voice[st])
         push!(t_c, n_c)
         st += n_c
         n_c = 0
         i = st
        end

    end
    return t_c, n_p
end

function get_all_notes(s)
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
    notes_s = []
    for i = 1:nv
        push!(notes_s,n_vs[i][:,3] )
    end
    a_n = convert(Array{Int64,1},sort(unique(vcat(notes_s...))))
    return a_n
end

function get_lin_w(piece, factor)
    notas = get_all_notes(piece)
    p_min = notas[1]
    p_max = notas[end]
    tam = length(notas)
    n_we = zeros(tam)
    for i = 1:tam
        n_we[i] = 1 + (factor-1) * (notas[i] - p_max)/(p_min-p_max)
    end
    return notas, n_we
end


function find_key_seq(durs, ptcs; r=1, h=sqrt(2/15), all_keys=all_keys, pos_all_keys=pos_all_keys)

    spi_ix = get_cfpitch(ptcs)
    spi_p = map(x-> get_pitch(x, r, h), spi_ix)

    cv_i = map((x,y)-> x*y, durs, spi_p) / sum(durs)

    c_i = sum(cv_i)

    d_to_keys = map(x-> euclidean(c_i, x)^2, pos_all_keys)

    ranking = sortperm(d_to_keys)
    return c_i, [all_keys[ranking] d_to_keys[ranking]]
end

function find_key_seq_bass(durs, ptcs, r, h, all_keys, pos_all_keys, notas, n_we)
    ii = vcat(map(x-> findall(y-> y==x, notas), ptcs)...)
    #println(ptcs)
    b_wei = n_we[ii]
    spi_ix = get_cfpitch(ptcs)
    spi_p = map(x-> get_pitch(x, r, h), spi_ix)
    #println(length(spi_p))
    #println(b_wei,'\n',durs)
    t_ws = map((x,y)-> x*y, durs, b_wei) #total weights
    #println(length(t_ws))
    cv_i = map((x,y)-> x*y, t_ws, spi_p) / sum(t_ws)

    c_i = sum(cv_i)

    d_to_keys = map(x-> euclidean(c_i, x)^2, pos_all_keys)

    ranking = sortperm(d_to_keys)
    return c_i, [all_keys[ranking] d_to_keys[ranking]]
end


function plot_path(pitches, minor_keys, major_keys, c_i_walk)
    ppits =hcat(pitches...)
    pits_x = ppits[1,:]
    pits_y = ppits[2,:]
    pits_z = ppits[3,:]

    ma_c = hcat(major_chords...)
    ma_c_x = ma_c[1,:]
    ma_c_y = ma_c[2,:]
    ma_c_z = ma_c[3,:]

    mi_c = hcat(minor_chords...)
    mi_c_x = mi_c[1,:]
    mi_c_y = mi_c[2,:]
    mi_c_z = mi_c[3,:]

    ma_k = hcat(major_keys...)
    ma_k_x = ma_k[1,:]
    ma_k_y = ma_k[2,:]
    ma_k_z = ma_k[3,:]

    mi_k = hcat(minor_keys...)
    mi_k_x = mi_k[1,:]
    mi_k_y = mi_k[2,:]
    mi_k_z = mi_k[3,:]

    ci_w = hcat(c_i_walk...)
    ci_w_x = ci_w[1,:]
    ci_w_y = ci_w[2,:]
    ci_w_z = ci_w[3,:]
    pt_path = []
    push!(pt_path,plot(pits_x, pits_y, pits_z, m=:o, label="Pitches"))
    #plot!(ma_c_x, ma_c_y, ma_c_z, m=:o, label="Major Chords")
    #plot!(mi_c_x, mi_c_y, mi_c_z, m=:o, label="Minor Chords")
    push!(pt_path,plot!(ma_k_x, ma_k_y, ma_k_z, m=:star,ms=3, label="Major Keys"))
    push!(pt_path,plot!(mi_k_x, mi_k_y, mi_k_z, m=:star,ms=3, label="Minor Keys"))
    push!(pt_path,plot!(ci_w_x, ci_w_y, ci_w_z, m=:o, ms=4, label="Walk", color=:red))
    return pt_path
end

function get_key_pvoices(ch_v, r, h)
    nv = size(ch_v)[1]
    n_chunks = size(ch_v[1], 1)
    voice_keys = []
    for v = 1:nv
        durs = []
        pits = []
        for i = 1:n_chunks
            if isempty(ch_v[v][i])
                push!(durs, 0)
                push!(pits, 0)
            else
                push!(durs,ch_v[v][i][:,2]-ch_v[v][i][:,1])
                push!(pits,ch_v[v][i][:,3])
            end
        end

        chunk_keys = []

        for i = 1:n_chunks
            if durs[i] == 0
                push!(chunk_keys, "No Key/Notes")
            else
                push!(chunk_keys, find_key_seq(durs[i], pits[i],r,h, all_keys, pos_all_keys)[2])
            end
        end
        push!(voice_keys, chunk_keys)
    end

    return voice_keys
end

function get_key_allvoices(ch_v, r, h)
    nv = size(ch_v)[1]
    n_chunks = size(ch_v[1], 1)
    chunk_keys = []
    for i = 1:n_chunks
        durs = []
        pits = []
        for v = 1:nv
            if isempty(ch_v[v][i])
                push!(durs, 0)
                push!(pits, 0)
            else
                push!(durs,ch_v[v][i][:,2]-ch_v[v][i][:,1])
                push!(pits,ch_v[v][i][:,3])
            end
        end
        all_durs = vcat(durs...)
        all_pits = vcat(pits...)
        if sum(all_durs) == 0
            push!(chunk_keys, "No Key/Notes")
        else
            push!(chunk_keys, find_key_seq(all_durs, all_pits,r,h, all_keys, pos_all_keys))
        end
    end

    return chunk_keys
end
function get_key_allvoices(ch_v, r, h, notas, n_we)
    nv = size(ch_v)[1]
    n_chunks = size(ch_v[1], 1)
    chunk_keys = []
    for i = 1:n_chunks
        durs = []
        pits = []
        for v = 1:nv
            if isempty(ch_v[v][i])
                push!(durs, 0)
                push!(pits, 0)
            else
                push!(durs,ch_v[v][i][:,2]-ch_v[v][i][:,1])
                push!(pits,ch_v[v][i][:,3])
            end
        end
        rpits = findall(x->x!=0, pits)

        if isempty(rpits)
            push!(chunk_keys, "No Key/Notes")
        else
            all_durs = vcat(durs[rpits]...)
            all_pits = vcat(pits[rpits]...)
            push!(chunk_keys, find_key_seq_bass(all_durs, all_pits,r,h, all_keys, pos_all_keys, notas, n_we))
        end
    end

    return chunk_keys
end
function get_key_allvoices(ch_v, r, h, f_w)
    nv = size(ch_v)[1]
    n_chunks = size(ch_v[1], 1)
    chunk_keys = []
    for i = 1:n_chunks
        durs = []
        pits = []
        for v = 1:nv
            if isempty(ch_v[v][i])
                push!(durs, 0)
                push!(pits, 0)
            else
                push!(durs,ch_v[v][i][:,2]-ch_v[v][i][:,1])
                push!(pits,ch_v[v][i][:,3])
            end
        end
        rpits = findall(x->x!=0, pits)

        if isempty(rpits)
            push!(chunk_keys, "No Key/Notes")
        else
            all_durs = vcat(durs[rpits]...)
            all_pits = vcat(pits[rpits]...)
            notas, n_we = get_local_lin_w(all_pits, f_w)
            push!(chunk_keys, find_key_seq_bass(all_durs, all_pits,r,h, all_keys, pos_all_keys, notas, n_we))
        end
    end

    return chunk_keys
end

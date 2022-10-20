"""
Script layout
generates a synthetic promotoe based on promoter database 
"""

cd("/home/bensiv/Dropbox/202207-RemilkProject/Localized_output/")

# activate project enviroment
using Pkg
Pkg.activate("/home/bensiv/Dropbox/202207-RemilkProject/")

# get command line arguments from the user
using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--output_dir", "-o"
            help = "outpput directory"
            default = pwd()
        "--promoter", "-p"
            help = "input promoter file as fasta"
            required = true
        "--regulatory_elements", "-r"
            help = "input regulatory_elements file as tsv"
            required = true
        "--interactive", "-i"
            help = "save the commandline arguments as a yaml file for interactive running"
            action = :store_true
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        println("  $arg  =>  $val")
    end

    return parsed_args
end

parsed_args = main()

using YAML
include("stop.jl")

if parsed_args["interactive"]
    YAML.write_file("parsed_args.yml", parsed_args)
    stop("stopped and saved the arguments as parsed_args.yml")
    # run `parsed_args = YAML.load_file("parsed_args.yml")` to retrive it in the script
end

# precompile packages
using CSV, DataFrames, DataFramesMeta
using BioSequences, FASTX


# load promoters
promoters = open(FASTA.Reader, parsed_args["promoter"])

# access the promotor and store them as a dictionary
promotersDict = Dict()
for record in promoters
    push!(promotersDict, "$(FASTA.identifier(record))" => FASTA.sequence(record))
end

# load regulatory elements as a dataframe
regulatory_elements = CSV.read(parsed_args["regulatory_elements"], DataFrame)
# unique_RE = @linq regulatory_elements |> unique(:Sequence) |> transform(:Sequence = LongDNA{4}.(:Sequence))
regulatory_elements = @linq regulatory_elements |> transform(:Sequence = LongDNA{4}.(:Sequence))

# # find all the occurences of regulatory elements in the promoters
# FeaturesTable = DataFrame(Promoter = [], Motif = [], Occurences = [], Direction = [])
# for key in keys(promotersDict)
#     for motif in eachrow(unique_RE)
#         sense = findall(motif.Sequence, string(promotersDict[key]))
#         anti_sense = findall(string(reverse_complement(LongDNA{4}(motif.Sequence))), string(promotersDict[key]))
#         occurences = [sense;anti_sense]
#         directions = [fill("+", length(sense)); fill("-", length(anti_sense))]
#         n = length(occurences)
#         if !isempty(occurences)
#             FeaturesTable = [FeaturesTable; DataFrame(Promoter = fill(key, n), Motif = fill(motif.Name, n), Occurences = occurences, Direction = directions)]
#             # push!(Features, ["$key", "$(motif.Name)", occurences, Directions])
#         end
#     end
# end


# find all the occurences of regulatory elements in the promoters
function find_features(promoters_dict, reg_df)
    FeaturesTable = DataFrame(Promoter = [], Motif = [], Occurences = [], Direction = [])

    for key in keys(promoters_dict)
        for motif in eachrow(reg_df)
            
            query = ExactSearchQuery(motif.Sequence, iscompatible)
            # rev_query = ExactSearchQuery(reverse_complement(motif.Sequence), iscompatible)
    
            sense = findall(query, promoters_dict[key])
            # anti_sense = findall(rev_query, promoters_dict[key])
            
            # occurences = [sense;anti_sense]
            occurences = sense
            # directions = [fill("+", length(sense)); fill("-", length(anti_sense))]
            directions = fill("+", length(sense))
            n = length(occurences)
            if !isempty(occurences)
                FeaturesTable = [FeaturesTable; DataFrame(Promoter = fill(key, n), Motif = fill(motif.Name, n), Occurences = occurences, Direction = directions)]
                # push!(Features, ["$key", "$(motif.Name)", occurences, Directions])
            end
        end
    end
    
    return FeaturesTable
end

FeaturesTable = find_features(promotersDict, unique_RE)


# using Pipe

# function motif_locations(motif)
#     locations = []
#     for key in keys(Features)
#         if motif in keys(Features[key])
#             push!(locations, Features[key][motif])
#         end
#     end
    
#     # return @pipe locations |> vcat(_...) |> first.(_) |> vcat(_...)        
#     return @pipe locations |> vcat(_...) |> collect.(_) |> vcat(_...)        
# end

function motif_locations(features_table, motif, direction)
    temp = @linq features_table |> where(:Motif .== motif .&& :Direction .== direction)
        
    return first.(temp.Occurences)
end


RYrepeat = motif_locations(FeaturesTable, "RY-repeat", "+")
TATAbox = motif_locations(FeaturesTable, "TATA-box1", "+")

using Plots, StatsPlots

density(TATAbox)
density(RYrepeat)

# all_features_found = unique(vcat([collect(keys(Features[x])) for x in keys(Features)]...))
# all_features_found = unique(FeaturesTable.Motif)
FeaturesTable_plus = @linq FeaturesTable |> where(:Direction .== "+")
all_features_found = unique(FeaturesTable_plus.Motif)

all_locations = []
for feature in all_features_found
    push!(all_locations, motif_locations(FeaturesTable, feature, "+"))
end

location_lengths = length.(all_locations)

using KernelDensity, Distributions

function FindPeak(data)
    den = kde(data)
    mval = findmax(den.density)
    return Int64(round(den.x[mval[2]])), den.density[mval[2]]
end

# FindPeak(RYrepeat)

peaks, den = [], []
for l in all_locations
    x = FindPeak(l)
    push!(peaks, x[1])
    push!(den, x[2])
end

features_optimal = DataFrame(feature = all_features_found, occure = location_lengths, peak = peaks, height = collect(den))
sort!(features_optimal, [:height,:occure], rev = true)

function FindAnotherPeak(data, from, to)
    subdata = data[data .> from .&& data .< to]
    den = kde(data)
    subden = kde(subdata)

    mval = findmax(subden.density)

    return Int64(round(subden.x[mval[2]])), den.density[mval[2]]
end

function FindMultiPeaks(motif, froms, tos, onlyplot)
    motif_locs = motif_locations(FeaturesTable, motif, "+")
    
    first_peak = FindPeak(motif_locs)
    motif_peaks = [first_peak[1]]
    motif_height = [first_peak[2]]

    for (f,t) in zip(froms, tos)
        p,h = FindAnotherPeak(motif_locs, f, t)
        push!(motif_peaks, p)
        push!(motif_height, h)
    end

    if onlyplot
        return density(motif_locs)
    else
        return hcat(motif_peaks, motif_height)
    end

end

TATAbox_den = kde(TATAbox)
RYrepeat_den = kde(RYrepeat)

density([TATAbox, RYrepeat], label = ["TATA-box" "RYrepeat"], linewidth = 3, grid = false, ylabel = "motif density", xlabel = "location (bp)", title = "motif density across multiple promoters")

plot([TATAbox_den.density RYrepeat_den.density], label = ["TATA-box" "RYrepeat"], linewidth = 3, grid = false, ylabel = "motif density", xlabel = "location (bp)", title = "motif density across multiple promoters")

plot(TATAbox_den.density, label = "TATA-box", linewidth = 3, grid = false, ylabel = "motif density", xlabel = "location (bp)", title = "motif density across multiple promoters")
plot!(RYrepeat_den.density, label = "RYrepeat", linewidth = 3, grid = false, ylabel = "motif density", xlabel = "location (bp)", title = "motif density across multiple promoters")

plot(TATAbox_den.density)
plot(RYrepeat_den.density)

savefig("motif_density.png")

# FindMultiPeaks("TATA-box1", [0,600], [200,1000], false)
# FindMultiPeaks("Light_res1", [155], [170], false)
# FindMultiPeaks("CAAT4", [1500], [3000], false)
# FindMultiPeaks("TATA-box2", [1800], [3000], false)
# FindMultiPeaks("CAAT5", [800], [1200], false)
# FindMultiPeaks("RAV1-A binding site motif", [0], [150], false)
# FindMultiPeaks("GATA promoter motif [LRE]", [1500], [3000], false)
FindMultiPeaks("RY-repeat promoter motif", [200], [600], false)
# FindMultiPeaks("MeJA_res1", [155], [200], false)
# FindMultiPeaks("CAAT3", [190],[300], false)
# FindMultiPeaks("W-box promoter motif", [0,800],[300,1200], false)
# FindMultiPeaks("ABA_res1", [165],[200], false)
# FindMultiPeaks("LFY consensus binding site motif", [700],[1000], false)
# FindMultiPeaks("MYB4 binding site motif", [0],[250], false)
# FindMultiPeaks("DPBF1&2 binding site motif", [1500],[2000], false)
FindMultiPeaks(features_optimal.feature[50], [],[], true)
findall(features_optimal.occure .> 1)

peak_to_add = copy(features_optimal)[50,:]
peak_to_add.peak = 187
peak_to_add.height = 0.0110962
push!(features_optimal, peak_to_add)

"""
How to find all the maximums automatic

find the first max.
search the min upstream from it.
search the next max upstream from it.
and repeat until the min equels max.
repeat downstream.
"""

mutable struct local_max
    location::Int64
    magnitude::Float64
end


function FindAllPeaks(features_table, motif, direction)
    # find all occurences
    motif_locs = motif_locations(features_table, motif, direction)

    den = kde(motif_locs)
    mval = findmax(den.density)
    max_peak = local_max(Int64(round(den.x[mval[2]])), den.density[mval[2]])
    
    all_peaks = local_max[max_peak]
    # bounds = []

    corrent_index = mval[2]

    # search upstream
    while all(den.density[corrent_index] .> den.density[corrent_index+1:end])
        corrent_index += 1
        if corrent_index == length(den.density)+1
            break
        end
    end

    if corrent_index < length(den.density)
        new_mval = findmax(den.density[corrent_index+1:end])
        # push!(bounds, (corrent_index+1:corrent_index+new_mval[2]-1))
        corrent_index += new_mval[2]
        next_peak = local_max(Int64(round(den.x[corrent_index])), den.density[corrent_index])
        push!(all_peaks, next_peak)
    end
    corrent_index = mval[2]

    # search downtream
    while all(den.density[corrent_index] .> den.density[1:corrent_index-1])
        corrent_index -= 1
        if corrent_index == 0
            break
        end
    end

    if corrent_index > 0
        new_mval = findmax(den.density[1:corrent_index-1])
        # push!(bounds, (new_mval[2]+1:corrent_index-1))
        corrent_index = new_mval[2]
        next_peak = local_max(Int64(round(den.x[corrent_index])), den.density[corrent_index])
        push!(all_peaks, next_peak)
    end
    
    # order the peaks list by magnitude
    decs_order = sortperm([p.magnitude for p in all_peaks], rev = true)
    all_peaks = all_peaks[decs_order]

    return all_peaks
end


features_optimal = DataFrame(feature = all_features_found, occure = location_lengths)

feature_peaks = DataFrame(feature = String[], location = Int64[], magnitude = Float64[])
for feature in FeaturesTable.Motif
    peaks = FindAllPeaks(FeaturesTable, feature, "+")
    feature_peaks = [feature_peaks; DataFrame(feature = fill(feature,length(peaks)), location = [p.location for p in  peaks], magnitude = [p.magnitude for p in  peaks])]
end

features_optimal = leftjoin(features_optimal, feature_peaks, on = :feature)

features_optimal = @pipe features_optimal |> groupby(_, [:feature, :location]) |> combine(_, :magnitude => maximum => :magnitude, :occure => maximum => :occure) |> sort(_, [:magnitude,:occure], rev = true)

features_optimal.peak2 = [1835; fill(missing, nrow(features_optimal)-1)]
# features_optimal.peak3 = [93; fill(missing, nrow(features_optimal)-1)]

features_optimal.peak2[21] = 382
# features_optimal.peak3[8] = 2352

# CSV.write("features_optimal.csv", features_optimal)
CSV.write("Psfeatures_optimal.csv", features_optimal)

features_optimal[1:12,:]
# ==========================================================================================================

# assemble promoter

features_optimal = CSV.read("Gmfeatures_optimal.csv", DataFrame)

rename!(features_optimal, :length => :coverage)

template_promoter = LongDNA{4}(fill(DNA_N, 1001))

assembly = template_promoter

unique_RE.length = length.(unique_RE.Sequence)
rename!(unique_RE, :Name => :feature)

features_optimal = leftjoin(features_optimal, unique_RE, on = :feature)

sort!(features_optimal, :coverage, rev = true)

features_optimal[2,"peak2"] = 858

from = features_optimal[2,"peak1"]
to = from + length(unique_RE[unique_RE.Name .== features_optimal[2,"feature"], "Sequence"][1]) - 1
assembly[from:to] = LongDNA{4}(unique_RE[unique_RE.Name .== features_optimal[1,"feature"], "Sequence"][1])


for row in eachrow(features_optimal)
    from = row.peak1
    to = from + row.length - 1
    assembly[from:to] = LongDNA{4}(row.Sequence)
end

for row in eachrow(features_optimal)
    from = row.peak1
    to = from + row.length - 1
    if assembly[from:to] != LongDNA{4}(row.Sequence)
        println(row.feature)
    end
end
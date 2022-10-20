"""
Script layout
generates a synthetic promotoe based on promoter database 
"""

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
            default = "PLACE_plus_audited.tsv"
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
# regulatory_elements = @linq regulatory_elements |> unique(:Sequence) |> transform(:Sequence = LongDNA{4}.(:Sequence))
regulatory_elements = @linq regulatory_elements |> transform(:Sequence = LongDNA{4}.(:Sequence))


# find all the occurences of regulatory elements in the promoters
function find_features(promoters_dict, reg_df)
    FeaturesTable = DataFrame(Promoter = [], Motif = [], Length = [], Occurences = [], Direction = [])

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
                FeaturesTable = [FeaturesTable; DataFrame(Promoter = fill(key, n), Motif = fill(motif.Name, n), Length = fill(motif.Length, n), Occurences = occurences, Direction = directions)]
                # push!(Features, ["$key", "$(motif.Name)", occurences, Directions])
            end
        end
    end
    
    return FeaturesTable
end

FeaturesTable = find_features(promotersDict, regulatory_elements)

function motif_locations(features_table, motif, direction)
    temp = @linq features_table |> where(:Motif .== motif .&& :Direction .== direction)
        
    return first.(temp.Occurences)
end


FeaturesTable_plus = @linq FeaturesTable |> where(:Direction .== "+")
all_features_found = unique(FeaturesTable_plus.Motif)

all_locations = []
for feature in all_features_found
    push!(all_locations, motif_locations(FeaturesTable, feature, "+"))
end

location_lengths = length.(all_locations)

using KernelDensity, Distributions

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
    # motif_locs = motif_locations(FeaturesTable, "CAATBOX1", "+")

    den = kde(motif_locs)
    mval = findmax(den.density)
    max_peak = local_max(Int64(round(den.x[mval[2]])), den.density[mval[2]])
    
    # all_peaks = local_max[max_peak]
    all_peaks_df = DataFrame(
        mval = mval[2], 
        location = max_peak.location, 
        magnitude = max_peak.magnitude, 
        upper_bound = Union{Missing, UnitRange{Int64}}[missing], 
        lower_bound = Union{Missing, UnitRange{Int64}}[missing], 
        searched = false
    )

    # corrent_index = mval[2]

    while !all(all_peaks_df.searched)
        temp_df = @linq all_peaks_df |> where(:searched .== false) |> first

        corrent_index = temp_df.mval

        # search upstream
        while all(den.density[corrent_index] .> den.density[corrent_index+1:end])
            corrent_index += 1
            if corrent_index == length(den.density)+1
                break
            end
        end

        if corrent_index < length(den.density) .&& !(corrent_index in all_peaks_df.mval)
            if ismissing(temp_df.upper_bound)
                new_mval = findmax(den.density[corrent_index+1:end])
            else
                new_mval = findmax(den.density[temp_df.upper_bound[1]])
            end
            bounds = sort([new_mval[2]+1,corrent_index-1])
            # push!(bounds, (corrent_index+1:corrent_index+new_mval[2]-1))
            corrent_index += new_mval[2]
            next_peak = local_max(Int64(round(den.x[corrent_index])), den.density[corrent_index])
            # push!(all_peaks, next_peak)
            # tojoin = DataFrame(mval = corrent_index, location = next_peak.location, magnitude = next_peak.magnitude, upper_bound = bounds[1]+1:bounds[2]-1, lower_bound = missing, searched = false)
            # all_peaks_df = innerjoin(all_peaks_df, tojoin, on = :mval)
            if length(bounds[1]+1:bounds[2]-1) != 0
                push!(all_peaks_df, [corrent_index, next_peak.location, next_peak.magnitude, bounds[1]+1:bounds[2]-1, missing, false])
            else
                push!(all_peaks_df, [corrent_index, next_peak.location, next_peak.magnitude, missing, missing, false])
            end
        end

        corrent_index = first(temp_df.mval)

        # search downtream
        while all(den.density[corrent_index] .> den.density[1:corrent_index-1])
            corrent_index -= 1
            if corrent_index == 0
                break
            end
        end

        if corrent_index > 0 .&& !(corrent_index in all_peaks_df.mval)
            if ismissing(temp_df.lower_bound)
                new_mval = findmax(den.density[1:corrent_index-1])
            else
                new_mval = findmax(den.density[temp_df.lower_bound[1]])
            end
            bounds = sort([new_mval[2]+1,corrent_index-1])
            # push!(bounds, (new_mval[2]+1:corrent_index-1))
            corrent_index = new_mval[2]
            next_peak = local_max(Int64(round(den.x[corrent_index])), den.density[corrent_index])
            # push!(all_peaks, next_peak)
            # tojoin = DataFrame(mval = corrent_index, location = next_peak.location, magnitude = next_peak.magnitude, upper_bound = missing, lower_bound = bounds[1]+1:bounds[2]-1, searched = false)
            # all_peaks_df = innerjoin(all_peaks_df, tojoin, on = :mval)
            if length(bounds[1]+1:bounds[2]-1) != 0
                push!(all_peaks_df, [corrent_index, next_peak.location, next_peak.magnitude, missing, bounds[1]+1:bounds[2]-1, false])
            else
                push!(all_peaks_df, [corrent_index, next_peak.location, next_peak.magnitude, missing, missing, false])
            end
        end
        
        # order the peaks list by magnitude
        # decs_order = sortperm([p.magnitude for p in all_peaks], rev = true)
        # all_peaks = all_peaks[decs_order]
        all_peaks_df.searched[all_peaks_df.location .== temp_df.location] .= true
        # all_peaks_df.searched[in(all_peaks_df.location).(temp_df.location)] .= true
        sort!(all_peaks_df, :magnitude, rev = true)
    end

    return all_peaks_df
end


features_optimal = DataFrame(feature = all_features_found, occur = location_lengths)
features_length = @linq FeaturesTable |> rename(:Motif => :feature, :Length => :length) |> unique(:feature) |> select(:feature, :length)

features_optimal = innerjoin(features_optimal, features_length, on = :feature)
# FindAllPeaks(FeaturesTable, "RY-repeat", "+")

feature_peaks = DataFrame(feature = String[], location = Int64[], magnitude = Float64[])
for feature in unique(FeaturesTable.Motif)
    peaks = FindAllPeaks(FeaturesTable, feature, "+")
    # append!(feature_peaks, DataFrame(feature = fill(feature,length(peaks)), location = [p.location for p in peaks], magnitude = [p.magnitude for p in  peaks]))
    append!(feature_peaks, DataFrame(feature = fill(feature,nrow(peaks)), location = peaks.location , magnitude = peaks.magnitude))
end

features_optimal = leftjoin(features_optimal, feature_peaks, on = :feature)

using Pipe

features_optimal = @pipe features_optimal |> 
    groupby(_, [:feature, :location]) |>
    combine(_, :magnitude => maximum => :magnitude, :occur => maximum => :occur, :length => mean => :length) |> 
    transform(_, [:occur, :magnitude] => ByRow(^) => :score) |>
    sort(_, :score, rev = true)

CSV.write(join([parsed_args["output_dir"], "features_optimal.tsv"],"/"), features_optimal, delim = "\t")
# ==========================================================================================================

# assemble promoter

# features_optimal = CSV.read("features_optimal.csv", DataFrame)

rename!(regulatory_elements, names(regulatory_elements) .=> ["feature", "sequence", "length", "accesion"])

features_optimal = innerjoin(features_optimal, regulatory_elements, on = :feature)
features_optimal.sequence = LongDNA{4}.(features_optimal.sequence)

maxlength = maximum(length.(values(promotersDict)))
features_optimal = @linq features_optimal |> where(:location .> 0) |> where(:location .< maxlength)

sort!(features_optimal, :score)

template_promoter = LongDNA{4}(fill(DNA_N, maximum(features_optimal.location)+maximum(features_optimal.length)))

assembly = copy(template_promoter)

for feature in eachrow(features_optimal)
    assembly[feature.location:feature.location+feature.length-1] = feature.sequence
end


not_ambiguous = alphabet(DNA)[findall(isambiguous.(alphabet(DNA)) .== false)]
not_ambiguous = collect(not_ambiguous[findall(not_ambiguous .!= DNA_Gap)])

for (i,nucl) in enumerate(assembly)
    if isambiguous(nucl)
        # compatible nucleotides
        compatibles = not_ambiguous[iscompatible.(nucl, not_ambiguous)]
        assembly[i] = rand(compatibles)
    end
end

SynthProm = FASTA.Record("SynthProm", assembly)
promoters = FASTA.Writer(open(parsed_args["promoter"], "a"))
write(promoters, SynthProm)
close(promoters)

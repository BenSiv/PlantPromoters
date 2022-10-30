"""
This scripts contains functions to generate a minimal promoter sequence based on gegulatory elements
"""

# using JLD

# load("Features.jld")["Features"]

# FeatureLoc = Dict()
# for key in keys(Features)
#     FeatureRange = vcat(values(Features[key])...)
#     push!(FeatureLoc, "$key" => vcat([[f;] for f in FeatureRange]...))
# end 

# FeatureLocUnique = Dict()
# for key in keys(FeatureLoc)
#     push!(FeatureLocUnique, "$key" => sort(unique(FeatureLoc[key])))
# end 


# using DataStructures

# FeatureLocCount = Dict()
# for key in keys(FeatureLoc)
#     push!(FeatureLocCount, "$key" => DataStructures.counter(FeatureLoc[key]))
# end 


function Dilution(NumSerie, DilutionFactor, PromLength)
    if minimum(NumSerie) - DilutionFactor < 1
        Start = 1
    else
        Start = minimum(NumSerie) - DilutionFactor
    end

    if maximum(NumSerie) + DilutionFactor > PromLength
        End = PromLength
    else
        End = maximum(NumSerie) + DilutionFactor
    end
    
    AllRange = [Start:End;]

    Alignment = in(NumSerie).(AllRange)

    cnt = 0
    for (i,x) in enumerate(Alignment)

        if x == 0 && cnt == DilutionFactor
            Alignment[i-DilutionFactor:i] .= 1
            cnt += 1
        elseif x == 0 && i == length(Alignment) && cnt <= DilutionFactor
            Alignment[i-cnt:i] .= 1
        elseif x == 0
            cnt += 1
        elseif cnt <= DilutionFactor-1
            cnt = 0
            Alignment[i-cnt:i] .= 1
        else 
            cnt = 0
            Alignment[i-DilutionFactor:i] .= 1
        end

    end

    return AllRange[Alignment]
end


function GetLocations(Vec_of_feature_locations)
    FeatureRange = vcat(values(Vec_of_feature_locations)...)
    FeatureRangeConcat = vcat([[f;] for f in FeatureRange]...)
    FeatureLocUnique = sort(unique(FeatureRangeConcat))
    return FeatureLocUnique
end

function GetLocSeq(sequence, locations)
    seq = LongDNA{4}([sequence[i] for i in locations])
    
    return seq
end

# using Plots, StatsPlots

# density(FeatureLoc["At4g28520"])
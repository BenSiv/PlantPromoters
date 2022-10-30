"""
A module to construct a GenBank file
"""

using Dates 
using BioSequences

mutable struct Feature
    name::String
    location::UnitRange{Int}
    complement::Bool
end

mutable struct Record
    name::String
    sequence::LongDNA{4}
    features::Vector{Feature}
end

function addfeature!(record::Record, feature::Feature)
    push!(record.features, feature)
end

"""

SomeRecord = Record("PL999", dna"ATATATATATAT", Vector{Feature}[])
SomeFeature = Feature("TATA-Box", 2:5, false)

addfeature!(SomeRecord, SomeFeature)
"""


function writeGB(record::Record)

    ORIGIN = "1 "
    for n in 1:length(record.sequence)
        if n%60 == 0
            ORIGIN *= string(record.sequence[n])*"\n$(n+1) "
        elseif n%10 == 0
            ORIGIN *= string(record.sequence[n])*" "
        else
            ORIGIN *= string(record.sequence[n])
        end
    end

    FEATURES = ""
    for feature in record.features
        if feature.complement
            FEATURES *= 
            """
                misc_feature    complement($(first(feature.location))..$(last(feature.location)))
                                /label="$(feature.name)"
            """
        else
            FEATURES *= 
            """
                misc_feature    $(first(feature.location))..$(last(feature.location))
                                /label="$(feature.name)"
            """
        end
    end

    GenBank = 
    """
    LOCUS       $(record.name)   $(length(record.sequence)) bp ds-DNA     linear       $(Dates.today())
    DEFINITION  .
    FEATURES            Location/Qualifiers
    $FEATURES
    ORIGIN
        $ORIGIN
    //

    """

    return GenBank
end

function write_gb_records(promoters_dict, features_table)
    GBRecords = []
    for key in keys(promoters_dict)
        push!(GBRecords, Record(key, promoters_dict[key], Vector{Feature}[]))
        promoterDF = @linq features_table |> where(:Promoter .== key)
        for feature in unique(promoterDF.Motif)
            featureDF = @linq promoterDF |> where(:Motif .== feature)
            for occurence in eachrow(featureDF)
                if occurence.Direction == "+"
                    SomeFeature = Feature(feature, occurence.Occurences, false)
                    addfeature!(last(GBRecords), SomeFeature)
                else
                    SomeFeature = Feature(feature, occurence.Occurences, true)
                    addfeature!(last(GBRecords), SomeFeature)
                end
            end
        end
    end
    
    return GBRecords
end

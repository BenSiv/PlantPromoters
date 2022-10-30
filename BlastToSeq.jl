"""
This is a script that takes as input tsv as outputted from blastn (outfmt 6) and returns eather the sequence in the specify location or possible promoter regon.
"""

using Pkg
Pkg.activate("/home/bensiv/Dropbox/202207-RemilkProject/")

using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--input", "-i"
            help = "blastn results (outfmt 6) as a tsv file"
            required = true
        "--genome", "-g"
            help = "genome file for sequence extraction"
            required = true
        "--output", "-o"
            help = "outpput name as a fasta file"
            default = "output.fasta"
        "--promoter", "-p"
            help = "return 1000 bp upstream to the sequence + 200 bp into the sequence"
            action = :store_true
        "--upstream", "-u"
            help = "defines the number of bp to take upstream"
            arg_type = Int
            default = 1000
        "--intoseq", "-d"
            help = "defines the number of bp to take into the sequence"
            arg_type = Int
            default = 200
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

using CSV, DataFrames, FASTX, BioSequences

Blast = CSV.read(parsed_args["input"], DataFrame, delim = "\t")

if parsed_args["promoter"]
    Promoters = []
    for blast in eachrow(Blast)

        Genome = open(FASTA.Reader, parsed_args["genome"])

        for record in Genome
            if FASTA.identifier(record) == blast.Reference
                if blast.Strand == "+"
                    if blast.Ref_Start-parsed_args["intoseq"] < 1
                        push!(Promoters, FASTA.sequence(record)[1:blast.Ref_Start+parsed_args["upstream"]])    
                    else
                        push!(Promoters, FASTA.sequence(record)[blast.Ref_Start-parsed_args["intoseq"]:blast.Ref_Start+parsed_args["upstream"]]) 
                    end
                else
                    if blast.Ref_End+parsed_args["intoseq"] > length(FASTA.sequence(record))
                        push!(Promoters, FASTA.sequence(record)[blast.Ref_End-parsed_args["upstream"]:end]) 
                    else
                        push!(Promoters, FASTA.sequence(record)[blast.Ref_End-parsed_args["upstream"]:blast.Ref_End+parsed_args["intoseq"]]) 
                    end
                end
            end
        end
    end

    w = open(FASTA.Writer, parsed_args["output"])
    for (seq, i) in zip(Promoters, 1:length(Promoters))
        rec = FASTA.Record("promoter_$i", seq)
        write(w, rec)
    end
    close(w)

else

    Sequences = []
    for blast in eachrow(Blast)

        Genome = open(FASTA.Reader, parsed_args["genome"])

        for record in Genome
            if FASTA.identifier(record) == blast.Reference
                if blast.Strand == "+"
                    push!(Sequences, FASTA.sequence(record)[blast.Ref_Start:blast.Ref_End])
                else
                    push!(Sequences, FASTA.sequence(record)[blast.Ref_End:blast.Ref_Start]) 
                end
            end
        end
    end

    w = open(FASTA.Writer, parsed_args["output"])
    for (seq, i) in zip(Sequences, 1:length(Sequences))
        rec = FASTA.Record("sequence_$i", seq)
        write(w, rec)
    end
    close(w)

end

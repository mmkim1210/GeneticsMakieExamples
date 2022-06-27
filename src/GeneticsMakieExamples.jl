using GeneticsMakie, CairoMakie, CSV, DataFrames, SnpArrays, Arrow

CairoMakie.activate!(type = "png")
set_theme!(font = "Arial")

@info "Loading GENCODE v39 annotation for chromosome 15"
if !isfile(joinpath(@__DIR__, "../data/gencode.v19.annotation.chr15.parsed.gtf.arrow"))
    @time gencode = CSV.read(joinpath(@__DIR__, "../data/gencode.v39lift37.annotation.chr15.gtf.gz"), DataFrame,
        header = ["seqnames", "source", "feature", "start", "end", "score", "strand", "phase", "info"],
        delim = "\t", skipto = 6)
    GeneticsMakie.parsegtf!(gencode)
    select!(gencode, :seqnames, :feature, :start, :end, :strand, :gene_id, :gene_name, :gene_type, :transcript_id)
    @time Arrow.write(joinpath(@__DIR__, "../data/gencode.v39lift37.annotation.chr15.parsed.gtf.arrow"), gencode)
end
@time gencode = Arrow.Table(joinpath(@__DIR__, "../data/gencode.v39lift37.annotation.chr15.parsed.gtf.arrow"))|> DataFrame
@assert (117_418, 9) == size(gencode)

@info "Loading 1000 Genomes reference panel for chromosome 15"
kgp = SnpData(joinpath(@__DIR__, "../data/kgp.chr15"))
@assert (503, 200_311) == size(kgp)

@info "Loading GWAS results for chromosome 15"
phenos = ["scz", "bd", "asd", "adhd", "neuroticism", "alz", "menopause", "weight", "height"]
gwas = []
for p in phenos
    @time if !isfile(joinpath(@__DIR__, "../data/", p * ".gwas.arrow"))
        sumstat = CSV.read(joinpath(@__DIR__, "../data/$p.chr15.txt.gz"), DataFrame, comment = "##", missingstring = ["NA"])
        GeneticsMakie.mungesumstats!(sumstat)
        Arrow.write(joinpath(@__DIR__, "../data", p * ".gwas.arrow"), sumstat)
    end
    push!(gwas, DataFrame(Arrow.Table(joinpath(@__DIR__, "../data", p * ".gwas.arrow"))))
end

function subsetgwas(gwas, chr::AbstractString, range1::Real, range2::Real)
    gwas_subset = Vector{DataFrame}(undef, length(gwas))
    for i in 1:length(gwas)
        df = gwas[i]
        gwas_subset[i] = df[findall((df.CHR .== chr) .& (df.BP .>= range1) .& (df.BP .<= range2)), :]
    end
    gwas_subset
end

issig(P::AbstractVector; p = 5e-8) = any(P .< p)
issig(df::DataFrame; p = 5e-8) = issig(df.P; p = p)

ispath(joinpath(@__DIR__, "../figs")) || mkpath(joinpath(@__DIR__, "../figs"))
function locuszoom(genes)
    for gene in genes
        @info "Working on $gene gene"
        window = 1e6
        chr, start, stop = GeneticsMakie.findgene(gene, gencode)
        range1, range2 = start - window, stop + window
        @info "Subsetting GWAS results"
        @time gwas_subset = subsetgwas(gwas, chr, range1, range2)
        titles = [GeneticsMakie.gwas[p].title for p in phenos]
        @info "Plotting LocusZoom"
        n = length(titles)
        f = Figure(resolution = (306, 1500))
        axs = [Axis(f[i, 1]) for i in 1:(n + 1)]
        for i in 1:n
            if issig(gwas_subset[i])
                GeneticsMakie.plotlocus!(axs[i], chr, range1, range2, gwas_subset[i]; ld = kgp)
                if gwas_subset[i].BP[argmin(gwas_subset[i].P)] < (range1 + range2) / 2
                    Label(f[i, 1, Top()], "$(titles[i])", textsize = 6, halign = :right, padding = (0, 7.5, -5, 0))
                else
                    Label(f[i, 1, Top()], "$(titles[i])", textsize = 6, halign = :left, padding = (7.5, 0, -5, 0))
                end    
            else
                GeneticsMakie.plotlocus!(axs[i], chr, range1, range2, gwas_subset[i])
                Label(f[i, 1, Top()], "$(titles[i])", textsize = 6, halign = :left, padding = (7.5, 0, -5, 0))
            end
            rowsize!(f.layout, i, 30)
        end
        rs = GeneticsMakie.plotgenes!(axs[n + 1], chr, range1, range2, gencode; height = 0.1)
        rowsize!(f.layout, n + 1, rs)
        GeneticsMakie.labelgenome(f[n + 1, 1, Bottom()], chr, range1, range2)
        Colorbar(f[1:n, 2], limits = (0, 1), ticks = 0:1:1, height = 20,
            colormap = (:gray60, :red2), label = "LD", ticksize = 0, tickwidth = 0,
            tickalign = 0, ticklabelsize = 6, flip_vertical_label = true,
            labelsize = 6, width = 5, spinewidth = 0.5)
        Label(f[1:n, 0], text = "-log[p]", textsize = 6, rotation = pi / 2)
        for i in 1:(n + 1)
            vlines!(axs[i], start, color = (:gold, 0.5), linewidth = 0.5)
            vlines!(axs[i], stop, color = (:gold, 0.5), linewidth = 0.5)
        end
        for i in 1:n
            lines!(axs[i], [range1, range2], fill(-log(10, 5e-8), 2), color = (:purple, 0.5), linewidth = 0.5)
        end
        rowgap!(f.layout, 5)
        colgap!(f.layout, 5)
        resize_to_layout!(f)
        save(joinpath(@__DIR__, "../figs/$(gene)-locuszoom.png"), f, px_per_unit = 4)
    end
end

genes = ["CHRNA5"]
@time locuszoom(genes)

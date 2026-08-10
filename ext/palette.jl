const INERT_COLOR = colorant"#E7E7E7"
const DEFAULT_COLORS = [
    [colorant"#B2D0EA", colorant"#54A3E4", colorant"#1A78C6", colorant"#00549A"],
    [colorant"#F39F9D", colorant"#ED7E7C", colorant"#E75451", colorant"#BE2E2C"],
    [colorant"#FEE3B5", colorant"#FFC563", colorant"#FFB12F", colorant"#FFA000"],
    [colorant"#DBB9E4", colorant"#D99DE8", colorant"#B571C4", colorant"#8E4E9C"],
    [colorant"#AAF5FF", colorant"#86E4F1", colorant"#43CFE2", colorant"#20BCD2"],
    [colorant"#BCF3E1", colorant"#3DD4A4", colorant"#19B684", colorant"#008D60"],
    [colorant"#FFD4B6", colorant"#F7AF7D", colorant"#FF8835", colorant"#ED6304"],
    [colorant"#C5C9FE", colorant"#A9AFFF", colorant"#7C85FF", colorant"#4854FD"],
    [colorant"#F4FFC4", colorant"#E0F290", colorant"#B9D63E", colorant"#859F13"],
    [colorant"#FFB2C5", colorant"#FF8EA9", colorant"#FF5C83", colorant"#E32654"],
    [colorant"#95FABA", colorant"#4DD980", colorant"#21AB53", colorant"#008832"]
]

function species_palette(spcs_index::Int, n::Int)
    base = DEFAULT_COLORS[mod1(spcs_index, length(DEFAULT_COLORS))]
    return cgrad(base, n; categorical=true)
end

"""
    _resolve_colors(spcs, species_index, sys, site_color)

Work out which palette a species is drawn in and what color each of its sites gets.

`species_index` selects the palette; when it is `nothing` it is looked up in `sys`, falling
back to `1`. `site_color` is an optional callback `(species_index, site_index) -> color`
that overrides both the palette and the greying of inert sites.

Returns `(species_index, palette, sitecolors)`.
"""
function _resolve_colors(spcs, species_index, sys, site_color)
    n = nsites(spcs)
    if isnothing(species_index)
        species_index = isnothing(sys) ? 1 : something(findfirst(==(spcs), species(sys)), 1)
    end
    pal = species_palette(species_index, n)

    if isnothing(site_color)
        site_color = if isnothing(sys)
            (_, i) -> pal[i]
        else
            (_, i) -> isinert(sys, (species_index, i)) ? INERT_COLOR : pal[i]
        end
    end
    return species_index, pal, [site_color(species_index, i) for i in 1:n]
end

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

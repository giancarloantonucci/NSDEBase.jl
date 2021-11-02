@userplot PHASEPLOT
@recipe function f(h::PHASEPLOT)
    @series begin
        plotkind --> :phaseplot
        h.args[1]
    end
end
